#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stats_core.py — 独立イベント集約とクラスタ頑健な区間推定
=========================================================

★このモジュールが解決する問題★

volatility_study.py は 1518 件の検知を「1518個の独立した観測」として扱い、
Wilcoxon 検定とブートストラップにかけている。これは誤りである可能性が高い。

理由は2つ:

  ① フォワード窓の重複
     4h先を評価するのに検知が2時間おきに出れば、隣り合う検知は
     同じ値動きを2回数えている。24h先ならもっとひどい。

  ② 銘柄間の相関
     BTC が急騰すれば ETH も SOL も DOGE も動く。同じ時刻に4銘柄で
     検知が出たとき、それは4つの独立な証拠ではなく【1つの相場イベント】である。

Wilcoxon 検定もブートストラップも観測の独立を仮定している。
独立でないものを独立として数えると、実効サンプル数が水増しされ、
**信頼区間は不当に狭くなり、p値は不当に小さくなる**。

つまり「1.38倍 CI 1.30–1.47 p<0.001」という数字は、
効果の【大きさ】は正しくても【確からしさ】を過大に主張している可能性がある。
公開前に必ず潰しておく必要がある。

★対処★

  1. 時間的に重なる検知を1つの【イベントクラスタ】にまとめる（銘柄をまたいで）
  2. クラスタ単位に集約してから検定する（実効サンプル数が正直な値になる）
  3. 区間推定は移動ブロック・ブートストラップで行う
     （クラスタ同士も時間的に相関するため、単純リサンプリングでは不十分）

依存ライブラリなし。既存の統計定義は volatility_study.py から借りて一本化する
（同じ検定を2箇所に書くと、片方だけ直して片方が腐る。実際に track_record.py で
 その事故が起きた）。
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from volatility_study import (TestResult, median, percentile,
                              two_sided_p_from_z, wilcoxon_signed_rank)

__all__ = ["build_clusters", "build_clusters_keyed", "cluster_reduce",
           "ClusterSummary", "moving_block_bootstrap_median",
           "analyze_clustered", "required_sample_size", "median",
           "wilcoxon_signed_rank", "ci_contains", "ci_excludes_above",
           "effect_verdict", "control_verdict"]

TF_MS = 3_600_000          # 1時間足のミリ秒
DEFAULT_ITERS = 2000


# ============================================================================
# ⓪ 判定ヘルパー — 判定ロジックはここ1箇所にしか書かない
# ============================================================================
# ★なぜ共通化するか★
# 陰性対照の判定を「1.0を超えていないか」だけで書いた結果、
# 1.0を有意に【下回る】異常（0.90倍 p<0.001）を「OK」と誤表示した。
# 同じ判定を複数箇所に書けば、必ずどれかが片側のまま取り残される。
# ゲートを追加するときは、必ずこの2つのどちらかを呼ぶこと。

def ci_contains(ci: tuple[float, float], value: float = 1.0) -> bool:
    """信頼区間が value を含むか。「差があるとは言えない」の判定。"""
    lo, hi = ci
    return lo <= value <= hi


def ci_excludes_above(ci: tuple[float, float], value: float = 1.0) -> bool:
    """信頼区間が完全に value より上にあるか。「効果あり」の判定。"""
    return ci[0] > value


def effect_verdict(ci: tuple[float, float], n: int,
                   min_n: int = 30, null: float = 1.0) -> str:
    """
    効果の有無の判定。'effect' / 'none' / 'insufficient'。
    片側だけ見て「効果あり」と言わないための唯一の入口。
    """
    if n < min_n:
        return "insufficient"
    return "effect" if ci_excludes_above(ci, null) else "none"


def control_verdict(ci: tuple[float, float], n: int,
                    min_n: int = 30, null: float = 1.0) -> str:
    """
    陰性対照（プラセボ）の判定。'clean' / 'biased' / 'insufficient'。

    ★両側で見る★
    上振れも下振れも異常。CIが null を含むことが唯一の合格条件。
    「1.0を超えていないから正常」は誤り。
    """
    if n < min_n:
        return "insufficient"
    return "clean" if ci_contains(ci, null) else "biased"


# ============================================================================
# ① イベントクラスタの構築
# ============================================================================
def build_clusters(timestamps_ms: list[int], horizon_bars: int,
                   tf_ms: int = TF_MS) -> list[list[int]]:
    """
    時間的に重なる検知を1つのイベントにまとめる。銘柄はまたぐ。

    引数の timestamps_ms は【時系列順にソート済み】であること。
    戻り値は、元のリストのインデックスをクラスタごとにまとめたもの。

    ★まとめ方（貪欲ブロッキング）★
    先頭の検知を起点とし、起点から horizon_bars 本ぶんの時間内に入る検知を
    すべて同じクラスタに入れる。範囲外に出た最初の検知が次の起点になる。

    単一連結（隣接する検知を次々に連結する）を使わない理由:
    検知が密な銘柄では連鎖が止まらず、全期間が1クラスタに潰れてしまう。
    貪欲ブロッキングなら、クラスタの窓が互いにほぼ重ならないことが保証され、
    かつクラスタ数が過度に減らない。区間推定としては保守的すぎず甘すぎない。

    ★銘柄をまたぐ理由★
    同時刻の BTC と ETH の検知は、独立な2つの証拠ではなく1つの相場イベント。
    銘柄ごとにクラスタを作ると、この相関をまったく潰せない。
    """
    if not timestamps_ms:
        return []
    span = max(1, horizon_bars) * tf_ms
    clusters: list[list[int]] = []
    start_ts = timestamps_ms[0]
    current = [0]
    for i in range(1, len(timestamps_ms)):
        ts = timestamps_ms[i]
        if ts - start_ts < span:
            current.append(i)
        else:
            clusters.append(current)
            current = [i]
            start_ts = ts
    clusters.append(current)
    return clusters


def build_clusters_keyed(timestamps_ms: list[int], keys: list,
                         window_bars: int = 1,
                         tf_ms: int = TF_MS) -> list[list[int]]:
    """
    事前登録 v2.0 §2 の独立イベント定義。

        ① 銘柄×方向 で window_bars の窓の重複検知を1件に集約
        ② 方向別に、同じ窓で銘柄横断クラスタ化

    keys には各検知の「方向」を渡す（BULLISH / BEARISH）。
    同一方向・同一時間窓の検知は、銘柄が違っても1イベントとして数える。

    ★方向で分ける影響について★
    同時刻に BTC が BULLISH、ETH が BEARISH で検知された場合、この定義では
    2イベントと数える。より保守的な立場（同時刻はすべて1イベント）もあり得るが、
    事前登録で採用した定義はこちらなので、凍結後は変更しない。
    変えたくなったら v3.0 を作ってカウントをやり直すこと。
    """
    if not timestamps_ms:
        return []
    span = max(1, window_bars) * tf_ms
    by_key: dict[object, list[int]] = {}
    for i, k in enumerate(keys):
        by_key.setdefault(k, []).append(i)

    clusters: list[list[int]] = []
    for _k, idxs in by_key.items():
        idxs.sort(key=lambda i: timestamps_ms[i])
        start = timestamps_ms[idxs[0]]
        cur = [idxs[0]]
        for i in idxs[1:]:
            if timestamps_ms[i] - start < span:
                cur.append(i)
            else:
                clusters.append(cur)
                cur = [i]
                start = timestamps_ms[i]
        clusters.append(cur)
    clusters.sort(key=lambda c: min(timestamps_ms[i] for i in c))
    return clusters


def cluster_reduce(values: list[float], clusters: list[list[int]]) -> list[float]:
    """
    クラスタ内の値を中央値で1つにまとめる。

    平均ではなく中央値を使うのは、log比の分布が歪んでおり、
    クラスタ内に1件だけ極端な値があるとクラスタ全体を代表してしまうため。
    """
    return [median([values[i] for i in idx]) for idx in clusters]


# ============================================================================
# ② 移動ブロック・ブートストラップ
# ============================================================================
def moving_block_bootstrap_median(values: list[float], rng: random.Random,
                                  block: int | None = None,
                                  iters: int = DEFAULT_ITERS
                                  ) -> tuple[float, float]:
    """
    時系列順に並んだ値の中央値の95%信頼区間（対数スケールのまま返す）。

    ★なぜ単純ブートストラップでは足りないか★
    クラスタにまとめても、隣り合うクラスタは同じ相場つきの中にある。
    ボラティリティ・レジームは数週間続くので、クラスタ同士も相関する。
    1点ずつ独立に復元抽出すると、この相関を壊してしまい区間がまた狭くなる。

    連続した block 個をひとかたまりとして抜き出して並べ直せば、
    かたまりの中の相関構造は保存される。
    ブロック長の既定値は n^(1/3)（標準的な選び方）。
    """
    n = len(values)
    if n < 10:
        return (0.0, 0.0)
    if block is None:
        block = max(2, int(round(n ** (1 / 3))))
    block = min(block, n)
    n_blocks = math.ceil(n / block)

    meds = []
    for _ in range(iters):
        sample: list[float] = []
        for _b in range(n_blocks):
            s = rng.randrange(0, n - block + 1)
            sample.extend(values[s:s + block])
        meds.append(median(sample[:n]))
    meds.sort()
    return (percentile(meds, 0.025), percentile(meds, 0.975))


# ============================================================================
# ③ クラスタ頑健な検定のまとめ
# ============================================================================
@dataclass
class ClusterSummary:
    """生の観測とクラスタ集約後の結果を並べて持つ。差そのものが報告対象。"""
    label: str = ""
    n_raw: int = 0
    n_clusters: int = 0
    ratio_raw: float = 0.0
    ci_raw: tuple[float, float] = (0.0, 0.0)
    p_raw: float = 1.0
    ratio: float = 0.0                       # クラスタ集約後
    ci: tuple[float, float] = (0.0, 0.0)
    p: float = 1.0
    block: int = 0

    @property
    def inflation(self) -> float:
        """水増し率。生の件数がクラスタ数の何倍あったか。"""
        return self.n_raw / self.n_clusters if self.n_clusters else 0.0

    @property
    def ci_widening(self) -> float:
        """信頼区間が何倍に広がったか。1.0 に近ければ元の主張は無事。"""
        w0 = self.ci_raw[1] - self.ci_raw[0]
        w1 = self.ci[1] - self.ci[0]
        return w1 / w0 if w0 > 0 else 0.0

    @property
    def survives(self) -> bool:
        """クラスタ補正後もなお「1.0倍を超える」と言えるか。"""
        return (self.n_clusters >= 30 and self.p < 0.05 and self.ci[0] > 1.0)

    @property
    def tag(self) -> str:
        if self.n_clusters < 30:
            return "判定不能"
        if not self.survives:
            return "補正後は有意でない"
        return "補正後も有意"


def analyze_clustered(log_ratios: list[float], timestamps_ms: list[int],
                      horizon_bars: int, rng: random.Random,
                      label: str = "", iters: int = DEFAULT_ITERS,
                      tf_ms: int = TF_MS) -> ClusterSummary:
    """
    log比の列を、生のまま／クラスタ集約後 の両方で検定して並べて返す。

    log_ratios と timestamps_ms は同じ順序（時系列順）であること。
    倍率は exp して返すので、呼び出し側はそのまま「◯◯倍」と読める。
    """
    s = ClusterSummary(label=label, n_raw=len(log_ratios))
    if not log_ratios:
        return s

    # --- 生の観測（volatility_study.py と同じ計算。比較のために残す）---
    s.ratio_raw = math.exp(median(log_ratios))
    s.p_raw = wilcoxon_signed_rank(log_ratios).p
    n = len(log_ratios)
    meds = []
    for _ in range(iters):
        meds.append(median([log_ratios[rng.randrange(n)] for _ in range(n)]))
    meds.sort()
    s.ci_raw = (math.exp(percentile(meds, 0.025)),
                math.exp(percentile(meds, 0.975)))

    # --- クラスタ集約 ---
    clusters = build_clusters(timestamps_ms, horizon_bars, tf_ms)
    reduced = cluster_reduce(log_ratios, clusters)
    s.n_clusters = len(reduced)
    s.block = max(2, int(round(len(reduced) ** (1 / 3))))
    s.ratio = math.exp(median(reduced))
    s.p = wilcoxon_signed_rank(reduced).p
    lo, hi = moving_block_bootstrap_median(reduced, rng, s.block, iters)
    s.ci = (math.exp(lo), math.exp(hi))
    return s


# ============================================================================
# ④ 検出力：何件あれば主張できるのか
# ============================================================================
def required_sample_size(p1: float, p0: float = 0.5, alpha: float = 0.05,
                         power: float = 0.80) -> int:
    """
    両側検定で、勝率 p1 を検出するのに必要な【独立】サンプル数。

    「100件で p<0.05 が出た」は、検出力が足りていないか独立性が崩れているか
    のどちらかである、と判断するための基準値。
    正規近似（arcsine 変換なしの素直な二項近似）。
    """
    if not (0 < p1 < 1) or p1 == p0:
        return 0
    z_a = 1.959963984540054 if abs(alpha - 0.05) < 1e-9 else _z_from_p(alpha / 2)
    z_b = 0.8416212335729143 if abs(power - 0.80) < 1e-9 else _z_from_p(1 - power)
    num = (z_a * math.sqrt(p0 * (1 - p0)) + z_b * math.sqrt(p1 * (1 - p1))) ** 2
    return math.ceil(num / (p1 - p0) ** 2)


def _z_from_p(p: float) -> float:
    """両側p値に対応する標準正規の分位点（二分法。scipy 不要）。"""
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if two_sided_p_from_z(mid) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ============================================================================
# 自己検査
# ============================================================================
if __name__ == "__main__":
    rng = random.Random(20260811)

    print("── build_clusters の挙動確認 ─────────────────────────")
    ts = [0, TF_MS, 2 * TF_MS, 10 * TF_MS, 10 * TF_MS, 11 * TF_MS]
    for h in (1, 4, 24):
        cl = build_clusters(ts, h)
        print(f"  horizon={h:>2}本  検知6件 → クラスタ {len(cl)} 個  {cl}")

    print()
    print("── 独立なデータでは補正してもほぼ変わらないこと ─────")
    n = 800
    indep = [rng.gauss(0.25, 1.0) for _ in range(n)]
    ts_far = [i * 100 * TF_MS for i in range(n)]      # 十分に離れている
    s = analyze_clustered(indep, ts_far, 1, rng, "独立", iters=400)
    print(f"  生 n={s.n_raw} → クラスタ {s.n_clusters}（水増し {s.inflation:.2f}倍）")
    print(f"  CI 幅の変化: {s.ci_widening:.2f}倍  ← 1.0付近なら正常")

    print()
    print("── 重複したデータでは区間が広がること ───────────────")
    ts_near = [i * TF_MS for i in range(n)]           # 1本ずつ隣接
    s2 = analyze_clustered(indep, ts_near, 24, rng, "重複", iters=400)
    print(f"  生 n={s2.n_raw} → クラスタ {s2.n_clusters}（水増し {s2.inflation:.2f}倍）")
    print(f"  CI 幅の変化: {s2.ci_widening:.2f}倍  ← 1.0より大きければ正常")

    print()
    print("── 必要サンプル数（両側α=0.05, 検出力80%）───────────")
    for p in (0.55, 0.575, 0.60):
        print(f"  勝率 {p*100:.1f}% を検出するのに必要な独立件数: "
              f"{required_sample_size(p)} 件")
