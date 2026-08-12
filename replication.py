#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replication.py — 事前登録 v2.1 の確認的再現プロトコル実装
==========================================================

★このファイルは pre_registration_v2.md §2 の定義を機械的に実装したものです★
凍結後にここを編集してはいけません。編集が必要になったら v3.0 を作り、
カウントをリセットして両方を残します。

──────────────────────────────────────────────────────────────
探索的 vs 確認的
──────────────────────────────────────────────────────────────
13ヶ月・1518件で得た 1.33倍 / 1.50倍 は【探索的】です。
3つの horizon と2つのベースラインを見た後に「1h・HAR-RV比」を選んでいるので、
どれだけ厳密に計算しても事後選択の産物である可能性を排除できません。

本モジュールが出す確認的な数字は、**凍結時刻以降に発生した検知だけ**を使います。
過去データからホールドアウトを切り出すことは禁止です（既に見てしまっているため）。

    python replication.py --exploratory   # 過去データ。探索的ラベル付きで参考出力
    python replication.py                 # 確認的。FREEZE.json 以降のみ
    python replication.py --selftest      # 自己検証（CI で必ず回すこと）

──────────────────────────────────────────────────────────────
主要エンドポイント（§2 の完全定義）
──────────────────────────────────────────────────────────────
  horizon        1時間（検知足の確定時刻を起点。未確定足は使わない）
  実現ボラ       次足から1時間ぶんの対数リターン二乗和の平方根
  第1段          log RV_fwd = a + b1 logRV(1h) + b2 logRV(24h) + b3 logRV(168h)
                 → 対照足のみで OLS 学習し、全足の HAR_pred を得る
  第2段（主要）  log RV_fwd = a + b·log(HAR_pred) + c·D_signal
  主要指標       exp(c)
  標準誤差       日次クラスタ頑健（Liang-Zeger）と
                 移動ブロック・ブートストラップの【広い方】を採用
  有意水準       α = 0.05（両側）
  判定           exp(c) の 95%CI 下限 > 1.0 なら再現成功。CIが1.0を含めば失敗
  独立イベント   1時間窓・銘柄横断・方向は区別しない（v2.1）
  判定件数       TARGET_EVENTS に到達した日に1回だけ。中間判定はしない
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import stats_core
from volatility_study import load_bars, load_signals, median, percentile

# ---- §2 で凍結する定数。実行時オプションで変えられないようにする ----
HORIZON_BARS = 1                      # 1時間
HAR_LAGS = (1, 24, 168)               # log RV(1h) / (24h) / (168h)
WARMUP = max(HAR_LAGS)                # 168本
EVENT_WINDOW_BARS = 1                 # 独立イベントの窓
TARGET_EVENTS = 267                   # 判定件数（power_calc.py で確定）
# ★267 の根拠★
# 点推定 1.312倍 ではなく【CI下限 1.215倍 を真値と仮定】して検出力80%となる件数。
# 探索段階の推定値は 3 horizon × 2ベースライン × 複数仕様を見た後に選ばれており、
# 系統的に上振れする（勝者の呪い）。点推定で設計すると構造的に under-power になり、
# 失敗したとき「効果が無い」のか「検出力不足」のか区別がつかなくなる。
# さらに §7 により再現失敗は【LP掲載・課金者通知・返金提示】を自動発動させるため、
# 検出力不足は統計上の不備ではなく自社プロダクトを誤って停止させる事業リスクである。
# 150件では CI下限が真値でも偽陰性率 44.3%（コイン投げ）だった。
ALPHA_Z = 1.959963984540054           # 両側 95%
BOOTSTRAP_ITERS = 2000
BLOCK_DAYS = 3                        # 移動ブロックの既定ブロック長（日）
MS_DAY = 86_400_000

FREEZE_FILE = Path("FREEZE.json")


# ============================================================================
# 線形代数（numpy 不要）
# ============================================================================
def solve(a: list[list[float]], b: list[float]) -> list[float] | None:
    """部分ピボット付きガウス消去。k は 3〜4 程度しか来ない。"""
    k = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            return None
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        for j in range(col, k + 1):
            m[col][j] /= pv
        for r in range(k):
            if r != col and m[r][col]:
                f = m[r][col]
                for j in range(col, k + 1):
                    m[r][j] -= f * m[col][j]
    return [m[i][k] for i in range(k)]


def inverse(a: list[list[float]]) -> list[list[float]] | None:
    k = len(a)
    cols = []
    for j in range(k):
        e = [1.0 if i == j else 0.0 for i in range(k)]
        c = solve(a, e)
        if c is None:
            return None
        cols.append(c)
    return [[cols[j][i] for j in range(k)] for i in range(k)]


def ols_fit(x: list[list[float]], y: list[float]) -> list[float] | None:
    """切片は呼び出し側が x に含めること。"""
    if not x or len(x) != len(y):
        return None
    k = len(x[0])
    if len(y) <= k:
        return None
    xtx = [[0.0] * k for _ in range(k)]
    xty = [0.0] * k
    for row, yi in zip(x, y):
        for i in range(k):
            ri = row[i]
            for j in range(k):
                xtx[i][j] += ri * row[j]
            xty[i] += ri * yi
    return solve(xtx, xty)


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


# ============================================================================
# 特徴量
# ============================================================================
@dataclass
class Row:
    ts: int
    symbol: str
    day: int                  # 日次クラスタのキー
    log_fwd: float
    har_x: list[float]        # [1, logRV1, logRV24, logRV168]
    is_signal: bool = False
    direction: str = ""


def build_rows(symbol: str, bars: list[list],
               sig_at: dict[int, str]) -> list[Row]:
    """
    ★先読みバイアスの排除★
    説明変数は足 i までの終値のみ。目的変数は i+1（1時間先）。
    足 i 自身の高値・安値・出来高は使わない。
    """
    n = len(bars)
    closes = [float(b[4]) for b in bars]
    s = [0.0] * n
    for i in range(1, n):
        if closes[i] > 0 and closes[i - 1] > 0:
            s[i] = s[i - 1] + math.log(closes[i] / closes[i - 1]) ** 2
        else:
            s[i] = s[i - 1]

    def rv(a: int, b: int) -> float:
        if a < 0 or b >= n or b <= a:
            return 0.0
        return math.sqrt(max(0.0, s[b] - s[a])) * 100.0

    out: list[Row] = []
    for i in range(WARMUP, n - HORIZON_BARS):
        fwd = rv(i, i + HORIZON_BARS)
        lags = [rv(i - L, i) for L in HAR_LAGS]
        if fwd <= 0 or min(lags) <= 0:
            continue
        ts = int(bars[i][0])
        out.append(Row(
            ts=ts, symbol=symbol, day=ts // MS_DAY,
            log_fwd=math.log(fwd),
            har_x=[1.0] + [math.log(v) for v in lags],
            is_signal=i in sig_at, direction=sig_at.get(i, ""),
        ))
    return out


# ============================================================================
# 推定（§2 の2段階回帰）
# ============================================================================
@dataclass
class Estimate:
    n_rows: int = 0
    n_signals: int = 0
    n_events: int = 0
    n_days: int = 0
    coef: list[float] = field(default_factory=list)   # [a, b, c]
    ratio: float = 0.0                                # exp(c)
    ci_robust: tuple[float, float] = (0.0, 0.0)
    ci_boot: tuple[float, float] = (0.0, 0.0)
    ci: tuple[float, float] = (0.0, 0.0)              # 広い方
    se_robust: float = 0.0
    har_b: float = 0.0                                # 校正係数（§5 診断1）
    har_b_ci: tuple[float, float] = (0.0, 0.0)
    ci_source: str = ""

    @property
    def verdict(self) -> str:
        return stats_core.effect_verdict(self.ci, self.n_events)

    @property
    def control_verdict(self) -> str:
        return stats_core.control_verdict(self.ci, self.n_events)


def _daily_stats(x: list[list[float]], y: list[float], days: list[int]):
    """日ごとの X'X と X'y を先に畳んでおく（ブートストラップを軽くするため）。"""
    k = len(x[0])
    acc: dict[int, tuple[list[list[float]], list[float]]] = {}
    for row, yi, d in zip(x, y, days):
        xtx, xty = acc.setdefault(
            d, ([[0.0] * k for _ in range(k)], [0.0] * k))
        for i in range(k):
            ri = row[i]
            for j in range(k):
                xtx[i][j] += ri * row[j]
            xty[i] += ri * yi
    return acc


def _cluster_robust_se(x, y, days, beta, k) -> list[float] | None:
    """
    Liang-Zeger の日次クラスタ頑健標準誤差。
    V = (X'X)^-1 (Σ_g X_g'u_g u_g'X_g) (X'X)^-1、小標本補正つき。
    """
    xtx = [[0.0] * k for _ in range(k)]
    for row in x:
        for i in range(k):
            ri = row[i]
            for j in range(k):
                xtx[i][j] += ri * row[j]
    inv = inverse(xtx)
    if inv is None:
        return None

    scores: dict[int, list[float]] = {}
    for row, yi, d in zip(x, y, days):
        u = yi - dot(beta, row)
        sc = scores.setdefault(d, [0.0] * k)
        for i in range(k):
            sc[i] += row[i] * u

    meat = [[0.0] * k for _ in range(k)]
    for sc in scores.values():
        for i in range(k):
            for j in range(k):
                meat[i][j] += sc[i] * sc[j]

    g, n = len(scores), len(y)
    if g < 2 or n <= k:
        return None
    corr = (g / (g - 1)) * ((n - 1) / (n - k))
    v = [[sum(inv[i][a] * meat[a][b] * inv[b][j] for a in range(k)
              for b in range(k)) * corr for j in range(k)] for i in range(k)]
    return [math.sqrt(max(0.0, v[i][i])) for i in range(k)]


def _block_bootstrap_ci(daily, beta_idx: int, rng: random.Random,
                        block: int, iters: int) -> tuple[float, float]:
    """
    日を時系列順に並べ、連続 block 日をひとかたまりとして復元抽出する。

    ★日ごとの十分統計量を先に畳んである理由★
    毎回 X'X を全行から作り直すと 4銘柄×9000本×2000回で現実的な時間に終わらない。
    日単位の X'X と X'y を足し合わせるだけなら、1回あたり数千回の演算で済む。
    """
    keys = sorted(daily)
    nd = len(keys)
    if nd < 10:
        return (0.0, 0.0)
    block = max(1, min(block, nd))
    nb = math.ceil(nd / block)
    k = len(daily[keys[0]][1])

    vals = []
    for _ in range(iters):
        xtx = [[0.0] * k for _ in range(k)]
        xty = [0.0] * k
        for _b in range(nb):
            st = rng.randrange(0, nd - block + 1)
            for d in keys[st:st + block]:
                gx, gy = daily[d]
                for i in range(k):
                    for j in range(k):
                        xtx[i][j] += gx[i][j]
                    xty[i] += gy[i]
        b = solve(xtx, xty)
        if b is not None:
            vals.append(b[beta_idx])
    if len(vals) < 100:
        return (0.0, 0.0)
    vals.sort()
    return (percentile(vals, 0.025), percentile(vals, 0.975))


def estimate(rows: list[Row], rng: random.Random,
             fit_on: str = "control", block_days: int = BLOCK_DAYS,
             iters: int = BOOTSTRAP_ITERS) -> Estimate:
    """§2 の2段階回帰を実行する。"""
    est = Estimate(n_rows=len(rows))
    if len(rows) < 200:
        return est

    # --- 第1段: HAR-RV を対照足のみで学習 ---
    train = [r for r in rows if (not r.is_signal or fit_on == "all")]
    if len(train) < 200:
        return est
    har = ols_fit([r.har_x for r in train], [r.log_fwd for r in train])
    if har is None:
        return est

    # --- 第2段: log RV_fwd = a + b·log(HAR_pred) + c·D ---
    x, y, days = [], [], []
    for r in rows:
        pred = dot(har, r.har_x)          # 既に対数スケール
        x.append([1.0, pred, 1.0 if r.is_signal else 0.0])
        y.append(r.log_fwd)
        days.append(r.day)
    beta = ols_fit(x, y)
    if beta is None:
        return est

    est.coef = beta
    est.har_b = beta[1]
    est.ratio = math.exp(beta[2])
    est.n_signals = sum(1 for r in rows if r.is_signal)
    est.n_days = len(set(days))

    # ★方向を区別せずクラスタ化する（v2.1）★
    # クラスタ化ルールはシグナルではなく【アウトカムの依存構造】に合わせる。
    # 本プロトコルのアウトカムは符号を持たない実現ボラであり、市場全体を動かす
    # ショックが来れば BTC も ETH も同時に高ボラになる。片方が BULLISH、
    # もう片方が BEARISH とタグ付けされていても、その後1時間のボラは同じ
    # ショックに駆動されており独立ではない。方向で分けて2イベントと数えると
    # 標本の独立性を過大評価する。
    # （方向性リターンを見ていた旧プロダクトでは、方向で分けるのが正しかった。
    #   アウトカムが変われば正しいクラスタ化も変わる。）
    sig_rows = [r for r in rows if r.is_signal]
    est.n_events = len(stats_core.build_clusters(
        sorted(r.ts for r in sig_rows), EVENT_WINDOW_BARS))

    se = _cluster_robust_se(x, y, days, beta, 3)
    if se:
        est.se_robust = se[2]
        est.ci_robust = (math.exp(beta[2] - ALPHA_Z * se[2]),
                         math.exp(beta[2] + ALPHA_Z * se[2]))
        est.har_b_ci = (beta[1] - ALPHA_Z * se[1], beta[1] + ALPHA_Z * se[1])

    daily = _daily_stats(x, y, days)
    lo, hi = _block_bootstrap_ci(daily, 2, rng, block_days, iters)
    if hi > lo:
        est.ci_boot = (math.exp(lo), math.exp(hi))

    # --- 広い方を採用（§2）---
    cands = [c for c in (est.ci_robust, est.ci_boot) if c[1] > c[0]]
    if cands:
        widest = max(cands, key=lambda c: c[1] - c[0])
        est.ci = widest
        est.ci_source = ("cluster-robust" if widest == est.ci_robust
                         else "moving-block bootstrap")
    return est


# ============================================================================
# データ組み立て
# ============================================================================
def collect_rows(csv_path: Path, exchange: str | None, timeframe: str,
                 since_ms: int | None = None,
                 placebo_rng: random.Random | None = None,
                 contaminate: bool = False) -> list[Row]:
    """
    placebo_rng を渡すと、検知足を同数のランダムな対照足に置き換える。

    ★contaminate=True は自己検証専用★
    「本物の検知足を対照群に残す」汚染をわざと起こす。
    このとき陰性対照が【必ず失敗する】ことを確認するためのスイッチであり、
    本番経路からは絶対に True にならない。
    """
    signals = load_signals(csv_path)
    if not signals:
        return []
    exchange = exchange or (signals[0].exchange or "bybit")
    out: list[Row] = []
    for symbol in sorted({s.symbol for s in signals}):
        bars = load_bars(exchange, symbol, timeframe)
        if not bars:
            continue
        pos = {int(b[0]): i for i, b in enumerate(bars)}
        real = {pos[s.ts_ms]: s.direction for s in signals
                if s.symbol == symbol and s.ts_ms in pos}
        rows = build_rows(symbol, bars, real)

        if placebo_rng is not None:
            real_idx = {r.ts for r in rows if r.is_signal}
            pool = [r for r in rows if not r.is_signal]
            k = min(len(real_idx), len(pool))
            picked = set(id(r) for r in placebo_rng.sample(pool, k)) if k else set()
            dirs = [r.direction for r in rows if r.is_signal]
            di = 0
            for r in rows:
                was_real = r.is_signal
                r.is_signal = id(r) in picked
                if r.is_signal:
                    r.direction = dirs[di % len(dirs)] if dirs else "BULLISH"
                    di += 1
                elif was_real:
                    r.direction = ""
            if not contaminate:
                # 本物の検知足は対照群からも外す（残すと基準が押し上がる）
                rows = [r for r in rows
                        if r.is_signal or r.ts not in real_idx]
        if since_ms is not None:
            rows = [r for r in rows if r.ts >= since_ms]
        out.extend(rows)
    out.sort(key=lambda r: r.ts)
    return out


# ============================================================================
# 自己検証（§6-2 / §4 の必須要件）
# ============================================================================
def selftest() -> int:
    """
    合成データで、この実装が
      ① 真の効果がゼロなら非有意にする
      ② 本物の効果は検出する
      ③ 汚染された入力では陰性対照が【必ず失敗する】
    ことを確認する。③が無い状態で確認的判定を実行してはいけない。
    """
    rng = random.Random(12345)
    ok = True

    def synth(n_days: int, effect: float, cluster: bool) -> list[Row]:
        rows: list[Row] = []
        level = 0.0
        for d in range(n_days):
            level = 0.9 * level + rng.gauss(0, 0.3) if cluster else rng.gauss(0, 0.3)
            for h in range(24):
                ts = (d * 24 + h) * 3_600_000
                base = level + rng.gauss(0, 0.4)
                is_sig = rng.random() < 0.05
                lag = level + rng.gauss(0, 0.2)
                rows.append(Row(
                    ts=ts, symbol="X", day=ts // MS_DAY,
                    log_fwd=base + (math.log(effect) if is_sig else 0.0),
                    har_x=[1.0, lag, lag + rng.gauss(0, 0.05),
                           lag + rng.gauss(0, 0.05)],
                    is_signal=is_sig,
                    direction="BULLISH" if rng.random() < 0.5 else "BEARISH"))
        return rows

    print("── ① 効果ゼロのデータで非有意になるか ──────────────")
    e = estimate(synth(400, 1.0, True), rng, iters=400)
    v = stats_core.control_verdict(e.ci, e.n_events)
    print(f"   推定 {e.ratio:.3f}倍 CI {e.ci[0]:.3f}-{e.ci[1]:.3f} "
          f"({e.ci_source}) → {v}")
    ok &= (v == "clean")
    print(f"   {'✅ 合格' if v == 'clean' else '❌ 不合格：偽陽性を出している'}")

    print("── ② 本物の効果（1.5倍）を検出できるか ─────────────")
    e2 = estimate(synth(400, 1.5, True), rng, iters=400)
    v2 = stats_core.effect_verdict(e2.ci, e2.n_events)
    print(f"   推定 {e2.ratio:.3f}倍 CI {e2.ci[0]:.3f}-{e2.ci[1]:.3f} → {v2}")
    ok &= (v2 == "effect")
    print(f"   {'✅ 合格' if v2 == 'effect' else '❌ 不合格：本物を見逃している'}")

    print("── ③ 汚染された入力で陰性対照が失敗するか ──────────")
    print("     （対照群に本物の高ボラ足を残すと、偽シグナルが低く出るはず）")
    rows = synth(400, 2.0, True)
    real = [r for r in rows if r.is_signal]
    pool = [r for r in rows if not r.is_signal]
    for r in real:
        r.is_signal = False                      # 本物を対照群に残す＝汚染
    for r in rng.sample(pool, len(real)):
        r.is_signal = True                       # 偽シグナル
    e3 = estimate(rows, rng, iters=400)
    v3 = stats_core.control_verdict(e3.ci, e3.n_events)
    print(f"   推定 {e3.ratio:.3f}倍 CI {e3.ci[0]:.3f}-{e3.ci[1]:.3f} → {v3}")
    ok &= (v3 == "biased")
    print(f"   {'✅ 合格：汚染を検出' if v3 == 'biased' else '❌ 不合格：汚染を見逃す。この状態で判定してはいけない'}")

    print("── ④ 判定ヘルパーが両側で動くか ────────────────────")
    cases = [((1.10, 1.30), "biased", "上振れ"), ((0.80, 0.95), "biased", "下振れ"),
             ((0.95, 1.05), "clean", "正常")]
    for ci, want, name in cases:
        got = stats_core.control_verdict(ci, 100)
        print(f"   {name:<6} CI{ci} → {got} "
              f"{'✅' if got == want else '❌'}")
        ok &= (got == want)

    print()
    print("=" * 60)
    print("  自己検証: " + ("✅ 全項目合格" if ok else "❌ 不合格。判定を実行しないこと"))
    print("=" * 60)
    return 0 if ok else 1


# ============================================================================
# main
# ============================================================================
def load_freeze() -> dict:
    if not FREEZE_FILE.exists():
        return {}
    try:
        return json.loads(FREEZE_FILE.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def main() -> None:
    p = argparse.ArgumentParser(description="事前登録 v2.0 の確認的再現")
    p.add_argument("--csv", type=Path, default=Path("data/spike_log.csv"))
    p.add_argument("--exchange", default=None)
    p.add_argument("--timeframe", default="1h")
    p.add_argument("--exploratory", action="store_true",
                   help="過去データ全体で実行する（探索的。証拠にはならない）")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--out", type=Path, default=Path("replication_status.json"))
    p.add_argument("--seed", type=int, default=20260811)
    args = p.parse_args()

    if args.selftest:
        sys.exit(selftest())

    rng = random.Random(args.seed)
    freeze = load_freeze()
    since = None
    mode = "exploratory"

    if not args.exploratory:
        if not freeze.get("frozen_at"):
            sys.exit("[ERROR] FREEZE.json がありません。\n"
                     "        確認的判定は凍結後にしか実行できません。\n"
                     "        参考値が見たい場合は --exploratory を付けてください。")
        since = int(datetime.strptime(freeze["frozen_at"], "%Y-%m-%dT%H:%M:%SZ")
                    .replace(tzinfo=timezone.utc).timestamp() * 1000)
        mode = "confirmatory"

    rows = collect_rows(args.csv, args.exchange, args.timeframe, since)
    if not rows:
        sys.exit("[ERROR] 対象データがありません。")

    est = estimate(rows, rng)
    prng = random.Random(args.seed + 1)
    pl = estimate(collect_rows(args.csv, args.exchange, args.timeframe, since,
                               placebo_rng=prng), rng)

    reached = est.n_events >= TARGET_EVENTS
    pl_verdict = pl.control_verdict

    print()
    print("=" * 72)
    print(f"  事前登録 v2.0 — {'確認的再現' if mode == 'confirmatory' else '探索的（証拠ではない）'}")
    print("=" * 72)
    if mode == "confirmatory":
        print(f"  凍結時刻: {freeze.get('frozen_at')}")
        print(f"  params_hash: {freeze.get('params_hash', '—')}")
    print(f"  独立イベント {est.n_events} / {TARGET_EVENTS} 件"
          f"（検知 {est.n_signals} 件 / 全足 {est.n_rows} 行 / {est.n_days} 日）")
    print(f"  exp(c) = {est.ratio:.3f}倍  CI {est.ci[0]:.3f}–{est.ci[1]:.3f}"
          f"（{est.ci_source} / 広い方を採用）")
    print(f"    ├ クラスタ頑健  {est.ci_robust[0]:.3f}–{est.ci_robust[1]:.3f}")
    print(f"    └ ブロックBS    {est.ci_boot[0]:.3f}–{est.ci_boot[1]:.3f}")
    print(f"  HAR校正係数 b = {est.har_b:.3f} "
          f"[{est.har_b_ci[0]:.3f}, {est.har_b_ci[1]:.3f}]"
          "  ← 1.0から大きく外れるなら §5 診断を参照")
    print(f"  陰性対照: {pl.ratio:.3f}倍 "
          f"[{pl.ci[0]:.3f}–{pl.ci[1]:.3f}] → {pl_verdict}")
    print()

    if mode == "exploratory":
        print("  ⚠️  これは探索的な参考値です。確認的な証拠ではありません。")
        print("     LP・マーケ文面では必ず Exploratory と明記してください。")
    elif not reached:
        print(f"  ⏳ 中間値です。判定は {TARGET_EVENTS} 件到達時に1回だけ行います。")
        print("     この数字を見て止める・続けるの判断をしてはいけません。")
    elif pl_verdict != "clean":
        print("  🛑 陰性対照が汚染しています。判定を保留し、原因を特定してください。")
    else:
        v = est.verdict
        if v == "effect":
            print(f"  ✅ 再現成功。exp(c) の95%CI下限 {est.ci[0]:.3f} > 1.0")
        else:
            print(f"  ❌ 再現失敗。CI {est.ci[0]:.3f}–{est.ci[1]:.3f} は 1.0 を含みます。")
            print("     §7 の報告義務を実行してください（LPトップに掲載・")
            print("     ゲートを replication_failed に・課金者へ通知・返金提示）。")
    print("=" * 72)
    print()

    payload = {
        "mode": mode,
        "generated": f"{datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}",
        "frozen_at": freeze.get("frozen_at", ""),
        "params_hash": freeze.get("params_hash", ""),
        "target_events": TARGET_EVENTS,
        "n_events": est.n_events,
        "n_signals": est.n_signals,
        "reached": reached,
        "interim": (mode == "confirmatory" and not reached),
        "ratio": round(est.ratio, 4),
        "ci": [round(est.ci[0], 4), round(est.ci[1], 4)],
        "ci_source": est.ci_source,
        "ci_robust": [round(est.ci_robust[0], 4), round(est.ci_robust[1], 4)],
        "ci_bootstrap": [round(est.ci_boot[0], 4), round(est.ci_boot[1], 4)],
        "har_calibration_b": round(est.har_b, 4),
        "har_calibration_ci": [round(est.har_b_ci[0], 4),
                               round(est.har_b_ci[1], 4)],
        "placebo": {"ratio": round(pl.ratio, 4),
                    "ci": [round(pl.ci[0], 4), round(pl.ci[1], 4)],
                    "verdict": pl_verdict},
        "verdict": (est.verdict if (mode == "confirmatory" and reached
                                    and pl_verdict == "clean") else "pending"),
    }
    args.out.write_text(json.dumps(payload, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(f"[DONE] {args.out}")


if __name__ == "__main__":
    main()
