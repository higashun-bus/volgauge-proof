#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
power_calc.py — 確認的再現に必要な独立イベント数を決める
=========================================================

★点推定を真値と仮定して検出力を設計してはいけない（勝者の呪い）★

探索段階の推定値は、3つの horizon × 2つのベースライン × 複数の仕様を見た後に
選ばれたものです。選択された推定値は系統的に上振れします
（winner's curse / regression to the mean）。
点推定で設計すると、再現テストは構造的に under-power になります。

そして under-power したテストが失敗したとき、
「効果が無かった」のか「検出力が足りなかった」のか区別がつきません。
この曖昧さが、このプロジェクトが積み上げた信頼を最も損ないます。

さらに事前登録 §7 により、再現失敗は
【LPトップへの掲載・課金者への通知・返金提示】という不可逆な対応を自動的に発動します。
検出力の不足は統計上の不備ではなく、**自分のプロダクトを誤って停止させる事業リスク**です。

したがって本ツールは【CI下限を真値と仮定したときに検出力80%】となる件数を返します。

使い方:
    python power_calc.py --mult 1.312 --lo 1.215 --hi 1.416 --n 872 --days 385
"""

from __future__ import annotations

import argparse
import math

Z_ALPHA = 1.959963984540054      # 両側 95%


def norm_cdf(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def power_at(n: int, c_true: float, se_ref: float, n_ref: int) -> float:
    """
    件数 n のときの検出力（両側 α=0.05）。

    SE は 1/sqrt(n) で縮むと仮定する。探索段階の (n_ref, se_ref) を基準に
    スケールさせる。両側なので反対側の棄却域もわずかに足す。
    """
    if n <= 0:
        return 0.0
    se = se_ref * math.sqrt(n_ref / n)
    t = c_true / se
    return norm_cdf(t - Z_ALPHA) + norm_cdf(-t - Z_ALPHA)


def main() -> None:
    p = argparse.ArgumentParser(
        description="CI下限を真値と仮定して必要な独立イベント数を決める")
    p.add_argument("--mult", type=float, required=True, help="探索段階の点推定（倍率）")
    p.add_argument("--lo", type=float, required=True, help="95%CI下限")
    p.add_argument("--hi", type=float, required=True, help="95%CI上限")
    p.add_argument("--n", type=int, required=True, help="探索段階の独立イベント数")
    p.add_argument("--days", type=int, required=True, help="探索段階の観測日数")
    p.add_argument("--power", type=float, default=0.80, help="目標検出力")
    args = p.parse_args()

    c_hat = math.log(args.mult)
    se = (math.log(args.hi) - math.log(args.lo)) / 2 / Z_ALPHA
    t = c_hat / se
    c_lo = math.log(args.lo)
    rate = args.n / args.days * 30.44        # 件/月

    print()
    print("=" * 70)
    print("  確認的再現に必要な独立イベント数")
    print("=" * 70)
    print(f"  探索段階: {args.mult:.4f}倍  CI [{args.lo:.4f}, {args.hi:.4f}]")
    print(f"            独立イベント {args.n} 件 / {args.days} 日")
    print(f"  c = ln(倍率) = {c_hat:.5f}   SE = {se:.5f}   t = {t:.3f}")
    print(f"  発生率 = {rate:.1f} 件/月")
    print()

    # 目標: CI下限が真値でも検出力 args.power
    need = 0
    for n in range(30, 5001):
        if power_at(n, c_lo, se, args.n) >= args.power:
            need = n
            break

    print("  ── 検出力の一覧 " + "─" * 48)
    print(f"  {'n':>6}{'点推定が真値':>16}{'CI下限が真値':>16}"
          f"{'偽陰性率':>12}{'到達':>10}")
    marks = sorted({150, 170, 200, 250, 300, 350, 400, need})
    for n in marks:
        p1 = power_at(n, c_hat, se, args.n)
        p2 = power_at(n, c_lo, se, args.n)
        mo = n / rate if rate > 0 else 0
        flag = "  ← 確定" if n == need else ""
        print(f"  {n:>6}{p1*100:>15.1f}%{p2*100:>15.1f}%"
              f"{(1-p2)*100:>11.1f}%{mo:>9.1f}ヶ月{flag}")
    print()
    print(f"  ✅ 判定件数: 独立イベント {need} 件")
    print(f"     （CI下限 {args.lo:.3f}倍 が真値でも検出力 "
          f"{power_at(need, c_lo, se, args.n)*100:.1f}%）")
    print(f"     到達予測: 約 {need/rate:.1f} ヶ月")
    print()
    print("  ⚠️  この件数は凍結対象です。到達前に増やすことも減らすことも")
    print("      事前登録違反です。判定を先送りする理由を後から作らないこと。")
    print("=" * 70)
    print()


if __name__ == "__main__":
    main()
