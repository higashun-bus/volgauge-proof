# Errata — pre_registration_v2.md

This file records errors found in the frozen pre-registration **after** it was frozen and
timestamped. It is a separate file by design. See "Why the frozen file is not edited" below.

- Pre-registration frozen at: `2026-08-11T06:43:19Z`
- Timestamped in Bitcoin block: `961982`
- Verdict at: `267` independent events
- This errata first published: `2026-08-17`

Nothing in this file changes the decision rule. Every entry states explicitly what it does
**not** affect, and that scope statement is the substance of the entry.

---

## E1 — §8 pairs the 1.33× figure with the wrong event count

**Location.** `pre_registration_v2.md`, §8 *LP 表示仕様* → *ヒーロー数値* (line 288).

**What it says:**

> **1.33×**（補正後の素の倍率、95% CI 1.22–1.45、**807 independent events**）

**What is correct.** The 1.33× multiple with 95% CI 1.22–1.45 comes from **853** independent
events, not 807. Source: `baseline_report.md`, §1, the 1h row —
`1h | 1513 raw | 853 independent | 1.38× uncorrected [1.30–1.47] | 1.33× corrected [1.22–1.45] | p<0.001`.

**Where 807 actually comes from.** It is the independent-event count of a **different row** in a
**different table** — the "no baseline (comparison against the per-asset median only)" model,
which yields **1.44×** over 807 events (`baseline_report.md`, §2). The two numbers were
transposed when §8 was drafted.

**Additional note against ourselves.** `807` appears on this project's own list of numbers that
must not be used (`PRODUCT_CLAIM_EN.md` §9, "forbidden numbers"), because it belongs to a
superseded draft. The frozen pre-registration therefore contains a figure that the project's own
publication rules prohibit. We are recording that rather than explaining it away.

**Direction of the error.** Against us. 807 is **lower** than the correct 853, so the frozen
document understates the sample size behind its own headline figure. The error makes the
evidence look weaker, not stronger. We note this because the direction is the first thing a
reader should be able to check, not because it makes the error acceptable.

**Discovered.** 2026-08-17, during a full machine re-check of every hash, block number and
event count appearing in any published artifact. Found by us, not reported to us.

### What E1 does not affect

This is the part that matters. E1 sits in **§8, which is display specification for the landing
page.** It does not touch any part of the protocol that decides the outcome.

| Section | Content | Affected by E1 |
|---|---|---|
| **§1** | The single hypothesis under confirmation | **No** |
| **§2** | Complete definition of the primary endpoint, including the independent-event clustering rule | **No** |
| **§3** | Sample size and stopping rule — the verdict at **267** independent events | **No** |
| **§4** | Negative controls (placebo) | **No** |
| **§5** | Pre-specified sensitivity analyses | **No** |
| **§7** | Reporting obligations, including what happens if the result fails | **No** |
| §8 | Landing-page display specification | **Yes — this entry** |

**The decision rule is untouched.** 267 is set in §3 and is computed by `power_calc.py` from the
CI lower bound 1.215 over 843 events — a different quantity from the 1.33× row and unrelated to
the 807/853 mix-up. Neither 807 nor 853 is an input to the stopping rule.

**No published claim used 807 with 1.33×.** Every downstream document pairs 1.33× with 853
(`PRODUCT_CLAIM_EN.md`, `MARKETING_LAUNCH_EN.md`). The error is confined to the frozen file.

### Why the frozen file is not edited

`pre_registration_v2.md` is one of four files whose SHA-256 is recorded in `FREEZE.json`, and
`FREEZE.json` itself is timestamped into Bitcoin block 961982. **Editing the file would change
its hash, break the chain, and destroy the proof that the criteria predate the data.** The
proof is worth more than the typo.

That is the whole point of freezing: the document is out of reach of the person it is about,
including when that person would like to fix something. **An errata file is the correct
instrument. Editing is not available and should not be.**

Verify for yourself that the frozen file still matches what was timestamped:

```
sha256  pre_registration_v2.md
  = 6a0fab6694562049004e38326f518a97c9bdaa1d1ca4be5e15587fa8c0f795b4
```

That value is recorded in `FREEZE.json`, whose own hash is
`84b754c27b655fff3917a0dc572214083070d829068bba251be9b3cc720229cc`, and whose existence before
block 961982 is provable from `FREEZE.json.ots` without trusting us. See the repository README
for three ways to check, one of which requires no software.

---

## Scope of this file

Entries are added when an error is found in the frozen pre-registration. Entries are never
removed or edited once published; corrections to an entry are appended as a new entry that
supersedes it, and the superseded entry stays.

**This file is not itself frozen.** It is intended to grow. Its hash is published on each
update so that any version of it can be pinned, but a reader who wants a guarantee about the
*criteria* should verify `pre_registration_v2.md` and `FREEZE.json`, not this file.

## Why an errata file exists at all

A pre-registration with no errata has not necessarily been checked. It may simply mean nobody
read it again. **This one was read again, four months before the verdict, and the error found
was in the direction that weakens our own claim.** Publishing it now rather than in December is
the difference between an audit and a defence.
