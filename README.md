# volgauge-proof

Pre-registration and analysis code for VolGauge, frozen before any of the data
that will decide the result was collected.

- **Frozen:** 2026-08-11T06:43:19Z
- **Timestamped in Bitcoin block:** 961982
- **Verdict at:** 267 independent events (~Dec 2026)
- **params_hash:** `dbf36471635b`

The same hash was posted at the time of freezing:
https://x.com/VolGauge/status/2087396150265749596

## What is here

| File | What it is |
|---|---|
| `pre_registration_v2.md` | The pre-registration itself |
| `FREEZE.json` | SHA-256 of every frozen file, the detection parameters, and their hash |
| `FREEZE.json.ots`, `pre_registration_v2.md.ots` | OpenTimestamps proofs, upgraded to Bitcoin |
| `stats_core.py`, `replication.py`, `power_calc.py` | The analysis code that will decide the result |

Published here: the analysis code that will decide the December result.
Not published: whale_alert.py (detection). Its parameters are fixed in
the pre-registration and locked by params_hash dbf36471635b. Changing
them resets the event count.

## Verify this yourself

**A. Easiest — no software required**

Download `FREEZE.json` and `FREEZE.json.ots`, then drop both into
<https://opentimestamps.org>. It verifies server-side.

**B. Without a Bitcoin node**

```
pip install opentimestamps-client
ots info FREEZE.json.ots
```

This shows `BitcoinBlockHeaderAttestation(961982)`. Look up block 961982 on any
block explorer and read its timestamp.

**C. Fully trustless — requires Bitcoin Core**

```
ots verify FREEZE.json.ots
```

This confirms `FREEZE.json` existed before Bitcoin block 961982
(2026-08-11T08:17:15Z). The hashes inside it identify the exact analysis code
and pre-registration that will decide the December result.

**Check the code against the pre-registration**

The `files` section of `FREEZE.json` records the SHA-256 of each analysis file.
Hash the `.py` files in this repository and compare. They must match, or the
pre-registration does not describe what is published here.
