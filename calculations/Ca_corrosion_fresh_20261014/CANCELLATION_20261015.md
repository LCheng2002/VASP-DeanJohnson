# CaFresh 23-job Slurm cancellation — 2026-10-15

**Execution:** GitHub Actions [run 38599214913](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38599214913), SUCCESS.

**Cluster time:** 2026-10-15 00:34:55 CST (China Standard Time, UTC+08:00).

The cluster script validated the production `~/Ca_corrosion_fresh_20261014/submission.json` full 23-entry exact job-name/Job-ID mapping, preserved all files and only passed **20 still-active jobs** to `scancel`. All 20 belonged to the same HPC account and had the exact `CaFresh_*` job names in the verified ledger.

**scancel invoked on these 20 jobs:**

```text
63682000 63682002 63682003 63682004 63682005 63682006
63682008 63682009 63682010 63682011 63682012 63682013
63682014 63682015 63682016 63682021 63682022 63682023
63682024 63682026
```

The script returned:
```text
SCANCEL_ACCEPTED|count=20
AFTER_ACTIVE|{}
CANCELLATION_VERIFIED|All_23_batch_ids_absent_from_squeue
OUTPUT_FILES_PRESERVED|No files modified or removed
```

Three reference jobs were already inactive and were **not** passed to `scancel`:
- `63681995 CaFresh_Mg_bulk`: `COMPLETED (0:0)` at 2026-10-14 22:21:44.
- `63681996 CaFresh_H2`: `COMPLETED (0:0)` at 2026-10-14 22:38:55.
- `63681997 CaFresh_H2O`: `FAILED (1:0)` at 2026-10-14 22:21:12.

The immediate sacct read still displayed some newly cancelled pending jobs as `PENDING`, a snapshot not yet updated by accounting. The independently verified **empty squeue for these exact IDs** confirms no matching jobs remained queued or running immediately after `scancel`.

**Scope:** Only `Ca_corrosion_fresh_20261014`, not `Ca_corrosion_20261014/verified_parents_20261014`, LaH3, MgNi, Mg2NiHx, or H2-dissociation NEB. No VASP inputs, output files, or existing completed jobs were altered. Future Ca experiments must be explicitly re-authorized rather than silently resubmitted.
