# Do not resubmit: duplicate Ca OH P0 workflow (2026-10-15)

This folder was prepared as a backup way to submit the same four archived Ca–OH **fixed-geometry Static dipole calibration** models. It created Slurm array **63693436**, but at that point the first submission **63686231** had already fully completed.

Reconciliation [Actions run 38669235325](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38669235325) verified all four original models with `CALIBRATION_COMPLETE` and `COMPLETED (0:0)`, and immediately cancelled all four duplicate tasks `63693436_0–3`. The original files were preserved.

**Canonical complete work is here:** [`Ca_OH_dipole_calibration_20261015/STATUS.md`](../Ca_OH_dipole_calibration_20261015/STATUS.md). Do not restart/reuse this duplicate directory or manually dispatch [`submit-ca-oh-p0-calibration-20261015.yml`](../../.github/workflows/submit-ca-oh-p0-calibration-20261015.yml). No new P1/P2 or previously cancelled CaFresh jobs were submitted by this reconciliation.

The duplicate P0 setup uses a slight variation in the in-plane DIPOL convention; because `IDIPOL=3`, the intention was to calibrate the slab-normal component. The original completed batch includes explicit NBANDS matching and should be retained as the source of truth.

**Next:** read original `~/Ca_OH_dipole_calibration_20261015/*/validation.json` and its P0 comparison before planning short relaxation or further corrosion chemistry.
