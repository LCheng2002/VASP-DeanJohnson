# Mg2NiHx / LaH3 extension batch — actual Slurm submission

**Status:** submitted successfully to Siyuan-1; **not yet finished**. The initial queue verification was performed on **2026-10-12 14:29:27 CST (China Standard Time)**.

| Slurm Job ID | Group | Scope | Verified at submission |
|---|---|---|---|
| 63585911 | parents | Array 0–4, maximum 2 concurrently: L_parent, L_MgH2_reference, Mg_parent, X025_A, X025_B | PENDING |
| 63585912 | H2 | Matching H2 molecule reference | PENDING |
| 63585913 | La | Five MgH2/LaH3 defect and transfer endpoints; starts only after parents succeed | PENDING (dependency) |
| 63585914 | Mg | Two pure-Mg H2/2H endpoints; starts after parents and H2 succeed | PENDING (dependency) |
| 63585915 | Hx | Two sequential Mg2NiHx hydrogen-transfer endpoints; starts only after parents succeed | PENDING (dependency) |

## Follow-up queue check — 2026-10-12 14:32:34 CST

Verified through the **read-only** [status workflow run 38213015746](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38213015746):

- **63585912 (H2): RUNNING** on `node014`; elapsed 00:02:42 at check.
- **63585911_0–_4 (all five parents): PENDING (Priority)**, not yet started.
- **63585913, 63585914, 63585915: PENDING (Dependency)**, expected until the five parent tasks (and H2 for Mg) satisfy their dependencies.
- All 15 model entries still reported `NO_FINAL_VALIDATION` at this snapshot. **No VASP final result is claimed.**

The original table above records the immediately-after-submission state, not this later queue state.

## Parent array concurrency cap removed — 2026-10-12 19:01:37 CST

The existing job array **63585911** was updated in place; its index roster remains **0–4** with five parent models.

- Before: `ArrayTaskId=0-4%2` (maximum two concurrent).
- Applied: `scontrol update JobId=63585911 ArrayTaskThrottle=0`.
- After: `ArrayTaskId=0-4` (no job-array-level concurrency throttle).
- Immediately afterward, tasks `63585911_0` through `63585911_4` were still **PENDING (None)**; available scheduler resources and system policy still determine start times.
- **No cancellation, resubmission, model-count change, dependency rewrite or relaxation of numerical convergence criteria.** The repository future-submission template also now specifies `0-4`.

Evidence: [successful GitHub Actions live Slurm update 38239984463](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38239984463).

The calculation directory is:
`~/MgNi_LaH3_extensions_20261012/`

The staged GitHub source/input directory is:
`~/VASP/github_jobs/MgNi_LaH3_extensions_repo_20261012/`

The parent array plus H2 totals six initial models. The nine dependent child structures **have not been calculated**; they will be constructed only after the applicable parent models pass their Relax/Static validation gates. Total authorized scope: 15 models; **no NEB**, DOS/Bader sweeps, AIMD, or extension to further Mg17La2 coverage.

**Submission evidence:** [GitHub Actions workflow run 38212693610](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38212693610). GitHub runner exited successfully, remote `sbatch --parsable` printed the actual five IDs, and read-only `sacct` confirmed the five records in PENDING state. This is **not** evidence that VASP has converged or that final output exists yet.

Only the approved runner, plan, and three explicitly pinned structure files were transmitted via SSH. POTCAR was not uploaded to GitHub, and the licensed potentials are assembled inside the HPC account. The driver will not overwrite an existing submitted batch if its submission record is present.

To check inside Siyuan, without submitting again:

```bash
python3 "$HOME/MgNi_LaH3_extensions_20261012/runner.py" --status
```

Notes: low-H `Mg2NiHx` seed is an `x=0.25` finite-cell proxy for the experimentally discussed `Mg2NiH0.3`, not a crystallographically refined `Mg2NiH0.3` structure. All structures await optimization. Under the physical settings of this batch (spin-polarized PBE without U), compare energies only within consistently parameterized groups. Refer to [README](README.md) for model design, limitations and references.
