# Mg2NiHx / LaH3 extension batch — actual Slurm submission

**Status:** submitted successfully to Siyuan-1; **not yet finished**. The initial queue verification was performed on **2026-10-12 14:29:27 CST (China Standard Time)**.

| Slurm Job ID | Group | Scope | Verified at submission |
|---|---|---|---|
| 63585911 | parents | Array 0–4, maximum 2 concurrently: L_parent, L_MgH2_reference, Mg_parent, X025_A, X025_B | PENDING |
| 63585912 | H2 | Matching H2 molecule reference | PENDING |
| 63585913 | La | Five MgH2/LaH3 defect and transfer endpoints; starts only after parents succeed | PENDING (dependency) |
| 63585914 | Mg | Two pure-Mg H2/2H endpoints; starts after parents and H2 succeed | PENDING (dependency) |
| 63585915 | Hx | Two sequential Mg2NiHx hydrogen-transfer endpoints; starts only after parents succeed | PENDING (dependency) |

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
