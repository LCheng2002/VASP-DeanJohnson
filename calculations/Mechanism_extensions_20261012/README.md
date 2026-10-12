# Mg2NiHx / LaH3 / pure-Mg mechanism extensions — submitted, results pending

## Actual status

**Submitted to Siyuan-1 on 2026-10-12 at 14:29:27 CST (Beijing time).** GitHub Actions run [38212693610](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38212693610) returned five genuine Slurm job IDs: **63585911** (5-parent array), **63585912** (H2), **63585913** (La branch), **63585914** (pure-Mg branch), and **63585915** (low-H Mg2NiHx branch). The Slurm records were **PENDING** immediately after submission; this does not imply that VASP has finished.

See [STATUS.md](STATUS.md) for the verified timestamped submission record, HPC paths and branch dependencies. No pseudopotentials, keys or licensed data were transferred to GitHub.

The previous GitHub-only input validation run **38205762176** passed and produced the six preview seeds; that preliminary test was separate from the later **actual Slurm submission**. All new VASP Relax and Static outputs remain subject to their convergence checks.

## Literature and the precise intermediate being tested

Čermák and Král, *Ageing of Mg–Ni–H hydrogen storage alloys*, International Journal of Hydrogen Energy 37 (2012), 14257–14264, DOI **10.1016/j.ijhydene.2012.07.049**, discuss residual low-H **Mg2NiH0.3** associated with MgH2 and its hydrogenation using hydrogen from adjoining MgH2. This is experimental phase-evolution evidence and a proposed interpretation, not an atom-resolved proof of a catalytic cycle.

Buzanov et al., *Preparation and Characterization of MgH2 Mechanocomposites with Mg2NiH0.3 + Mg2NiH4−δ Two-Phase Mixture*, Russian Journal of Inorganic Chemistry 63 (2018), 1529–1533, DOI **10.1134/S0036023618120057**, explicitly distinguish **Mg2NiHx (0<x<0.3)** from **Mg2NiH4−δ**.

Primary sources:
- https://www.sciencedirect.com/science/article/abs/pii/S0360319912016266
- https://www.ceitec.eu/ageing-of-mg-ni-h-hydrogen-storage-alloys/p110597
- https://link.springer.com/article/10.1134/S0036023618120057

Consequently, the previous rich-hydride two-vacancy model is NOT relabeled Mg2NiH0.3. This batch adds dilute H to the **Mg2Ni metal framework** instead of deleting almost all H from an Mg2NiH4 framework.

## Model roster: six actual seed inputs; nine dependent models

| Group | Model count including its parents | Scientific question |
|---|---:|---|
| MgH2/LaH3 and matched finite MgH2 cleavage control | 7 | How does the interface change selected H-vacancy endpoints, and can an existing MgH2-side H fill a La-side vacancy? |
| MgH2/Mg2NiHx low-H proxy | 4 | Can a low-H Mg2Ni domain accept one and then a second H from adjacent MgH2? |
| Pure Mg cleavage control and H2/2H endpoints | 3 | Prepare a genuine no-Ni finite-surface comparison for molecular/dissociated H endpoints |
| H2 reference | 1 | Supply a matching H chemical-potential reference |
| **Maximum total** | **15** | Relax→Static only; no NEB |

### Six generated seeds

| Name | Composition | Atoms | Source/meaning |
|---|---|---:|---|
| `L_parent` | Mg12La16H72 | 100 | Supplied MgH2(100)/LaH3(110) CONTCAR, to be reoptimized under the new settings |
| `L_MgH2_reference` | Mg12H24 | 36 | The MgH2 component cut from the same interface cell; finite, strained cleavage control |
| `Mg_parent` | Mg24 | 24 | Mg(10-10) cleavage component of the existing A interface; not an Ni-to-Mg substitution model |
| `X025_A` | Mg60Ni12H75 | 147 | 1×2 C-interface supercell with three dispersed background H candidates in the Mg2Ni domain |
| `X025_B` | Mg60Ni12H75 | 147 | Same cell/composition with a different background-H arrangement |
| `H2_reference` | H2 | 2 | 20 Å cubic box, initial H-H=0.75 Å |

The six seed POSCARs are **unrelaxed initial structures**, not results. The previous archived CONTCAR is a geometry source, not proof that this new spin-polarized batch has converged. Exact source paths and blob hashes are in `plan.json`.

### Why x=0.25 rather than exactly 0.3

The low-H interface contains 36 MgH2 formula units and 12 Mg2Ni formula units. Three background H atoms in the latter give **x=3/12=0.25**. The total composition is Mg60Ni12H75: 72 H initially belong to MgH2, and three to the low-H receiving domain.

This is a finite-cell low-H **proxy within the reported dilute-composition range**, not an exact Mg2NiH0.3 crystal structure, nor a refinement of experimentally measured H positions. The two arrangements sample distinct geometrical candidates rather than assuming a known ordered H sublattice. Agreement of two configurations does not establish a global minimum or configurational convergence.

`X_transfer_one` and `X_transfer_two` move existing MgH2-side H atoms into safe receiving candidates; they never add external H. Total Mg60Ni12H75 remains unchanged. Nominal receiving-domain H counts would become 4/12 and 5/12 only if the optimized structures confirm the transfer. These are uptake endpoints, **not claims of stable homogeneous alpha phases at x>0.3**.

### LaH3-dependent endpoints

After `L_parent` and the finite MgH2 reference pass validation:

- `L_vMg_near`: remove original H39 in the MgH2-side near-interface region.
- `L_vMg_far`: remove original H31 at the finite slab's other-side region.
- `L_vLa`: remove original La-derived H55.
- `L_MgH2_vMg`: remove the corresponding H39 from the matched MgH2 cleavage control.
- `L_transfer`: from the actual converged `L_vLa`, move an eligible original MgH2-side H into the translated vacancy candidate; conserve Mg12La16H71.

The fixed labels identify the original sites and must be checked against final coordination. Only one La-derived vacancy is sampled. The finite far-side position is not a bulk limit. Interface loss, excessive reconstruction, or failure to retain the intended transfer endpoint requires interpretation or redesign, not automatic mechanistic claims.

### Pure-Mg and low-H dependencies

After `Mg_parent` converges, `Mg_H2` and `Mg_2H` place a molecular or pre-separated pair above its exposed Mg region and allow relaxation. These are potential path endpoints, not dissociation barriers.

After both `X025_A/B` converge, the lower-energy candidate is selected only when TOTEN and sigma→0 rankings agree and the background H remains in the intended Ni-neighbor environment. `X_transfer_one` is then built and optimized; only after it passes can `X_transfer_two` be constructed. All H identities and actual final metal coordination must be inspected before interpreting a transfer as successful.

## Numerical settings and limitations

All NEW references, parents and descendants use spin-polarized PBE (**ISPIN=2**), no U, ENCUT=520 eV, ISMEAR=0, SIGMA=0.05, LREAL=.FALSE., EDIFF=1e-7, fixed cell and inherited atomic constraints. Starting magnetic moments are Ni=1.0 μB, La=0.1 μB and zero elsewhere; this is not an exhaustive magnetic-order search. H2 starts with zero moments.

Relax: IBRION=2, POTIM=0.10, NSW=400, EDIFFG=-0.015. Static: NSW=0 and the same force-evaluation settings. Acceptance requires final electronic convergence, a normal-end marker, conserved composition/cell/constraints and maximum free-atom force below 0.015 eV/Å for Relax and 0.02 eV/Å for Static. H2 uses EDIFF=1e-8 and force threshold 0.005 eV/Å.

K meshes: La group 3×5×1, Mg group 2×5×1, expanded low-H group 3×3×1, H2 Gamma only. These are initial group-consistent settings, **not claims of k-point, vacuum or slab-thickness convergence**. The thin cleavage controls and lack of an explicit dipole-convergence study limit quantitative interpretation, especially for small energy differences.

The prior Mg–Ni screening used non-spin-polarized settings. Do not mix old and new bare energies when claiming an effect of Ni, LaH3 or hydride state. This batch adds missing parent/defect controls; old Ni endpoint energies need a matching numerical comparison before a quantitative pure-Mg versus Ni effect is asserted.

The intended output is endpoint thermodynamics. No H2 dissociation or diffusion barrier, explicit gas-exit trajectory, high-temperature phase transformation, complete low-H→Mg2NiH4 conversion, or repeated catalytic cycle is calculated here. No new Mg17La2 fifth-H calculation, large three-phase junction or DOS/Bader sweep is included.

## Cluster-side execution and tracking

The previous manual-only instructions have been superseded by **successful GitHub Actions submission**. The input stage was deployed to `~/VASP/github_jobs/MgNi_LaH3_extensions_repo_20261012/` and VASP production outputs will be written within `~/MgNi_LaH3_extensions_20261012/`.

Job IDs and dependent steps are recorded in [STATUS.md](STATUS.md), along with the [successful GitHub Actions workflow run](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38212693610).

The following on-cluster command is read-only and does **not** submit again:

```bash
python3 "$HOME/MgNi_LaH3_extensions_20261012/runner.py" --status
```

All five Slurm submission calls were acknowledged and recorded. The parent array was **originally submitted** with a maximum of two concurrent members, but that restriction was **removed live on 2026-10-12 at 19:01 CST**, using `scontrol update JobId=63585911 ArrayTaskThrottle=0` ([verified Actions run 38239984463](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38239984463)). The live `ArrayTaskId` changed from `0-4%2` to `0-4`: all five parent tasks are eligible to start when resources permit. The future source template was also updated from `0-4%2` to `0-4`; **no Slurm job was cancelled or resubmitted**. The dependent La, Mg and low-H branches do not run until their upstream models have met the Relax/Static validation criteria. An upstream failure may cancel a dependent branch; this avoids producing unreliable continuation calculations. There is no automatic retry of failed physics or numerical checks.
