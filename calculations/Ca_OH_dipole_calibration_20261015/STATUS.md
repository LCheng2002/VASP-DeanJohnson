# Ca OH 四点偶极校准：真实提交记录

**已于 2026-10-15 00:55:14（北京时间，UTC+08:00）通过 GitHub Actions 提交思源一号。**

Slurm 数组 **63686231**，索引 **0–3**，共4个固定几何 Static。

最近本轮 Slurm 核验时间：**2026-10-15 00:55:24**。这是状态快照，不会自行更新。

| Slurm Job ID | 模型 | 组成 | 核查状态 |
|---|---|---|---|
| 63686231_0 | M0_OH | Mg96O1H1 | PENDING (Priority) |
| 63686231_1 | M1_OH | Mg95Ca1O1H1 | PENDING (Priority) |
| 63686231_2 | M0_OH_H2O | Mg96O2H3 | PENDING (Priority) |
| 63686231_3 | M1_OH_H2O | Mg95Ca1O2H3 | PENDING (Priority) |

`squeue` 已核对四个数组成员，`sacct` 返回 `63686231_[0-3]|PENDING`。四项均为 `NOT_YET_CALIBRATED`，**目前没有新的校准能量或受力结果**。不得把GitHub工作流SUCCESS当作VASP已经收敛。

- 成功提交及Slurm核查：[GitHub Actions run 38602080764](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38602080764)。
- 实际生产目录：`~/Ca_OH_dipole_calibration_20261015/`。
- 部署输入来源：`~/VASP/github_jobs/Ca_OH_dipole_calibration_20261015/`。
- 提交回执：生产目录 `submission.json`。
- GitHub artifact：11914880245，`Ca-OH-four-point-calibration-inputs-and-receipt-20261015`。
- 回收artifact ZIP SHA-256：`9365e6054aea723ad5d8e463436fcf2a8989e884ad0e263f8b779e0a96728892`。

真实回执：
```text
SUBMITTED|CaOH_P0|63686231|array=0-3
```

## 本轮只做什么

保留四个原始收敛CONTCAR中的坐标、晶胞、组成及32个固定底层原子，以统一规则的直接坐标DIPOL重新做单点。PBE、ISPIN=1、ISMEAR=1/SIGMA=0.2、ENCUT=520、Gamma 3×3×1、LREAL=F、EDIFF=1e-7及原数据集不改变；NBANDS固定为各旧Static的实际值。

统一 `LDIPOL=T; IDIPOL=3`，不做原子优化。每个作业64 MPI、128GB、4小时上限，不设额外数组并发限制。

校准将分别报告电子收敛、受力与真空电荷/电势检查。若新自由原子力>0.02 eV/Å，标注需要短程再弛豫，不将这一单点诊断误报为已获得新平衡构型；也不把它直接叫作VASP崩溃。没有自动提交弛豫、NEB或更多覆盖度模型。

## 原 O 赝势一致性问题已解决

第一次运行 38601441502 在提交前因原 O 数据集哈希未匹配而主动停止，**没有产生 Slurm Job ID，也没有运行VASP**。缓存O的SHA-256为a9b088…，并非这批旧结果的818f92…。

之后只在旧Ca作业63131486/63306398对应目录内查找，找到完全匹配的原O数据集；Mg/Ca_sv/O/H各块及最终组合POTCAR全部按原归档SHA-256通过检查，随后才提交63686231。没有改变预期赝势哈希来绕过检查，没有将POTCAR上传GitHub，没有覆盖旧赝势或旧结果。

## 与取消操作的边界

**本次没有恢复已取消的CaFresh 20个排队作业，没有重建23个新模型，也没有修改其他MgNi/LaH3/NEB任务。**

旧引用方案中“利用已提交的63682006/63682021和覆盖度任务”已被用户的取消决定取代：这些不再是待运行计算。P0之后，根据校准结果判断需不需要少量再弛豫，再在同一旧结构系列内构建第二分子水解离及OH覆盖候选；它们本轮没有自动排队。

只读查询：
```bash
python3 ~/Ca_OH_dipole_calibration_20261015/calibrate.py --status
```

## Latest verification — 2026-10-15 10:07 CST

Previous Slurm array **63686231** has now completed all four model tasks with **COMPLETED (0:0)**. Its four model directories each contain the completion marker `CALIBRATION_COMPLETE`. This was independently checked in [GitHub Actions run 38669235325](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38669235325).

| Model | Old task | Slurm |
|---|---|---|
| M0_OH | 63686231_0 | COMPLETED (0:0), 00:16:17 |
| M1_OH | 63686231_1 | COMPLETED (0:0), 00:19:08 |
| M0_OH_H2O | 63686231_2 | COMPLETED (0:0), 00:18:29 |
| M1_OH_H2O | 63686231_3 | COMPLETED (0:0), 00:20:01 |

A later GitHub workflow [38668900795](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38668900795) inadvertently submitted a second array **63693436** for the same four scientific comparisons, in separate directory `~/Ca_OH_P0_calibration_v1/`. Since all four original calibrations were already completed, the later redundant array was cancelled without modifying original outputs. `sacct` confirmed **63693436_0–_3 all CANCELLED** after about 2 minutes 49 seconds; `squeue` was empty for the duplicate array.

**Use 63686231 as the canonical completed P0 calculation. Do not reuse/restart the redundant 63693436 array.** Actual corrected energy values and geometry force gates should be read from the four original `validation.json` files, rather than inferred from Slurm COMPLETED alone; no P1/P2 continuation was authorized or submitted in this check. The cancelled old `CaFresh` 23-job batch remains cancelled.

## Fixed-geometry P0 numerical outcomes — verified 2026-10-15 10:13:46 CST

This result was read live from the **canonical** `~/Ca_OH_dipole_calibration_20261015/*/validation.json` through read-only [GitHub Actions run 38669730002, attempt 2](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38669730002/attempts/2). The first attempt timed out connecting by SSH and did not retrieve results; retry succeeded. No new VASP or Slurm job was submitted during this status check.

| Model | Old E0 (eV) | New E0 (eV) | ΔE0 (eV) | New F / TOTEN (eV) | New max free-force (eV/Å) | Geometry gate <0.020 |
|---|---:|---:|---:|---:|---:|---|
| M0_OH | -146.71669462 | -146.71541142 | +0.00128320 | -146.75752250 | 0.00981093 | PASS |
| M1_OH | -147.63723064 | -147.67346118 | -0.03623054 | -147.70856646 | **0.02194710** | **FAIL: requires short consistent relaxation** |
| M0_OH_H2O | -161.39391351 | -161.39047832 | +0.00343519 | -161.43234408 | 0.01693181 | PASS |
| M1_OH_H2O | -162.53184536 | -162.55005424 | -0.01820888 | -162.58564174 | **0.06526722** | **FAIL: requires short consistent relaxation** |

Every model passed last-SCF EDIFF, VASP normal termination, expected original PAW/NELECT/NBANDS/effective IPSIN/LDIPOL/IDIPOL checks and `CALIBRATION_COMPLETE` marker. **This marker indicates electronically validated fixed-geometry Static, not that all geometries meet the force gate.** `all_force_gates_passed=false`. Each atom-coordinate output matches its archived Relax CONTCAR within ~8.3e-6 Å; no ionic relaxation occurred in this P0 batch.

The corrected reference-cancelled H2O incremental adsorption **double difference** is:

```text
DeltaDeltaE0(new) = -0.20152616 eV
DeltaDeltaF (new) = -0.20225370 eV
DeltaDeltaE0(old) = -0.21739583 eV
DeltaDeltaF (old) = -0.21814632 eV
Shift in E0 double difference = +0.01586967 eV
```

Because the same isolated H2O reference cancels between M0 and M1, this **double difference** is meaningful within the paired fixed-geometry calculations; individual terms `E(OH+H2O)-E(OH)` are **not absolute water adsorption energies without subtracting E(H2O)**. Ca still stabilizes the OH-covered water endpoint relative to pure Mg by ~0.202 eV, but **this remains a fixed-geometry estimate**. Neither faster/slower water dissociation nor corrosion resistance can be concluded from it.

All four runs used effective `ISPIN=1`, `ISMEAR=1`, `SIGMA=0.2`, `LDIPOL=T`, `IDIPOL=3`, proper fractional dipole centers, and 3x3x1 Gamma k mesh. Vacuum-charge cut ratios were ~1e-8 to 7e-8 of density peak, consistent with a very small electronic charge at the correction cut. This supports the numerical dipole setup but does **not** replace a vacuum/dipole-center convergence study.

**Next physically justified calculation:** under the same calibrated dipole settings, do **only two short ionic relaxations** beginning with already completed `M1_OH` and `M1_OH_H2O` calibrated geometries, then evaluate paired Static energies. Do not bypass the force gate. Before final quantitative comparison, make sure M0 controls remain consistent with the chosen relaxed references; optionally re-evaluate the paired M0 statics if numerical settings beyond P0 change. The previously cancelled 23 `CaFresh` tasks remain cancelled. No P1/P2 jobs were automatically submitted.
