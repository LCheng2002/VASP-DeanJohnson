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
