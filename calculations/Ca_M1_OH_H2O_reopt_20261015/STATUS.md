# M1_OH_H2O 单独重弛豫：真实 Slurm 提交记录

2026-10-15 10:24:16 CST (UTC+08) 已通过 [GitHub Actions run 38670721598](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38670721598) 真正提交到思源一号；**只有一个新 Slurm Job ID：63693676**。

| 任务 | Job ID | 工作目录 | 首次 Slurm 状态 |
|---|---|---|---|
| M1_OH_H2O，先结构 Relax 再 Static | **63693676** | `~/Ca_M1_OH_H2O_reopt_20261015/` | **PENDING (Priority)**，2026-10-15 10:24:21 CST |

提交的校验流程已实际通过：
- P0 来源 `63686231_3` 的 `CALIBRATION_COMPLETE`、实际旧输入哈希和有效 PBE/PAW、ISPIN=1、ISMEAR=1、SIGMA=0.2、NBANDS=160、LDIPOL=T、IDIPOL=3 等；
- 母体 `F=-162.58564174 eV`、`E0=-162.55005424 eV`，且最大可移动原子受力 **0.06526721594 eV/Å**：确有几何重弛豫必要；
- 直接在集群内复用母体已授权 POTCAR，组合 SHA256 `2f62cba23027f134e5e5c49961c085dd38eededc9046d7bebf475fbded911d01`；POTCAR 不经 GitHub 上传；
- 原有 101 原子、Mg95CaO2H3、底部 32 原子固定，晶胞固定、k网格和 DIPOL 沿用校准值；
- 仅产生一个独立重算任务，`Relax (EDIFFG=-0.012, NSW=160)` → `Static`，受力严格验收 <0.015/<0.020 eV/Å；
- 已上传日志 Artifact：`CaM1_OH_H2O_reopt_submission_20261015` (GitHub artifact ID `11943537380`)；
- 首次 `squeue` 和 `sacct` 均确认 Job ID `63693676` 为 PENDING，**尚无 VASP 收敛结果**。

后续核查（只读，不再次提交）：

```bash
python3 ~/Ca_M1_OH_H2O_reopt_20261015/reopt.py --status
```

若 Relax 或 Static 的电子、自洽、晶胞/受限原子或自由原子力未达到既定条件，脚本拒绝写入整体完成标记，不会放宽验收、也不自动重试。收敛后应进一步核查新水是否仍然保持两个 O–H 键；终态解离或发生 H 转移时需重新分类。

**不要将本任务的新总能量直接与旧 M1_OH 或 M0/M1 对照的裸能量混用。** 本轮只重弛豫 `M1_OH_H2O`，要得到完全一致的 Ca 对照吸附能差，还应处理未通过几何受力阈值的 `M1_OH`，但该任务本次未经授权、未提交。

此前原始 P0 `63686231` 保留，重复 P0 `63693436` 已取消；23 项 `CaFresh` 旧批次仍处于取消状态。本次没有涉及 MgNi、LaH3、H2-NEB 其他项目。
