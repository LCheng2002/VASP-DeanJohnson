# M1_OH_H₂O 单模型重弛豫（2026-10-15）

## 用户指定范围

**仅重新弛豫 M1_OH_H₂O：Mg₉₅CaO₂H₃，101 原子。**

已完成且经核验的 P0 固定几何偶极校准 `63686231_3`：Static `F=-162.58564174 eV`，`energy(sigma→0)=-162.55005424 eV`，最大可移动原子受力 **0.0652672159 eV/Å**。电子收敛且正常结束，但旧几何未通过 0.02 eV/Å 的最终受力要求。来源：[已完成四点校准状态](../Ca_OH_dipole_calibration_20261015/STATUS.md)。

## 严格的科学与输入约束

- 从 `~/Ca_OH_dipole_calibration_20261015/M1_OH_H2O/POSCAR` 开始，不使用后来取消的 63693436，也不复算其他 M0/M1 结构。
- 与 P0 保持原 PBE、ENCUT=520 eV、ISPIN=1、ISMEAR=1、SIGMA=0.2、Γ 3×3×1 K 点、NBANDS=160、`LREAL=.FALSE.`、`LDIPOL=.TRUE.`、`IDIPOL=3` 及 **原始校准 DIPOL 分数坐标**。
- 按原模型保留 Mg95CaO2H3 原子顺序、底部 32 个固定金属原子、真空层及固定晶胞。 `ISIF=2` 不放开晶胞。
- **Relax** `IBRION=2, NSW=160, POTIM=0.10, EDIFFG=-0.012 eV/Å`，数值验收最大全自由原子力 `<0.015 eV/Å`。
- **Static** 从新 Relax 的收敛 CONTCAR 出发，保持相同物理参数、`EDIFF=1E-7, LREAL=.FALSE.`，受力验收 `<0.020 eV/Å`。
- 在新旧 M1_H₂O 相同理论与赝势参数下比较 Static 电子能量；检查水是否仍保留两条 O–H 键，及是否出现质子交换或解离。若优化后水已解离，应按终态重新分类。
- 真实反应能和 Ca 对照必须配对经过一致处理的 M1_OH 和 M0 对照。**本轮不会自动提交这些附加模型。**

## 安全、防重复与执行

- 原样复用 P0 集群已授权 POTCAR，精确总 SHA-256 `2f62cba23027f134e5e5c49961c085dd38eededc9046d7bebf475fbded911d01`，仅在思源一号节点内复制；不上传 POTCAR、凭据、波函数或电荷文件到 GitHub。
- 提交前校验原 `submission.json` 的 `63686231` Job ID、完成标记、真实收敛校验 JSON、P0 原输入哈希、原始组成及力；如任一不符则停止。
- 独立生产目录：`~/Ca_M1_OH_H2O_reopt_20261015/`。已有提交回执时绝不再提交；运行前检查同名活跃 Slurm 作业，避免重复计算。
- 单 Slurm 作业，64 MPI、128G、12 小时，`Relax→Static` 顺序执行。**没有 NEB，没有其他 CaFresh 任务重新激活。**
- 旧 23 个 CaFresh 批次、此前已取消的 `63693436`，以及 Mg–Ni/LaH₃ 其他项目的 Slurm 作业不受影响。

真实提交的 Slurm Job ID 以 `submission.json` 和 GitHub Actions 记录为准。

状态查询（只读）：

```bash
python3 ~/Ca_M1_OH_H2O_reopt_20261015/reopt.py --status
```

结果通过两个阶段检查后写入 `~/Ca_M1_OH_H2O_reopt_20261015/final_result.json`。
