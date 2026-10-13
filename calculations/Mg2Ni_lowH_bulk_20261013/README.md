# Mg₂NiHₓ 低氢占位筛选：小体相模型（2026-10-13）

## 目标与科学问题

本批是对之前两个 147 原子 `MgH₂/Mg₂NiH₀.₂₅` 界面结构优化超时的**独立小模型筛选**，并不恢复它们的旧任务，也不直接模拟 MgH₂ 向 Mg₂NiHₓ 的跨界面传 H。

**目的：** 从已收敛的 Mg₂Ni 金属晶胞出发，以 36–39 原子的体相模型找出少量 H 的合理局部间隙配置。比较 `x=0.167` 与 `x=0.25` 及同浓度下的两种四面体占位；待本批收敛后再决定是否建立 MgH₂/Mg₂NiHₓ 界面。

## 结构来源与模型

- 已完成 `Mg₂Ni` 体相 Relax/Static：`~/Mg2Ni_H_transfer_bulk_screen_20261009_supplement_v2/supplement_20261009/02_Mg2Ni/02_static/POSCAR`。
- 固定该源 POSCAR SHA-256：`8fc1354eadf5590cd5ce6a84a9ba4eb321779eadc4b7e6c29e28433ba1220be1`；原 18 原子为 `Mg₁₂Ni₆`，沿 `c` 轴扩展 2 倍成为 `Mg₂₄Ni₁₂`。
- 初始 H 放置在经 Mg/Ni 近邻数验证的四面体空隙中心，`T31` 指 H@Mg₃Ni，`T22` 指 H@Mg₂Ni₂。具体金属-金属空隙中心和周期坐标见 `plan.json`。
- **这些是与 Oshchapovsky 等文献 X3/X7 局部配位类型相似的几何候选位点，不是复制了论文 X3/X7 的精确实验或计算分数坐标。** 是否为稳定占位，以后续弛豫为准。

| 模型 | 组成 | 总原子数 | 名义 x | 候选 H 环境 |
|---|---|---:|---:|---|
| M0_Mg2Ni | Mg24Ni12 | 36 | 0 | 洁净 Mg₂Ni 体相对照 |
| M2_T31_T31 | Mg24Ni12H2 | 38 | 0.167 | 两个 Mg₃Ni 型间隙 |
| M2_T31_T22 | Mg24Ni12H2 | 38 | 0.167 | Mg₃Ni + Mg₂Ni₂ 型 |
| M3_T31_T31_T22 | Mg24Ni12H3 | 39 | 0.25 | 2 个 Mg₃Ni + 1 个 Mg₂Ni₂ 型 |
| M3_T31_T22_T22 | Mg24Ni12H3 | 39 | 0.25 | 1 个 Mg₃Ni + 2 个 Mg₂Ni₂ 型 |

初始构型按同一个 36 原子母体生成，两个/三个 H 相互间距均至少 **4.8 Å**，不预先形成 H₂。只有 `x=0.167` 同成分两个构型、`x=0.25` 同成分两个构型之间的裸总能量可直接比较。总能量不能跨不同 H 数直接排序。

## 数值设置（全部任务统一）

- VASP 6.3.2，PBE，`ENCUT=520 eV`，`ISPIN=2`，`LDAU=.FALSE.`。
- `ISMEAR=0; SIGMA=0.05; LREAL=.FALSE.; LASPH=.TRUE.; ADDGRID=.TRUE.`。
- 体相结构允许原子与晶胞松弛：`ISIF=3`，`IBRION=2`，`NSW=240`，`POTIM=0.12`。
- `EDIFF=1E-7`，优化 `EDIFFG=-0.015 eV/Å`，Static 受力重新检查 `<0.020 eV/Å`；`KPOINTS=Γ 5×5×2`。
- 初始自旋：Mg 和 H 为零，Ni 为 1.0 μB/原子。未做各类磁序全面搜索，不能将结果与旧非自旋极化模型的裸总能量混比。
- 每个模型独立 Relax→Static；32 MPI 进程，内存 128G，12h 上限，一个 Slurm 数组（5 个子任务，无并发上限）。
- **POTCAR 仅在思源一号授权用户目录读取和组装**，根据精确 SHA 校验 Mg/Ni/H PAW；不会提交到 GitHub。

## H₂ 参考与能量计算

复用此前已经收敛的同一 PBE/H PAW 气相 H₂：
`~/MgNi_LaH3_extensions_20261012/H2_reference`，`F(H₂)=-6.77240668` eV。生产前检查 H2 完成标记、电子与受力收敛、与 Mg/Ni/H 合并势中的 H POTCAR 完全一致。

每个含 n 个 H 的结构获得稳定终态后，可计算电子能层面的低氢插入：

\[
\Delta E_{\mathrm{ins}}(n)=E(\mathrm{Mg}_{24}\mathrm{Ni}_{12}\mathrm H_n)
-E(\mathrm{Mg}_{24}\mathrm{Ni}_{12})-\frac{n}{2}E(\mathrm H_2).
\]

还可在相同底层数据下评估从 2H→3H 加入第三个 H 的增量反应能：

\[
\Delta E_{2\rightarrow3}=E(\mathrm{Mg}_{24}\mathrm{Ni}_{12}\mathrm H_3)
-E(\mathrm{Mg}_{24}\mathrm{Ni}_{12}\mathrm H_2)-\frac12 E(\mathrm H_2).
\]

后者应报告使用了哪一种 2H/3H 结构以及 H 的最终占位。若多个模型最终落入不同局部极小值，不应在未做更多采样前声称全局基态。

**局限：** 5×5×2 网格和 2×c 小超胞是第一轮筛选，不代表 K 点、体相有限尺寸、构型熵已经收敛。低氢 Mg₂NiHₓ 是可能存在部分占位和无序的固溶体；有限有序超胞不能被直接认作唯一独立相。后续仍需要结合文献坐标、实验与相稳定性判断。

## 提交与安全边界

- 本次科学范围固定 5 个独立模型，不做大型界面续算、NEB、Mg₂NiH₂ 数据库独立晶型或自动穷举。
- 入口程序：`screen.py`，源晶胞与候选位点见 `plan.json`。
- GitHub 结构生成预览：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38418881569 ，已检查 5 个模型组成、四面体金属近邻及 H–H 避碰。这个工作流 **没有运行 VASP 或 sbatch**。
- 正式提交后，生产目录固定 `~/Mg2Ni_lowH_bulk_20261013/`，以 `submission.json` 中真实 Job ID 为准。若已提交，程序禁止重复提交；未收敛则保留诊断文件，不放宽验收、不自动增加新模型。
- Slurm 状态可在集群执行：`python3 ~/Mg2Ni_lowH_bulk_20261013/screen.py --status`（只读）。

参考文献：Oshchapovsky et al., *Materials Today Communications* 27 (2021) 102174, DOI [10.1016/j.mtcomm.2021.102174](https://doi.org/10.1016/j.mtcomm.2021.102174)。作者相关原始资料 DOI [10.17632/dn7nng9b45.1](https://doi.org/10.17632/dn7nng9b45.1)。
