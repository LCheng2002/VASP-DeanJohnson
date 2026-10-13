# Mg₂NiHₓ 低氢小体相筛选：Slurm 提交记录

**2026-10-13 23:07 CST（北京时间）成功提交并核验 Slurm 回执。**

- 作业数组：**63638134**；子任务索引 **0–4**；共 **5** 个 Relax→Static 模型；无人工并发上限。
- GitHub Actions 提交工作流：[run 38419250472](https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38419250472)，结论 SUCCESS。
- 提交目录：`~/Mg2Ni_lowH_bulk_20261013/`。
- 初次 `squeue` 为 **PENDING (Priority)**，`sacct` 为 `63638134_[0-4]|PENDING`。**目前不能宣称 VASP 收敛。**

| 子任务 | 模型 | 组成 | 初始四面体位点 |
|---|---|---|---|
| 63638134_0 | M0_Mg2Ni | Mg24Ni12 | 0H 参考 |
| 63638134_1 | M2_T31_T31 | Mg24Ni12H2 | 两个 Mg3Ni 型 |
| 63638134_2 | M2_T31_T22 | Mg24Ni12H2 | Mg3Ni + Mg2Ni2 型 |
| 63638134_3 | M3_T31_T31_T22 | Mg24Ni12H3 | 两个 Mg3Ni + 一个 Mg2Ni2 型 |
| 63638134_4 | M3_T31_T22_T22 | Mg24Ni12H3 | 一个 Mg3Ni + 两个 Mg2Ni2 型 |

生产前进行了以下真实核验：

- 旧 Mg₂Ni 18 原子 Relax+Static 参考：`02_static/POSCAR` SHA-256 `8fc1354eadf5590cd5ce6a84a9ba4eb321779eadc4b7e6c29e28433ba1220be1`，原参考最大力 **0.016505 eV/Å**。
- 已收敛气相 H₂：**F = −6.77240668 eV**，沿用同一 H 赝势。
- 许可 Mg/Ni/H POTCAR 在 HPC 本地组装，完整合并势 SHA-256 `72b607bd5434dcb71996bb13fd478af3592115762200d3c53e0f47c68b0a84ae`；**没有上传赝势文件到 GitHub**。
- 4 个含 H 模型的最短初始 H–H 距离均 **>4.90 Å**，局部金属四面体配位类型通过几何检查。
- 不重提旧两个 147 原子 MgH₂/Mg₂NiH₀.₂₅ 界面模型，也不影响已提交的 LaH₃、纯 Mg 和 H₂ 解离 NEB 任务。

数值设置：PBE、ENCUT=520 eV、ISPIN=2、LDAU=false、5×5×2 K 点、全晶胞体相 ISIF=3，EDIFF=1E-7，优化 EDIFFG=−0.015 eV/Å；32 MPI、128G、12h/子任务。Relax 和 Static 必须分别通过电子、离子、受力等数值验证，才标记 `VERIFIED_FINISHED`。

复核状态（只读，不重复提交）：

```bash
python3 ~/Mg2Ni_lowH_bulk_20261013/screen.py --status
```

本批候选 H 配位仅用于 **x≈0.167、0.25 的有限占位热力学筛选**，不能自动宣称已证明稳定 Mg₂NiH₀.₂₅ 物相。详细物理局限见 [README.md](README.md)，生成规则和固定参考见 [plan.json](plan.json)。
