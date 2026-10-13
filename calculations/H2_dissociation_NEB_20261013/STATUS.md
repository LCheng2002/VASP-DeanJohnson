# H2解离CI-NEB：提交与核验快照

实际提交：**2026-10-13 11:41:36 CST（北京时间）**。
最近本轮查询：**2026-10-13 11:43:15 CST**。这是状态快照，不会自行更新。

| Slurm Job ID | 任务 | 11:43:15 实际状态 |
|---|---|---|
| 63614945_0 | 纯Mg洁净参考，Relax→Static | PENDING (Priority) |
| 63614945_1 | Mg/Mg2Ni洁净参考，Relax→Static | PENDING (Priority) |
| 63614946 | H2参考，Relax→Static | **RUNNING**，node041，已运行21秒 |
| 63614947_0 | 纯Mg分子态/双H端点 | PENDING (Dependency) |
| 63614947_1 | Mg2Ni暴露区域分子态/双H端点 | PENDING (Dependency) |
| **63614948** | **纯Mg：NEB→CI-NEB** | **PENDING (Dependency)** |
| **63614949** | **Mg2Ni暴露区域：NEB→CI-NEB** | **PENDING (Dependency)** |

五次sbatch调用均有真实回执（其中两项为数组）；Slurm sacct/squeue均已核验。**提交成功不等于能垒已经算出**，两条路径目前均报告NOT_YET_VERIFIED。

- 提交工作流：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38342837099
- 后续只读核查：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38343023015
- 生产目录：`~/H2_dissociation_compare_20261013/`
- 程序和详细模型说明：本目录的README.md、plan.json、neb_core.py、neb_driver.py。

## 模型和参数摘要

纯Mg体系为Mg(10-10)有限薄片，24个基底原子，含H端点Mg24H2（26原子）。含Mg2Ni体系为现有Mg(10-10)/Mg2Ni(100)界面暴露的Mg2Ni区域，36Mg+6Ni，含H端点Mg36Ni6H2（44原子）；**不是独立单相Mg2Ni薄片或纯Ni表面**。

两个基底同面内尺寸约12.95437×5.18097 Å、每单元1个H2、相同3×6×1 K网格。均固定最底部8个Mg原子，其余自由；重新计算洁净表面与所有端点。共同物理设置：PBE、ENCUT=520、ISPIN=2、无U、ISMEAR=0/SIGMA=0.05、LREAL=.FALSE.、EDIFF=1e-7、LDIPOL/IDIPOL=3。端点力要求<0.01 eV/Å。

每条路径5个中间图像，加2个固定端点；NEB先到0.05 eV/Å，再开启LCLIMB到0.03 eV/Å。VTST FIRE (IOPT=7)，SPRING=-5。每条带使用一个节点60 MPI进程，分成5×12，各图像NCORE=4。每条NEB作业24小时上限；不会无限续算。

集群实际使用了VTST编译版本：`vasp/6.3.2-vtst-intel-2021.4.0`。实际二进制SHA-256：`454da6713051a916308a37843d23fb41cd31161285dcf09606a830a9f470a106`，已确认包含LCLIMB、ICHAIN、IOPT标识，后续计算还要求输出相应NEB标记。

## 文件核验

已回收23个白名单输入/驱动/回执文件，并按生产manifest逐项比对SHA-256（manifest、submission为单独记录文件）。包含3份实际生产起始POSCAR、4份只作来源记录的既有CONTCAR转换文件及控制程序；不含完整POTCAR或凭据。

使用回收的实际生产版本程序进行离线测试，两条预览路径均通过组成、H标签、周期边界、固定原子、IDPP与避碰检查。该测试没有执行VASP或sbatch。**正式4个端点与两条路径的生产图像仍需在母体/端点收敛后构建，不将预览当作收敛结构。**

最终将报告相对于分子吸附态及洁净表面+气相H2的两种能量参照，并给出反应能和逆向势垒。不存在能量峰或有效力检查失败时不报告为零势垒。尚未开展过渡态虚频、双向连接、零点能、温压或扩散计算；未修改原有LaH3、低氢Mg2NiHx、Mg17La2作业。
