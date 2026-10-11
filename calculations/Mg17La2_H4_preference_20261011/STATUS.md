# H4实际提交与队列核查

提交时间：2026-10-11 12:00:33 CST（北京时间）。
队列核查时间：2026-10-11 12:01:44 CST。

| Slurm数组任务 | 模型 | 核查状态 |
|---|---|---|
| 63554980_0 | H4_La_partner | PENDING (Priority) |
| 63554980_1 | H4_Mg_top | PENDING (Priority) |
| 63554980_2 | H4_Mg_bridge | PENDING (Priority) |

Slurm接受的数组为 `63554980_[0-2%2]`，三个模型，每个1节点64 MPI进程、12小时上限，最多同时运行2个。该状态是快照，不是自动更新。

实际执行：在思源一号读取并核验H3_addLa的CONTCAR/OUTCAR/INCAR/KPOINTS及POTCAR校验值，生成3个Mg68La8H4模型。前79原子和晶胞一致，只有新增H初始位置不同。三个候选全部通过避碰和组成检查；新增H与既有H的最短周期距离依次为3.463732、4.061716、3.309592 Å。

已提交不等于计算完成；上述时刻三个任务尚未开始VASP。完成后需要分别核验Relax/Static的最后电子收敛、力标准及最终4个H的配位。

真实目录：`~/TaskA_Mg17La2_H4_preference_20261011/`。

提交记录：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38110157617

队列核查记录：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38110263239

初期受阻的是额外的母结构导出工作流，未执行该导出。最终使用集群内直接建模/提交方式成功完成，无需重新上传母结构或POTCAR。
