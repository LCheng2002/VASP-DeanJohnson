# B组续算提交与队列快照

核查时间：2026-10-11 11:34:07 CST（北京时间）。

| 原Job ID | 模型 | 新Job ID | 新名称 | 本次核查状态 |
|---|---|---|---|---|
| 63535597 | B-H-Mg_Mg_Mg2NiH4 | 63554410 | mgni_Bretry_Mg | PENDING (Priority) |
| 63535598 | B-H-Ni_Mg_Mg2NiH4 | 63554411 | mgni_Bretry_Ni | PENDING (Priority) |

两个任务于2026-10-11 11:33提交，已收到真实Slurm Job ID。11:34时尚未获得计算节点，因此不能把GitHub工作流成功当作VASP续算完成。

新根目录：`~/Mg2Ni_Mg_H_thermo_validation_20261010/retry_20261011_B_v1/`。
每项1节点、64 MPI ranks、24小时上限；从原最后CONTCAR继续，POTIM=0.10、NSW=500、EDIFFG=-0.02。其余物理设置保留，Relax达标后才运行Static，原任务所有文件保持不变。

原始提交证据：
https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38108731427

原始Slurm队列核查证据：
https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38108796020

本文件记录上述时刻的快照，不是自动更新的实时状态。
