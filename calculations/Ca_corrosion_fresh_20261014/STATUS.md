# Ca抗腐蚀重新建模：真实提交记录

**已通过GitHub提交思源一号。** 23次sbatch调用均有真实回执。开始提交：2026-10-14 22:18:39 CST；全部提交完成：22:19:02；Slurm核查：**2026-10-14 22:19:08 CST（北京时间）**。

## 本次任务清单

| 模型 | Slurm Job ID | 核查时状态 |
|---|---|---|
| 新Mg体相参考 Mg_bulk | 63681995 | PENDING (Priority) |
| 新H2参考 | 63681996 | PENDING (Priority) |
| 新H2O参考 | 63681997 | PENDING (Priority) |
| M0_clean | 63682000 | PENDING (Dependency) |
| M0_water | 63682002 | PENDING (Dependency) |
| M0_OH | 63682003 | PENDING (Dependency) |
| M0_OH_H | 63682004 | PENDING (Dependency) |
| M0_OH_water | 63682005 | PENDING (Dependency) |
| M0_2OH_H | **63682006** | PENDING (Dependency) |
| M0_2OH_near | 63682008 | PENDING (Dependency) |
| M0_2OH_far | 63682009 | PENDING (Dependency) |
| M0_4OH_near | 63682010 | PENDING (Dependency) |
| M0_4OH_far | 63682011 | PENDING (Dependency) |
| M1_clean | 63682012 | PENDING (Dependency) |
| M1_water | 63682013 | PENDING (Dependency) |
| M1_OH | 63682014 | PENDING (Dependency) |
| M1_OH_H | 63682015 | PENDING (Dependency) |
| M1_OH_water | 63682016 | PENDING (Dependency) |
| M1_2OH_H | **63682021** | PENDING (Dependency) |
| M1_2OH_near | 63682022 | PENDING (Dependency) |
| M1_2OH_far | 63682023 | PENDING (Dependency) |
| M1_4OH_near | 63682024 | PENDING (Dependency) |
| M1_4OH_far | 63682026 | PENDING (Dependency) |

所有23项在上述核查时均为NO_FINAL_VALIDATION，**不代表VASP已经收敛**。3个新体相/气相初始结构已在集群生成；其余20个生产结构按各自新母体收敛结果生成，不使用旧结构。23份理想初始几何已在本地和GitHub独立生成、验碰并逐个核对POSCAR哈希；它们是预览，不是收敛结果。

## 本批与旧计算的区别

- 从理想hcp-Mg晶格和解析H2/H2O/OH几何重新建模；不导入旧POSCAR、CONTCAR、能量或任何已失败任务的依赖关系。
- Mg初始a=3.21、c=5.21 Å仅作初猜，先完成新Mg体相优化，再构造4×4×6层的Mg(0001)表面。
- M0为Mg96，M1为Mg95Ca1，Ca替换一个表层Mg；底部2层共32原子固定，真空20 Å。
- 23模型=13个重新建立的参照/基线+原计划10个第二分子水解离及OH覆盖构型。
- 每个任务只依赖其直接需要的新母体；不使用整体数组阻塞无关分支。
- 表面PBE、520 eV、ISPIN=2、无U、Gamma 3×3×1、ISMEAR=0/SIGMA=0.05、LREAL=.FALSE.，LDIPOL/IDIPOL=3；优化/静态自由原子力阈值分别0.015/0.020 eV/Å。
- 表面每任务64 MPI、128GB、24h；不自动重试，不增加NEB或AIMD，不声称已形成保护膜。

## 赝势问题已经解决

第一次尝试仅找到默认Mg赝势，在sbatch之前停止，没有新作业。第二次在集群的已有授权POTCAR中找到Mg、Ca_sv、O、H（仅检查5份POTCAR、8个目录），在集群本地为本批统一组装。**没有读取或导出旧结构与旧能量，也没有上传POTCAR或凭据到GitHub。**

生产目录：`~/Ca_corrosion_fresh_20261014/`。
提交回执：该目录下的`submission.json`。
查询命令（只读，不重提）：`python3 ~/Ca_corrosion_fresh_20261014/fresh.py --status`。

- 输入验证：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38578586431
- 成功提交及squeue/sacct核查：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38579248394
- 提交日志artifact：Ca-fresh-submission-status，11902563032。

本文件记录核查时刻快照，不会自行更新；没有修改其他Ca旧项目、MgNi、LaH3或NEB任务。
