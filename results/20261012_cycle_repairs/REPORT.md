# MgNi 吸放氢任务：未通过的含义、诊断与修复提交

诊断采集时间：2026-10-12 00:40:00 CST。修复提交时间：00:49:14 CST。状态核查：00:51:02 CST。以上均为北京时间。

## 一、未通过不等于所有VASP任务均崩溃

本次三个失败的计算中，两项属于“VASP已经正常算完，但脚本的最终静态受力检查未通过”，一项属于真正的内存不足终止。另一个回填作业因依赖失败取消，未运行VASP。

| 原作业 | 模型 | 确认原因 |
|---|---|---|
| 63565364_1 | D_S1b_split_Ni | Relax正常收敛；Static正常结束且最终电子自洽收敛；Static最大自由原子力0.0203883922 eV/Å，超过0.02验收线 |
| 63565364_5 | A_H2_Mg_parallel | Relax正常收敛；Static正常结束且最终电子自洽收敛；Static最大自由原子力0.0200247159 eV/Å，超过0.02验收线 |
| 63565365 | H2_reference | VASP初始化阶段内存不足，vasp.stderr明确记录oom_kill和Out Of Memory；没有有效电子迭代或收敛结构 |
| 63565366 | S2/S3回填链 | 依赖的63565364_1失败，按kill-on-invalid-dep配置自动取消；S2/S3均未构建 |

两项固体的Slurm FAILED/1:0来自计算后的Python质量检查，不是VASP自身运行崩溃。原能量作为待复核中间记录保留，不进入正式机理结论。0.02 eV/Å是本批自行设定的验收要求，不是“超过就发生程序错误”的阈值。Static的ionic_converged=false在NSW=0时正常，不是本次失败原因。

## 二、固体为何Relax过关、Static仍被阻断

| 模型 | Relax最大力 eV/Å | Static最大力 eV/Å | 比0.02超出 |
|---|---:|---:|---:|
| D_S1b_split_Ni | 0.0197545534 | 0.0203883922 | 1.942% |
| A_H2_Mg_parallel | 0.0191497434 | 0.0200247159 | 0.124% |

原Relax采用LREAL=Auto、EDIFF=1e-6，Static采用LREAL=.FALSE.、EDIFF=1e-7。Static在同一末态几何上以另一数值精度重新评估受力，两个模型的结果刚好跨过验收线。直接故障原因已经由错误日志的“Convergence check failed .../02_static”确认；数值设置的变化是需要针对的原因，但本次没有通过额外对照把影响全部单独归因于LREAL。

不是NSW跑满：两个原Relax分别83步和110步，远小于400步上限。两份原Static的vasp.stderr为空，OUTCAR有正常结束页脚和最终EDIFF收敛标记。

实际修复：从各自02_static/CONTCAR继续，在独立目录进行精修；Relax改为LREAL=.FALSE.、EDIFF=1e-7，与原Static一致；EDIFFG=-0.015，给最终0.02验收线留余量；NSW=200、POTIM=0.10。之后仍用原Static参数核验。泛函、赝势、ENCUT、KPOINTS、自旋、U、晶胞和原子组成均不改变。没有将验收线放宽到0.03或0.05。

## 三、H2为什么内存不足及修复

原H2为20 Å立方盒、ENCUT=520、4个MPI进程，但只申请总内存4 GB。日志明确写：

```
Detected 1 oom_kill event
Out Of Memory
```

OUTCAR在初始化时报告rank0约1512798 kBytes内存，其中FFT计划约999584 kBytes、网格约460800 kBytes。虽然只有两个原子，较大盒子和密集平面波网格仍占内存；原4 GB申请偏小。不能用随后出现的MPI启动超时表象取代明确的OOM证据。

修复只把Slurm总内存4G提高到16G，保留4个MPI进程、20 Å盒子、截断能、K点、H赝势和电子设置；保持H2的0.005 eV/Å目标。实际读取新Slurm脚本和sacct已经确认ReqMem=16G。未缩小盒子或降低ENCUT来绕过内存问题。

## 四、已经提交的新任务与实际状态

| 新Job ID | 修复内容 | 00:51:02核查状态 | 资源 |
|---|---|---|---|
| 63575040 | D_S1b_split_Ni精修→Static | PENDING (Priority) | 1节点、64 MPI进程、12小时上限 |
| 63575041 | A_H2_Mg_parallel精修→Static | PENDING (Priority) | 1节点、64 MPI进程、12小时上限 |
| 63575042 | H2内存修复后重算 | RUNNING，node195，已运行1分16秒 | 4 MPI进程、16G、30分钟上限 |
| 63575043 | 恢复S2/S3依赖链 | PENDING (Dependency) | 1节点、64 MPI进程、24小时上限 |

新回填依赖是afterok:63575040。已成功的S1a只复制为参考，不重算。S1b修复后，仍须核验两者并比较终态能量，再按原原子守恒、供H位置及避碰规则生成S2，S2通过后再生成S3。未修改物理问题、未新增NEB、未直接跳过前置验证。

00:51时H2处于运行，尚没有本轮通过验收的最终能量；两个精修仍排队，回填尚未生成。提交恢复不等于计算已经全部完成或保证收敛。

新目录：`~/MgNi_abs_des_cycle_20261011/repairs_20261012_v1/`。

原任务的输出、失败日志、执行锁与提交记录均保留；没有覆盖原输入、没有重复提交已成功模型，也未取消其他正在运行的任务。新建模输入和源CONTCAR/POTCAR均校验哈希；新目录、提交锁和逐项作业回执防重复。

## 五、追溯

- 原因核验：GitHub Actions run 38156380933，artifact cycle-failure-audit-20261012。
- 修复与提交：GitHub Actions run 38156995500，artifact cycle-repair-submission-20261012。
- 实际状态与参数：GitHub Actions run 38157157195，artifact cycle-repair-status-20261012。
- 后续只读检查：`.github/workflows/check-cycle-repairs-20261012.yml`。

上述报告是本次状态快照，不会自行更新。旧作业仍保留FAILED/CANCELLED历史记录，应使用上面的新Job ID追踪修复进度。

方法依据：VASP官方LREAL文档 https://vasp.at/wiki/index.php/LREAL 及EDIFFG文档 https://vasp.at/wiki/index.php/EDIFFG。用户计算的具体失败原因来自本次集群实际日志，不由这些通用文档推断。
