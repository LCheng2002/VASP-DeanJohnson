# 纯Mg与Mg2Ni暴露区域：H2解离能垒对比（CI-NEB）

## 1. 实际执行状态

**已于2026-10-13 11:41:36（北京时间）通过GitHub Actions提交思源一号。** 11:41:46核查时，参考态任务排队，端点与两条NEB路径等待依赖；当前没有新的解离能垒数值。

| Slurm Job ID | 实际任务 | 条件/当时状态 |
|---|---|---|
| 63614945_0、63614945_1 | 纯Mg与Mg/Mg2Ni洁净表面的统一参数重优化、Static | PENDING (Priority) |
| 63614946 | H2气相参考的Relax→Static | PENDING (Priority) |
| 63614947_0、63614947_1 | 每个表面的分子态与解离态端点，分别Relax→Static | PENDING (Dependency)，等待两洁净表面 |
| **63614948** | **纯Mg：5-image NEB→CI-NEB** | 等待纯Mg端点和H2参考通过 |
| **63614949** | **Mg2Ni暴露区域：5-image NEB→CI-NEB** | 等待对应端点和H2参考通过 |

- 提交证据：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38342837099
- 只读核查工作流：`.github/workflows/check-h2-neb-20261013.yml`
- 原始端点与环境审计：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38341672691
- 生产目录：`~/H2_dissociation_compare_20261013/`
- 驱动部署目录：`~/VASP/github_jobs/H2_NEB_driver_20261013/`

已得到五项sbatch回执，其中两个为数组，因此不等于只运行五个几何结构。实际科学范围为两条路径、两份洁净表面参考、一份H2参考、四个端点；每条路径5个运动图像。尚未收敛的端点不会被强行插值。原有LaH3、Mg2NiHx、Mg17La2和其他计算均未修改、取消或重提。

## 2. 比较的准确对象

| 系统 | 母体 | 含H端点 | 本轮定位 |
|---|---|---|---|
| Mg | 已收敛的独立Mg(10-10)有限薄片，24 Mg | Mg24H2，26原子 | 真正不含Ni的表面对照 |
| Mg2Ni_region | 既有Mg(10-10)/Mg2Ni(100)模型，36 Mg+6 Ni | Mg36Ni6H2，44原子 | Mg支撑的暴露Mg2Ni区域，**不是纯Ni、也不是独立单相Mg2Ni薄片** |

这是延续用户既有吸附结果的第一轮局部路径比较。可判断选定Mg2Ni暴露区域相对于选定纯Mg表面的能垒差异，不应扩写为所有晶面、所有Mg2Ni材料的普遍结论，也不能单独分离化学成分、应变及支撑效应。

两系统采用相同的面内晶胞（约12.95437×5.18097 Å）、相同每表面单元一分子H2的覆盖度和同一3×6×1 K网格。原有限薄片几何厚度不同：纯Mg约6.166 Å，Mg/Mg2Ni约11.169 Å；相应名义真空约21.048和16.045 Å。尚未完成厚度或真空收敛研究。

## 3. 结构来源与初终态

所有四个来源CONTCAR均在11:26的只读审计中确认存在，最后电子自洽和正常结束标记齐全；具体路径与SHA-256记录在plan.json。

- 纯Mg母体：此前机理补充批次已完成的Mg_parent/02_static/CONTCAR；其吸附分支尚未开始，不将尚不存在的Mg_H2/Mg_2H当作端点。
- Ni侧洁净母体：A0_Mg_Mg2Ni_clean/02_static/CONTCAR。
- Ni侧分子候选：A_H2_Ni_vertical的已收敛CONTCAR；H-H约0.753 Å。
- Ni侧产物候选：A_H2_Ni_parallel的已收敛解离终态；H-H约2.269 Å。

为避免自由薄片整体漂移，两系统均居中并固定底部8个Mg原子，其余金属和H放开；初终态的固定位置必须完全一致。新增LDIPOL/IDIPOL=3、统一自旋及更密K网格后，旧能量不能直接沿用，因此两份洁净表面及全部端点均重新优化并计算Static。

Ni侧端点在新洁净母体上，按原Ni41局域坐标搬移既有H对，随后完全优化可移动原子。纯Mg分子候选在暴露Mg邻近平行放置0.75 Å的H2，双H候选预分开至2.00 Å，分别优化。若分子解离/逸出、产物复合，或端点发生过大重构，停止对应NEB链；不能将未解离或未稳定的构型强行称为目标端点。

**Ni侧路径可能包含H2转动、靠近表面和断键，不是人为固定取向的纯键长扫描。** 只有一条局部路径，不是全局最小势垒搜索。

## 4. 同一计算口径

| 参数 | 设置 |
|---|---|
| 程序 | `vasp/6.3.2-vtst-intel-2021.4.0` |
| 泛函、自旋 | PBE，ISPIN=2，无U；Ni初始1 μB，其余0；H2总自旋0 |
| 截断能 | ENCUT=520 eV |
| 展宽 | ISMEAR=0，SIGMA=0.05 |
| KPOINTS | 固体Gamma 3×6×1；H2 Gamma单点 |
| 电子收敛 | 固体EDIFF=1e-7；H2=1e-8，NELM=200 |
| 投影、精度 | LREAL=.FALSE.，PREC=Accurate，LASPH、ADDGRID打开 |
| 薄片修正 | LDIPOL=.TRUE.，IDIPOL=3，DIPOL=(0.5,0.5,0.5) |
| 端点/参考优化 | IBRION=2，POTIM=0.10，NSW=300，EDIFFG=-0.01 |
| 固体端点验收 | 最后SCF、正常结束、Relax离子收敛；最终可移动原子力<0.01 eV/Å |
| H2参考 | 20 Å盒子，4 MPI、16 GB、NBANDS=8；力<0.005 eV/Å |
| 固体NBANDS | Mg及其含H端点48；Mg2Ni相关及其含H端点96，路径全程相同 |

NBANDS的设置同时兼容端点64 MPI/NCORE=4与NEB每图像12 MPI/NCORE=4，避免不同MPI划分导致能带数自动取整发生变化。新H2参考在同一VTST二进制与H赝势下独立核验，不直接引用旧参考能量。

普通VASP模块的二进制未发现VTST参数标识，集群另外提供VTST模块。本次真正提交之前，驱动要求所加载程序包含LCLIMB、ICHAIN、IOPT标识并记录二进制SHA-256；每个计算节点启动时再次检查该二进制一致性。没有把普通VASP无条件称为CI-NEB。

## 5. 路径算法与资源

1. 必须先取得两端点的收敛几何。
2. 对相同H原子的编号交换、周期边界跨越作最小位移匹配；不允许换掉金属原子映射。
3. 使用H-H轴转动/伸长的几何初猜，再做有界IDPP预优化；它不是DFT轨迹。最终DFT图像放开全部非固定原子，不保留IDPP边界。
4. 每条路径生成00–06七个位置，其中00、06固定，01–05为五个中间图像。
5. 先非爬山NEB预优化到0.05 eV/Å，再CI-NEB到0.03 eV/Å。采用VTST FIRE：IBRION=3、POTIM=0、IOPT=7、SPRING=-5；250/400步上限。
6. 每条NEB用**一个节点60 MPI进程**，分为5图像×12进程，NCORE=4；不是错误地将64个进程平均分给5图像。每条路径24小时上限；参考态数组各64 MPI，12小时；端点数组各64 MPI，24小时。
7. 校验使用VTST有效力指标与CI最高图像受力，不把非最高图像的真实切向力误当成NEB未收敛。日志标记缺失或校验失败时保留文件并停止，不给出可靠势垒声明。

先前两种初始路径的离线几何测试已通过H守恒、周期匹配、固定原子一致、IDPP和避碰检查；纯Mg预览H-H约0.75→0.958→1.167→1.375→1.584→1.793→2.00 Å，Ni区域约0.753→1.005→1.255→1.514→1.773→2.031→2.269 Å。**这些只是未重新优化端点的预览，不是已计算的反应轨迹，也不是正式生产NEB图像。** 生产图像将在新端点收敛后重建。

## 6. 结果与结论口径

最终分别报告：

- 相对分子吸附态的正向电子能垒：Ea_ads = E_TS − E_H2_ads。
- 相对气相参考的过渡态高度：Ea_gas = E_TS − E_clean − E_H2。
- 解离反应能与逆向复合能垒，检查详细平衡关系。
- 两条路径采用相同参考定义的势垒差，而不是跨体系直接比较E_TS。
- 同时存储TOTEN和energy(sigma→0)，核对展宽影响。

没有内部能量峰、电子不收敛、有效力未达标、图像发生不合理穿越，均不得报告为“零势垒”或“已证明促进解离”。即使CI-NEB通过，本轮尚未做完整过渡态虚频和双向下降连通性检查，也没有零点能、温压修正和速率计算。不能直接把电子势垒等同于实验表观活化能。

## 7. 文件与追踪

生产目录结构：

```
H2_dissociation_compare_20261013/
  Mg/{clean,molecule,product,neb_coarse,ci_neb}/
  Mg2Ni_region/{clean,molecule,product,neb_coarse,ci_neb}/
  H2/clean/
  sources/
  submission.json
  manifest.json
```

尚未达到依赖条件的目录不会提前伪造已优化文件。原始已完成计算不覆盖；提交前后均保留锁、源结构哈希和逐项作业回执。断连导致回执不明确时不盲目重复提交；不设置后台自动重试。

只读检查：
`python3 ~/H2_dissociation_compare_20261013/neb_driver.py --status`

完整POTCAR只在集群许可目录和本批生产目录内使用；GitHub和下载包中只有非许可输入、校验值及日志，不含POTCAR或SSH凭据。

## 方法文档

- VASP NEB：https://vasp.at/wiki/Nudged_elastic_bands
- VTST CI-NEB：https://henkelmangroup.github.io/vtsttools/neb.html
- VTST FIRE/IOPT：https://henkelmangroup.github.io/vtsttools/optimizers.html

这些文档说明方法；本项目结构、作业状态和能量只能由本次实际计算记录证明。
