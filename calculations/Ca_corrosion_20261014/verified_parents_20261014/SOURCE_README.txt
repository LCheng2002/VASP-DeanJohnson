Ca 腐蚀计算：四个 OH/H2O 模型 Relax/Static 归档
归档日期：2026-10-14

模型：M0_OH、M1_OH、M0_OH_H2O、M1_OH_H2O。
每个模型包含 calculation_inputs/（实际 Relax 输入 POSCAR、INCAR.relax、INCAR.static、KPOINTS、POTCAR、metadata.json）以及 results/（01_relax/OUTCAR、01_relax/CONTCAR、02_static/OUTCAR）。
M0_OH_H2O 和 M1_OH_H2O 的最终合格结果来自 continuation；归档输入 POSCAR 是续算起始结构。M0_OH 和 M1_OH 使用原始模型输入。
四个 Relax OUTCAR 均含 reached required accuracy；四个 Static OUTCAR 均正常结束并包含 energy(sigma->0)。
本地结果副本没有静态阶段的 CONTCAR、vasprun.xml、WAVECAR、CHGCAR 等文件，因此未将其列为已归档文件。
来源：oh_coverage_63131486_20261002/{results_remote,continuation}/。原始计算目录未修改。

最终 Static energy(sigma->0), eV:
M0_OH: -146.71669462 (original relax batch)
M1_OH: -147.63723064 (original relax batch)
M0_OH_H2O: -161.39391351 (continuation)
M1_OH_H2O: -162.53184536 (continuation)

SHA-256 清单见 SHA256SUMS.txt。
