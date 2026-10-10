# Mg 单点静态能量（提交示例）

- 输入结构来自 [GitHub 的 hcp Mg 两原子 POSCAR](../../structures/bulk/Mg/POSCAR)，由先前 CIF 转换，尚非 DFT 优化结构。
- 用于一次初步静态能量计算，固定晶胞/原子（NSW=0）；此任务不进行 relax。
- INCAR 沿用此前研究常用设置 ENCUT=520、ISPIN=2、ISMEAR=1、SIGMA=0.2、EDIFF=1E-6；纯 Mg 不设置 U。
- Gamma 18×18×12 仅作为本次初步试算网格；如需可靠反应能或与过去的 -3.013855 eV Mg 参考数比较，必须核实之前的 POTCAR、KPOINTS、INCAR 和晶胞定义并开展收敛检查。
- 计算资源：思源一号 64c512g，单节点共享使用 8 个 MPI ranks，40 分钟上限，16GB 内存。
- 目标目录：`~/VASP/github_jobs/Mg_bulk_static_20261010`。
- **POTCAR 不在仓库内**：SSH 自动化仅使用个人账户已获授权的 `~/VASP/potcar_profiles/Mg/POTCAR`；若不存在但有可用的 VASPKIT，则在集群端尝试 `vaspkit -task 103`。要求 POTCAR 为单一 Mg 的 PBE 赝势，否则禁止提交。
- Slurm Job ID 出现在 GitHub Actions 日志的 `SUBMITTED_JOB_ID|...`，并保存在 HPC 的 `.slurm_jobid` 文件中。
- 自动重复运行时不会重复提交同一个任务（目录中的 `.submit_lock` / `.slurm_jobid` 保护）。
- GitHub Actions 未执行成功或没有返回 Job ID **不等于已提交**。登录 HPC 可用 `squeue -u "$USER"` 检查。
