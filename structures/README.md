# VASP 原子结构库

该目录由用户提供的 `归档.zip` 归档整理而来，包含 **8 个体相结构** 和 **5 个界面模型**。
源压缩包 SHA-256：`4a0598c7734b8799da4ddbf4743569501320c1fcc5bf5072e3c7d3de43ac25f3`。

## 体相结构

| 结构 | 元素及单胞原子数 | 总原子数 | 文件 |
|---|---|---:|---|
| [Mg2Ni](bulk/Mg2Ni/) | Mg:12, Ni:6 | 18 | CIF + POSCAR |
| [MgH2](bulk/MgH2/) | Mg:2, H:4 | 6 | CIF + POSCAR |
| [CaH2](bulk/CaH2/) | Ca:4, H:8 | 12 | CIF + POSCAR |
| [Mg](bulk/Mg/) | Mg:2 | 2 | CIF + POSCAR |
| [LaH3](bulk/LaH3/) | La:4, H:12 | 16 | CIF + POSCAR |
| [CaMg2](bulk/CaMg2/) | Ca:4, Mg:8 | 12 | CIF + POSCAR |
| [LaH2](bulk/LaH2/) | La:4, H:8 | 12 | CIF + POSCAR |
| [Mg2NiH4](bulk/Mg2NiH4/) | Mg:16, Ni:8, H:32 | 56 | CIF + POSCAR |

`bulk/*/*.cif` 是源文件；对应的 `bulk/*/POSCAR` 根据同一个 CIF 的显式 P1 原子位点、晶胞参数和分数坐标转换得到，并非新的 DFT 优化结果。构建计算任务前，应检查结构合理性以及 POTCAR 元素顺序。

## 界面结构

| 模型 | 原子数 | 结构文件 |
|---|---:|---|
| [Mg_LaH3_Mg1120_LaH3_110](interfaces/Mg_LaH3_Mg1120_LaH3_110/) | 88 | POSCAR, CONTCAR |
| [Mg_Mg2Ni_Mg10-10_Mg2Ni_100](interfaces/Mg_Mg2Ni_Mg10-10_Mg2Ni_100/) | 42 | POSCAR, CONTCAR |
| [MgH2_Mg2Ni](interfaces/MgH2_Mg2Ni/) | 72 | POSCAR |
| [MgH2_Mg2NiH4](interfaces/MgH2_Mg2NiH4/) | 68 | POSCAR, CONTCAR |
| [MgH2_LaH3](interfaces/MgH2_LaH3/) | 100 | POSCAR, CONTCAR |

- `POSCAR`：源归档中提供的初始/待计算结构。
- `CONTCAR`：源归档中同目录存在的输出结构。**由于没有对应完整的 OUTCAR/vasprun.xml 收敛证据，此处不宣称已经收敛。**
- `MgH2_Mg2Ni` 仅提供 `POSCAR`，源归档中不存在对应的 `CONTCAR`。
- 部分 `CONTCAR` 尾部包含速度区块，按原文件保留。
- 原始界面文件的 `INCAR`、`KPOINTS` 和其他运行文件未作为基础结构入库，避免将历史计算参数误当成当前任务默认参数。

## 安全与完整性

- 不包含 VASP 许可限制的 `POTCAR`，也不包含源文件中的 TIFF 图像、`.DS_Store` 和 `__MACOSX` 资源叉文件。
- 不包含 `WAVECAR`、`CHGCAR` 或其他大型 VASP 计算输出。
- 使用 `manifest.json` 查看各文件的 SHA-256 校验值、单胞组成及源归档相对路径。
- 源归档中**没有 Mg17La2 的独立体相结构**，后续需从已有计算目录补充，不能通过当前文件推断或虚构。
- 体相文件名 `CaMg2` 沿用源 CIF 命名，在化学式上与 `Mg2Ca` 等价。

上传此结构库本身**不会触发思源一号的 Slurm 提交**。
