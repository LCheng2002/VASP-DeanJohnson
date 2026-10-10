# VASP 原子结构库

该结构库最初由用户提供的 `归档.zip` 整理而来，包含 **8 个体相结构** 和 **5 个界面模型**；现另补充 **1 组 Mg/Mg17La2 界面**及 **2 组匹配的参考薄片**。
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

## 补充：Mg/Mg17La2 界面与参考薄片（`01_Mg_Mg17La2.zip`）

在此前的 8 个体相、5 组界面基础上，另加入以下三个 **原始参考/待计算模型**：

| 模型 | 组成 | 原子数 | 文件 | 结构类型 |
|---|---|---:|---|---|
| [Mg 参考薄片](surfaces/Mg_cleavage_reference/) | Mg54 | 54 | POSCAR, KPOINTS | 带约 18 Å 真空层的 cleavage slab |
| [Mg17La2 参考薄片](surfaces/Mg17La2_cleavage_reference/) | Mg68La8 | 76 | POSCAR, KPOINTS | 带约 18 Å 真空层的 cleavage slab，**非体相** |
| [Mg/Mg17La2 界面](interfaces/Mg_Mg17La2/) | Mg122La8 | 130 | POSCAR, KPOINTS, region_by_atom.json, metadata.json | 单个界面与约 15 Å 真空层 |

- [新增源文件清单、SHA-256 和结构自检信息](sources/01_Mg_Mg17La2/manifest.json)。
- `region_by_atom.json` 共含 130 个逐原子标签：54 个 Mg 区域，76 个 Mg17La2 区域；标签顺序与所提供的 POSCAR 原子列表对应。
- 三组模型源文件均采用 Cartesian 坐标和 Gamma 3×3×1 KPOINTS；K 点仅为归档原有设置，**未验证收敛性**，不是新的默认计算配置。
- 本压缩包没有 CONTCAR 或 OUTCAR，因此不声称完成了结构优化或计算收敛。
- 原归档包含的三个 POTCAR 文件因 VASP 许可限制未上传；没有任何计算任务被触发。
- **目前仍缺 Mg17La2 的真实三维体相结构**；这里保留的是 Mg17La2 切割薄片以及 Mg/Mg17La2 界面，不能用薄片直接代替体相进行计算。
