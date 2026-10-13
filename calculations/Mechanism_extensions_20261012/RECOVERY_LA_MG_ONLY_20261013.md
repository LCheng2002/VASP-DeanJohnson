# MgH₂/LaH₃ 和纯 Mg 分支恢复记录（2026-10-13）

## 本次真实提交

原母体作业 `63585911` 中只有 3 个子任务通过了数值收敛验收：

| 母体 | Slurm | 最终结果 |
|---|---|---|
| `L_parent`，MgH₂/LaH₃ | `63585911_0` | COMPLETED / VERIFIED_FINISHED |
| `L_MgH2_reference`，MgH₂对照 | `63585911_1` | COMPLETED / VERIFIED_FINISHED |
| `Mg_parent`，纯Mg | `63585911_2` | COMPLETED / VERIFIED_FINISHED |
| `H2_reference`，气相H₂ | `63585912` | COMPLETED / VERIFIED_FINISHED |
| `X025_A` | `63585911_3` | TIMEOUT，未验收 |
| `X025_B` | `63585911_4` | TIMEOUT，未验收 |

原 `63585913`、`63585914`、`63585915` 均因整体母体数组未成功而被 Slurm 取消。原始 `63585915` **保持 CANCELLED，未重新运行**。

本次仅恢复有独立有效前置母体的 La 和纯 Mg 分支：

| 新 Slurm Job ID | 计算对象 | 子模型 | 初次提交后 sacct 状态 |
|---|---|---|---|
| **63636883** | MgH₂/LaH₃ | `L_vMg_near`、`L_vMg_far`、`L_vLa`、`L_MgH2_vMg`、`L_transfer`（5个，依次执行） | **PENDING** |
| **63636884** | 纯 Mg | `Mg_H2`、`Mg_2H`（2个，依次执行） | **PENDING** |

提交日志：https://github.com/LCheng2002/VASP-DeanJohnson/actions/runs/38410924942

前一次 GitHub Actions 执行曾因 `sacct JobIDRaw` 格式的解析未匹配数组成员而在提交前主动停止，**没有产生新 Job ID**。修正为 Slurm 的 `JobID` 格式并重新核验后才实际提交。

真实提交输出：
```text
SUBMITTED|La|63636883
SUBMITTED|Mg|63636884
SUBMISSION_COMPLETE|{"jobs":{"La":"63636883","Mg":"63636884"},"Hx_submitted":false}
```

提交记录保存在集群：
`~/MgNi_LaH3_extensions_20261012/recovery_La_Mg_only_20261013/submission.json`

执行期间没有修改原始 VASP 输入、没有提交或续算 `X025_A`、`X025_B`，没有恢复 Mg₂NiHₓ 的后续任务，没有重跑 `L_parent`、`L_MgH2_reference`、`Mg_parent` 或 `H2_reference`。后续结构只有在对应父模型已经 VERIFIED 的条件下由原计算驱动构建。原作业依赖已失效，恢复后的两个新 Job ID 不继承全数组依赖。

**状态边界：** 这是提交后的瞬时快照，不代表两组子模型已经收敛。真实结果以各目录的 `validation.json` 和 `VERIFIED_FINISHED` 标记为准。

## Mg₂NiHₓ 的后续（仅讨论，不提交）

用户明确要求暂停低氢 Mg₂NiHₓ 的 VASP 作业。数据库截图仅能确认存在多个候选化学式，不能仅凭化学式确定晶体结构是否为 Mg₂Ni 的同一氢化反应链。

- `Mg₂NiH₂` 与 `Mg₂NiH₄` 保持 Mg:Ni=2:1，是优先考虑的同组成氢化程度对照；但 `x=2` 不是低氢 `x≈0.3`，且是否可作为独立中间相必须根据稳定性和金属骨架验证。
- `MgNiH`、`MgNiH₂`、`MgNiH₃` 的 Mg:Ni=1:1，不应被冒充 Mg₂NiHₓ；可另行作为不同化学组成竞争相，需包含 Mg 化学势或相平衡反应。
- 现有 147 原子界面 Mg₂NiH₀.₂₅ 构型是人工插 H 的试探模型，两项均在 24 小时离子优化超时。因此未来若需重启，优先用更小 Mg₂Ni 金属母体超胞及文献四面体低氢占位（X3/X7）做 1–2 个结构测试，待能量稳定后再延伸界面。没有获得用户许可前不提交。
