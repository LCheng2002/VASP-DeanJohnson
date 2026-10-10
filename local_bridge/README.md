# MacBook 免费提交思源一号：第一步（单个 Mg 静态计算）

本目录提供 **不使用 GitHub Actions** 的第一阶段试运行。由你自己的 MacBook 使用现有 SSH 证书连接思源一号，运行 `sbatch`。本脚本不会索取或复制你的 SSH 私钥、密码短语、GitHub Secrets，也不会将授权 POTCAR 上传 GitHub。

## 1. 在 MacBook 安装/检查 GitHub 登录

GitHub 与 ChatGPT 连接成功，不代表 **MacBook 上的 git** 已经登录 GitHub。推荐从本地终端执行：

```bash
git --version
gh --version
```

如 `gh` 未安装且你使用 Homebrew：`brew install gh`。然后：

```bash
gh auth login
gh auth setup-git
git clone https://github.com/LCheng2002/vasp-hpc-control.git "$HOME/vasp-hpc-control"
cd "$HOME/vasp-hpc-control"
```

如果已有克隆，改用 `cd "$HOME/vasp-hpc-control" && git pull --ff-only`。也可以通过 **GitHub Desktop** 克隆私人仓库，之后在终端进入该克隆目录。

## 2. 解锁现有 SSH 私钥

```bash
ssh-add --apple-use-keychain "$HOME/.ssh/sjtu_hpc_actions"
ssh-add -l
```

密码短语只在自己电脑终端输入。脚本使用 `BatchMode=yes`，不会停在远程交互式密码提示；如果自动认证不成功，会直接停止。

## 3. 只进行连接预检（不会消耗 Slurm 计算资源）

在仓库根目录执行：

```bash
bash local_bridge/submit_once.sh check
```

预期出现：用户 `deanfilan`、`sbatch` 路径、`sydata` SFTP 的工作目录、`CHECK_OK`。

如果提示 POTCAR 不可用，先检查思源一号 `~/VASP/potcar_profiles/Mg/POTCAR` 是否有合规授权的 Mg PBE 赝势，或确认 VASPKIT 103 可用于此账户。POTCAR 不上传 GitHub。

## 4. 手动提交这一个 Mg 静态任务

在确认资源与输入参数后执行：

```bash
bash local_bridge/submit_once.sh submit
```

只有返回 `SUBMITTED_JOB_ID|数字` 才代表已通过 `sbatch` 提交。脚本依次检查远端是否已经有 Job ID/锁、使用 `sydata` 的 SFTP 上传 4 个输入、在 `sylogin` 准备 POTCAR 并提交。

## 5. 查看作业

```bash
bash local_bridge/submit_once.sh status
```

实际作业运行在思源一号，MacBook 随后关机也不影响已提交的 Slurm 作业；GitHub Actions 不参与此次操作。

## 约束

- **只支持固定的 Mg 两原子静态任务** `Mg_bulk_static_20261010`。已固定校验 POSCAR、INCAR、KPOINTS、job.slurm 的 Git blob SHA-1，以防止不经确认就执行更改后的作业文件。
- Slurm 模块与赝势来源需要通过第一次真实任务核验。已有结构来自 CIF，不是 DFT relax 后结构；本次 K 点为试算设置。
- 远端目录：`~/VASP/github_jobs/Mg_bulk_static_20261010/`。
- `.submit_lock` 为防重复保护；如果作业结果不明确，**不要自行删锁并直接重试**。应先用 `squeue`、`sacct` 排查。
- 当前只实现 **手动提交 + 手动查询**。首次成功取得 Slurm Job ID 后，再加 macOS `launchd` 五分钟轮询与 GitHub 结果同步更稳妥。
- SSH 证书到期仍需按学校规定续签；此方案不绕过任何认证政策。
