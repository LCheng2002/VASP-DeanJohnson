# Siyuan SSH connection setup

The initial workflow is **manual and read-only**: [Test Siyuan HPC SSH](../.github/workflows/test-siyuan-ssh.yml). It does not submit VASP jobs, upload input files, or use `sbatch`.

## 1. Verify local certificate authentication

Use this on your Mac (replace the username only if your HPC account differs):

```bash
ssh -o IdentityAgent=none -o IdentitiesOnly=yes -o BatchMode=no \
  -i "$HOME/.ssh/sjtu_hpc_actions" \
  -o CertificateFile="$HOME/.ssh/sjtu_hpc_actions-cert.pub" \
  deanfilan@sylogin.hpc.sjtu.edu.cn
```

If the private key has a passphrase, SSH may prompt for it. Do not share any private key or passphrase in ChatGPT or commit them into a repository.

## 2. Populate GitHub Actions credentials

Open `Settings → Secrets and variables → Actions` in this **private repository**. Add these Repository **Secrets**:

| Secret | Exact value |
|---|---|
| `HPC_SSH_KEY` | Complete local file contents: `~/.ssh/sjtu_hpc_actions` (PRIVATE key; not `.pub`) |
| `HPC_SSH_CERT` | Complete local file contents: `~/.ssh/sjtu_hpc_actions-cert.pub` |
| `HPC_SSH_PASSPHRASE` | Passphrase used when generating the private key |
| `HPC_KNOWN_HOSTS` | Validated OpenSSH known_hosts entries for **both** `sylogin.hpc.sjtu.edu.cn` and `sydata.hpc.sjtu.edu.cn` |

Under `Variables`, add:

| Variable | Value |
|---|---|
| `HPC_USER` | `deanfilan` |

### Known-host verification

On your Mac you can collect *candidate* server keys:

```bash
ssh-keyscan -t ed25519 sylogin.hpc.sjtu.edu.cn sydata.hpc.sjtu.edu.cn > "$HOME/hpc_known_hosts_candidates"
ssh-keygen -lf "$HOME/hpc_known_hosts_candidates"
```

**Before** using these in `HPC_KNOWN_HOSTS`, independently verify their SHA256 fingerprints with the HPC administrator / platform documentation, and ensure they match the host keys you have previously trusted. A raw `ssh-keyscan` response alone is not authentication.

The workflow enforces `StrictHostKeyChecking=yes` and will fail rather than silently accept an unknown or changed host key.

## 3. Run the connection test

1. Open repository **Actions**.
2. Select **Test Siyuan HPC SSH**.
3. Choose **Run workflow** on `main`.
4. Open the `connection-test` job and inspect the output.

A successful run must show the return values of `whoami`, `hostname`, `command -v sbatch`, an SFTP `pwd` for `sydata`, and `SSH and SFTP connection tests passed.`

The initial test makes no computations and does not transfer VASP input files.

## 4. If it fails

- **Missing GitHub secret or variable**: Check the exact names in Step 2.
- **Error loading key / incorrect passphrase**: Confirm the secret is the **private** key and its original passphrase.
- **Permission denied after offering certificate**: The HPC server refused authentication; verify the certificate validity, principal, login policy and origin restrictions.
- **Connection timeout**: The GitHub-hosted runner is off-campus; ask the HPC administrator whether non-interactive cloud-runner SSH access is permitted. Do not bypass network or dynamic authorization policies.
- **Host key verification failed**: Re-check `HPC_KNOWN_HOSTS` against an authoritative fingerprint.
- **Certificate expired**: Obtain a renewed HPC-issued certificate.

GitHub-hosted Actions may have different access conditions from a Mac that has successfully logged in. A successful local SSH login is not proof the GitHub runner will be allowed.

## 5. Later: VASP submission

Only after this test passes should a separate workflow be enabled for copying **licensed POTCAR files from authorized HPC-local storage**, transferring inputs through `sydata`, submitting on `sylogin`, and recording the Slurm Job ID. Such a workflow needs input validation, explicit submission requests, duplicate-job prevention and careful privilege boundaries. No such submission has yet been configured.
