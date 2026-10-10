#!/usr/bin/env bash
# Free MacBook relay: GitHub clone -> sydata SFTP -> sylogin sbatch.
# Supports only one pinned Mg static job. Do not run on the HPC login node.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
TASK="Mg_bulk_static_20261010"
CALC="calculations/$TASK"
USER_NAME="${HPC_USER:-deanfilan}"
KEY="${HPC_SSH_IDENTITY:-$HOME/.ssh/sjtu_hpc_actions}"
CERT="${HPC_SSH_CERT_FILE:-$HOME/.ssh/sjtu_hpc_actions-cert.pub}"
MODE="${1:-check}"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "ERROR: This local script must run on MacBook, not the cluster." >&2
  exit 2
fi
[[ "$USER_NAME" =~ ^[a-zA-Z0-9._-]+$ ]] || { echo "Invalid HPC_USER" >&2; exit 2; }
for p in "$KEY" "$CERT" "$HOME/.ssh/known_hosts"; do
  test -s "$p" || { echo "Missing: $p" >&2; exit 2; }
done
for p in POSCAR INCAR KPOINTS job.slurm; do
  test -s "$CALC/$p" || { echo "Missing calculation file: $p" >&2; exit 2; }
done

# Pin the inputs for this one trial. Never execute silently changed job scripts.
check_sha() {
  p="$1"
  expected_sha="$2"
  actual_sha="$(git hash-object "$CALC/$p")"
  [[ "$actual_sha" == "$expected_sha" ]] || {
    echo "Input $p was changed since approval; refusing to submit." >&2
    exit 3
  }
}
check_sha POSCAR c33c8efd5116272ebfea023d21b4d9ceb7fe4fa3
check_sha INCAR e94cc193cbefd3abe2eb1533ab9649259f475b38
check_sha KPOINTS 058220a09049879e8a43e546470302356f4a0e10
check_sha job.slurm a1aca2cb014b28a27637900c1331f507b70de5fe

ssh_opts=(-F /dev/null
          -o BatchMode=yes
          -o IdentitiesOnly=yes
          -o StrictHostKeyChecking=yes
          -o ConnectTimeout=20
          -i "$KEY"
          -o "CertificateFile=$CERT")
login="$USER_NAME@sylogin.hpc.sjtu.edu.cn"
data="$USER_NAME@sydata.hpc.sjtu.edu.cn"

remote() { ssh "${ssh_opts[@]}" "$login" 'bash -s'; }
data_sftp() { sftp "${ssh_opts[@]}" -b - "$data"; }

case "$MODE" in
  check)
    echo "Checking non-interactive SSH and Slurm (no job submitted)..."
    remote <<'REMOTE_CHECK'
set -e
printf 'USER|%s\n' "$(whoami)"
printf 'HOST|%s\n' "$(hostname)"
command -v sbatch
if test -s "$HOME/VASP/potcar_profiles/Mg/POTCAR"; then
  echo 'POTCAR|Mg profile exists'
elif command -v vaspkit >/dev/null 2>&1; then
  echo 'POTCAR|VASPKIT found; generation has not been tested'
else
  echo 'POTCAR|Not found: configure licensed Mg POTCAR or VASPKIT first'
fi
REMOTE_CHECK
    printf 'pwd\n' | data_sftp
    echo "CHECK_OK|SSH/SFTP reached Siyuan, no sbatch executed"
    ;;

  submit)
    echo "Checking whether this exact task was already attempted..."
    remote <<'REMOTE_GUARD'
set -e
d="$HOME/VASP/github_jobs/Mg_bulk_static_20261010"
if test -s "$d/.slurm_jobid"; then
  printf 'ALREADY_SUBMITTED_JOB_ID|'; cat "$d/.slurm_jobid"
  exit 20
fi
if test -d "$d/.submit_lock" || test -e "$d/OUTCAR"; then
  echo 'EXISTING_JOB_OR_LOCK|Do not retry until verified manually'
  exit 21
fi
mkdir -p "$d"
REMOTE_GUARD

    echo "Uploading four input files via sydata SFTP..."
    data_sftp <<'SFTP_BATCH'
cd VASP/github_jobs/Mg_bulk_static_20261010
put calculations/Mg_bulk_static_20261010/POSCAR POSCAR
put calculations/Mg_bulk_static_20261010/INCAR INCAR
put calculations/Mg_bulk_static_20261010/KPOINTS KPOINTS
put calculations/Mg_bulk_static_20261010/job.slurm job.slurm
SFTP_BATCH

    echo "Preparing POTCAR on HPC and submitting one Slurm job..."
    remote <<'REMOTE_SUBMIT'
set -euo pipefail
cd "$HOME/VASP/github_jobs/Mg_bulk_static_20261010"
if test -s .slurm_jobid || test -d .submit_lock; then
  echo 'ALREADY_SUBMITTED_OR_LOCKED|Stopped without resubmission'
  exit 22
fi
for f in POSCAR INCAR KPOINTS job.slurm; do
  test -s "$f" || { echo "MISSING_INPUT|$f"; exit 23; }
done
if test -s "$HOME/VASP/potcar_profiles/Mg/POTCAR"; then
  cp "$HOME/VASP/potcar_profiles/Mg/POTCAR" POTCAR
elif command -v vaspkit >/dev/null 2>&1; then
  if ! timeout 120 vaspkit -task 103 </dev/null > vaspkit_potcar.log 2>&1; then
    echo 'VASPKIT_FAILED|See vaspkit_potcar.log'; exit 24
  fi
else
  echo 'POTCAR_UNAVAILABLE|Provide licensed Mg POTCAR or configure VASPKIT'; exit 25
fi
if ! test -s POTCAR || test "$(grep -c TITEL POTCAR)" -ne 1 ||
   ! grep -Eq 'TITEL[[:space:]]*=[[:space:]]*PAW_PBE[[:space:]]+Mg[[:space:]]' POTCAR; then
  echo 'POTCAR_MISMATCH|Expected one PAW_PBE Mg (not Mg_pv)'; exit 26
fi
# Persist the lock even if sbatch's result is ambiguous: never blindly resubmit.
mkdir .submit_lock
jid="$(sbatch --parsable job.slurm)"
printf '%s\n' "$jid" > .slurm_jobid.tmp
mv .slurm_jobid.tmp .slurm_jobid
printf 'SUBMITTED_JOB_ID|%s\n' "$jid"
REMOTE_SUBMIT
    ;;

  status)
    remote <<'REMOTE_STATUS'
set -euo pipefail
d="$HOME/VASP/github_jobs/Mg_bulk_static_20261010"
if ! test -s "$d/.slurm_jobid"; then
  echo 'NO_JOB_ID|No confirmed Slurm submission'; exit 0
fi
cd "$d"
jid="$(cat .slurm_jobid)"
printf 'JOB_ID|%s\n' "$jid"
jid="${jid%%;*}"
squeue -j "$jid" -o '%.18i %.20j %.12T %.10M %R' || true
sacct -X -j "$jid" --format=JobID,JobName%20,State,Elapsed,ExitCode || true
if test -s OUTCAR; then
  grep 'free  energy   TOTEN' OUTCAR | tail -n 1 || true
fi
REMOTE_STATUS
    ;;

  *)
    echo "Usage: bash local_bridge/submit_once.sh {check|submit|status}" >&2
    exit 2
    ;;
esac
