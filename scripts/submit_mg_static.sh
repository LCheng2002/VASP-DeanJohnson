#!/usr/bin/env bash
set -euo pipefail

TASK="Mg_bulk_static_20261010"
LOCAL_DIR="calculations/$TASK"
REMOTE_DIR="VASP/github_jobs/$TASK"

for v in HPC_USER HPC_SSH_KEY HPC_SSH_CERT HPC_SSH_PASSPHRASE HPC_KNOWN_HOSTS; do
  if [[ -z "${!v:-}" ]]; then
    echo "::error::Missing secret or variable: $v"
    exit 1
  fi
done

for f in POSCAR INCAR KPOINTS job.slurm; do
  test -s "$LOCAL_DIR/$f" || { echo "::error::Missing input: $f"; exit 1; }
done
awk 'NR==6 {if ($1 != "Mg" || NF != 1) exit 1}
     NR==7 {if ($1 != 2 || NF != 1) exit 1}' "$LOCAL_DIR/POSCAR"

umask 077
mkdir -p "$HOME/.ssh"
chmod 700 "$HOME/.ssh"
printf '%s\n' "$HPC_SSH_KEY" > "$HOME/.ssh/hpc_automation"
printf '%s\n' "$HPC_SSH_CERT" > "$HOME/.ssh/hpc_automation-cert.pub"
printf '%s\n' "$HPC_KNOWN_HOSTS" > "$HOME/.ssh/known_hosts"
chmod 600 "$HOME/.ssh/hpc_automation" "$HOME/.ssh/known_hosts"
chmod 644 "$HOME/.ssh/hpc_automation-cert.pub"

cat > "$RUNNER_TEMP/hpc-askpass.sh" <<'ASKPASS'
#!/bin/sh
printf '%s' "$HPC_SSH_PASSPHRASE"
ASKPASS
chmod 700 "$RUNNER_TEMP/hpc-askpass.sh"

eval "$(ssh-agent -s)"
cleanup() {
  ssh-agent -k >/dev/null 2>&1 || true
  rm -f "$HOME/.ssh/hpc_automation" "$HOME/.ssh/hpc_automation-cert.pub" "$RUNNER_TEMP/hpc-askpass.sh"
}
trap cleanup EXIT

DISPLAY=:0 SSH_ASKPASS="$RUNNER_TEMP/hpc-askpass.sh" SSH_ASKPASS_REQUIRE=force \
  ssh-add "$HOME/.ssh/hpc_automation" </dev/null

ssh_opts=(
  -F /dev/null
  -o BatchMode=yes
  -o IdentitiesOnly=yes
  -o StrictHostKeyChecking=yes
  -o ConnectTimeout=15
  -i "$HOME/.ssh/hpc_automation"
  -o CertificateFile="$HOME/.ssh/hpc_automation-cert.pub"
)

echo "STAGE|Preparing fresh HPC calculation directory"
ssh "${ssh_opts[@]}" "$HPC_USER@sylogin.hpc.sjtu.edu.cn" \
  'set -e; d="$HOME/VASP/github_jobs/Mg_bulk_static_20261010"; if test -e "$d/.slurm_jobid" || test -d "$d/.submit_lock"; then echo "ALREADY_ATTEMPTED: $d"; exit 20; fi; mkdir -p "$d"'

echo "STAGE|Uploading four text inputs via sydata"
rsync -av \
  -e "ssh -F /dev/null -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=15 -i $HOME/.ssh/hpc_automation -o CertificateFile=$HOME/.ssh/hpc_automation-cert.pub" \
  -- "$LOCAL_DIR/POSCAR" "$LOCAL_DIR/INCAR" "$LOCAL_DIR/KPOINTS" "$LOCAL_DIR/job.slurm" \
  "$HPC_USER@sydata.hpc.sjtu.edu.cn:$REMOTE_DIR/"

echo "STAGE|Generating/locating licensed POTCAR and submitting Slurm"
ssh "${ssh_opts[@]}" "$HPC_USER@sylogin.hpc.sjtu.edu.cn" 'bash -s' <<'REMOTE'
set -euo pipefail
cd "$HOME/VASP/github_jobs/Mg_bulk_static_20261010"
if test -s .slurm_jobid; then
  echo "ALREADY_SUBMITTED_JOB_ID|$(cat .slurm_jobid)"
  exit 0
fi
if ! mkdir .submit_lock 2>/dev/null; then
  echo "SUBMISSION_LOCK_EXISTS|No resubmission attempted."
  exit 25
fi
submitted=0
cleanup_lock() {
  if [[ "$submitted" -eq 0 ]]; then
    rmdir .submit_lock 2>/dev/null || true
  fi
}
trap cleanup_lock EXIT
for f in POSCAR INCAR KPOINTS job.slurm; do
  test -s "$f" || { echo "MISSING_INPUT|$f"; exit 26; }
done

if test -s "$HOME/VASP/potcar_profiles/Mg/POTCAR"; then
  cp "$HOME/VASP/potcar_profiles/Mg/POTCAR" POTCAR
  echo "POTCAR_SOURCE|authorized user-local profile"
elif command -v vaspkit >/dev/null 2>&1; then
  echo "POTCAR_SOURCE|vaspkit -task 103"
  if ! timeout 120 vaspkit -task 103 </dev/null > vaspkit_potcar.log 2>&1; then
    echo "POTCAR_GENERATION_FAILED|Check ~/.vaspkit, PBE_PATH and vaspkit_potcar.log on HPC"
    exit 27
  fi
else
  echo "POTCAR_UNAVAILABLE|Provide licensed ~/VASP/potcar_profiles/Mg/POTCAR or configure VASPKIT on HPC"
  exit 28
fi

test -s POTCAR || { echo "POTCAR_MISSING"; exit 29; }
test "$(grep -c 'TITEL' POTCAR)" -eq 1 || { echo "POTCAR_ELEMENT_COUNT_INVALID"; exit 30; }
grep -Eq 'TITEL[[:space:]]*=[[:space:]]*PAW_PBE[[:space:]]+Mg[[:space:]]' POTCAR || {
  echo "POTCAR_TYPE_MISMATCH|Expected single PAW_PBE Mg potential (not Mg_pv)."
  exit 31
}
grep 'TITEL' POTCAR | sed 's/^/POTCAR_TITLE|/'

jid=$(sbatch --parsable job.slurm)
submitted=1
printf '%s\n' "$jid" > .slurm_jobid.tmp
mv .slurm_jobid.tmp .slurm_jobid
echo "SUBMITTED_JOB_ID|$jid"
REMOTE
