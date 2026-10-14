#!/usr/bin/env python3
"""Cancel ONLY the explicitly authorized 2026-10-14 CaFresh Slurm batch.

Guardrails:
- Verify exact production submission.json name-to-job mapping and full 23-item set.
- Only cancel currently active Slurm jobs belonging to current user and named CaFresh_*.
- Never cancel jobs in unrelated Ca, MgNi, LaH3, NEB or previous studies.
- Read-only post-check accounting; preserve all existing output and input directories.
"""
import datetime
import getpass
import json
import re
import subprocess
import time
from pathlib import Path

ROOT = Path.home() / "Ca_corrosion_fresh_20261014"
EXPECTED = {
    "Mg_bulk": "63681995",
    "H2_ref": "63681996",
    "H2O_ref": "63681997",
    "M0_clean": "63682000",
    "M0_water": "63682002",
    "M0_OH": "63682003",
    "M0_OH_H": "63682004",
    "M0_OH_water": "63682005",
    "M0_2OH_H": "63682006",
    "M0_2OH_near": "63682008",
    "M0_2OH_far": "63682009",
    "M0_4OH_near": "63682010",
    "M0_4OH_far": "63682011",
    "M1_clean": "63682012",
    "M1_water": "63682013",
    "M1_OH": "63682014",
    "M1_OH_H": "63682015",
    "M1_OH_water": "63682016",
    "M1_2OH_H": "63682021",
    "M1_2OH_near": "63682022",
    "M1_2OH_far": "63682023",
    "M1_4OH_near": "63682024",
    "M1_4OH_far": "63682026",
}
assert len(EXPECTED) == 23
NOW = lambda: datetime.datetime.now().astimezone().isoformat()


def call(args, timeout=40):
    p = subprocess.run(args, text=True, capture_output=True, timeout=timeout)
    if p.returncode:
        raise RuntimeError(f"COMMAND_FAILED|{args[0]}|rc={p.returncode}|{p.stderr.strip()[:500]}")
    return p.stdout


def active_queue():
    ids = ",".join(EXPECTED.values())
    raw = call(["squeue", "-h", "-j", ids, "-o", "%A|%u|%j|%T"])
    out = {}
    user = getpass.getuser()
    allowed_ids = set(EXPECTED.values())
    allowed_names = {"CaFresh_" + key for key in EXPECTED}
    for line in raw.splitlines():
        fields = [s.strip() for s in line.split("|")]
        if len(fields) != 4:
            raise RuntimeError("UNEXPECTED_SQUEUE_ROW|" + line)
        jid, owner, name, status = fields
        if jid not in allowed_ids or owner != user or name not in allowed_names:
            raise RuntimeError("UNEXPECTED_ACTIVE_JOB|" + line)
        out[jid] = dict(owner=owner, name=name, status=status)
    return out


def main():
    original = ROOT / "submission.json"
    if not original.is_file():
        raise RuntimeError("MISSING_PRODUCTION_LEDGER|refusing cancellation")
    ledger = json.loads(original.read_text())
    if ledger.get("complete") is not True or ledger.get("max_models") != 23:
        raise RuntimeError("INCOMPLETE_OR_WRONG_LEDGER|refusing cancellation")
    if ledger.get("jobs") != EXPECTED:
        raise RuntimeError("LEDGER_DOES_NOT_MATCH_ALL_23_EXPECTED_IDS|refusing cancellation")
    if not (ROOT / "manifest.json").is_file():
        raise RuntimeError("MISSING_PRODUCTION_MANIFEST|refusing cancellation")
    manifest = json.loads((ROOT / "manifest.json").read_text())
    if manifest.get("maximum_models") != 23 or manifest.get("preview"):
        raise RuntimeError("WRONG_PRODUCTION_MANIFEST|refusing cancellation")
    print("VERIFIED_BATCH|Ca_corrosion_fresh_20261014|count=23|" + NOW(), flush=True)
    active = active_queue()
    print("BEFORE_ACTIVE|" + json.dumps(active, ensure_ascii=False, sort_keys=True), flush=True)
    targets = sorted(active, key=int)
    if targets:
        print("APPLY_SCANCEL|count=" + str(len(targets)) + "|ids=" + ",".join(targets), flush=True)
        call(["scancel", *targets], timeout=45)
        print("SCANCEL_ACCEPTED|count=" + str(len(targets)), flush=True)
    else:
        print("NO_ACTIVE_JOBS|already finished or previously cancelled", flush=True)

    after = {}
    for attempt in range(6):
        after = active_queue()
        if not after:
            break
        time.sleep(2)
    print("AFTER_ACTIVE|" + json.dumps(after, ensure_ascii=False, sort_keys=True), flush=True)
    rows = call(["sacct", "-X", "-n", "-P", "-j", ",".join(EXPECTED.values()),
                 "-o", "JobID,JobName%35,State,ExitCode,Start,End"], timeout=45)
    print("ACCOUNTING_BEGIN")
    print(rows.strip(), flush=True)
    print("ACCOUNTING_END", flush=True)
    if after:
        print("NOT_YET_REMOVED_FROM_QUEUE|" + ",".join(sorted(after)), flush=True)
        raise RuntimeError("Some requested CaFresh jobs still active in Slurm after cancellation")
    print("CANCELLATION_VERIFIED|All_23_batch_ids_absent_from_squeue|" + NOW(), flush=True)
    print("OUTPUT_FILES_PRESERVED|No files modified or removed", flush=True)


if __name__ == "__main__":
    main()
