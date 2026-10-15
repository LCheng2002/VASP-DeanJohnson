#!/usr/bin/env python3
"""Checkpoint continuation for the two timed-out 2026-10-13 H2 NEB paths.

Preserves both original bands and all validated reference/endpoint calculations.
Starts two new coarse-NEB jobs from the five saved CONTCARs, with separate,
dependency-gated CI-NEB jobs only when the corresponding coarse band verifies.
POTCAR stays on the HPC filesystem. No automatic retries on a failed branch.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path.home() / "H2_dissociation_compare_20261013"
STAGE = ROOT / "resume_63614948_63614949_20261015"
sys.path.insert(0, str(ROOT))
import neb_core as c

JOBS = {"Mg": "63614948", "Mg2Ni_region": "63614949"}
MODULE = "vasp/6.3.2-vtst-intel-2021.4.0"
EXPECTED_BINARY_SHA = "454da6713051a916308a37843d23fb41cd31161285dcf09606a830a9f470a106"
SOURCE_DIR = {"Mg": ROOT / "Mg" / "neb_coarse",
              "Mg2Ni_region": ROOT / "Mg2Ni_region" / "neb_coarse"}

def need(condition, message):
    if not condition:
        raise RuntimeError("H2_NEB_CONTINUATION_STOP|" + message)

def now():
    return datetime.datetime.now().astimezone().isoformat()

def save(path, value):
    p = Path(path)
    temp = p.with_name(p.name + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temp.replace(p)

def call(argv, cwd=None, seconds=50):
    p = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, timeout=seconds)
    need(p.returncode == 0, "COMMAND_FAILED|" + " ".join(argv) + "|stderr=" + p.stderr[-450:])
    return p.stdout

def check_slurm():
    need((ROOT / "submission.json").is_file(), "Original ledger missing")
    ledger = json.loads((ROOT / "submission.json").read_text())
    need(ledger.get("all_submitted") is True, "Original submission incomplete")
    need(ledger.get("jobs", {}).get("Mg") == JOBS["Mg"] and
         ledger.get("jobs", {}).get("Mg2Ni") == JOBS["Mg2Ni_region"],
         "Original Slurm IDs do not match production manifest")
    ids = ",".join(JOBS.values())
    stdout = call(["sacct", "-X", "-n", "-P", "-j", ids,
                   "-o", "JobID,JobName%32,State,ExitCode,Elapsed,Timelimit,Start,End"])
    rows = {parts[0]: parts for line in stdout.splitlines()
            if (parts := line.split("|")) and parts[0] in JOBS.values()}
    need(len(rows) == 2, "Both original jobs must be independently present in sacct")
    for sysname, jid in JOBS.items():
        need(rows[jid][2] == "TIMEOUT", sysname + " is no longer confirmed TIMEOUT: " + rows[jid][2])
        need(rows[jid][1] == ("h2neb_Mg" if sysname == "Mg" else "h2neb_Mg2Ni"),
             "Original job name mismatch for " + jid)
    pending = call(["squeue", "-h", "-j", ids, "-o", "%i|%j|%T"])
    need(not pending.strip(), "Original job still in the Slurm queue")
    other = call(["squeue", "-h", "-u", os.environ.get("USER", "deanfilan"), "-o", "%i|%j|%T"])
    need(not any("h2neb_resume_" in line for line in other.splitlines()),
         "Another H2 NEB continuation is already queued")
    return rows

def check_inputs():
    need(ROOT.is_dir(), "Original production directory unavailable")
    manifest = json.loads((ROOT / "manifest.json").read_text())
    need(manifest.get("preview") is False, "Cannot restart a preview calculation")
    backend = manifest.get("backend", {})
    need(backend.get("sha256") == EXPECTED_BINARY_SHA and backend.get("module") == MODULE,
         "Original production binary differs")
    for name in ("neb_core.py", "neb_driver.py", "plan.json"):
        need(c.sha(ROOT / name) == manifest["input_hashes"][name],
             "Original pinned control file modified: " + name)
    for sysname in JOBS:
        base = ROOT / sysname
        for parent in ("clean", "molecule", "product"):
            need((base / parent / "VERIFIED_FINISHED").is_file(),
                 "Original verified " + sysname + "/" + parent + " missing")
            v = json.loads((base / parent / "static" / "validation.json").read_text())
            need(v.get("verified") is True, "Original parent numerical gate failed")
        need((base / "ENDPOINTS_VERIFIED").is_file(), "Original endpoints gate missing")
        need((base / ".neb_started").is_dir(), "Original NEB run was not started")
        need(not (base / "NEB_VERIFIED").exists(), "Original band is already verified")
        src = SOURCE_DIR[sysname]
        for file in ("INCAR", "KPOINTS", "POTCAR"):
            p = src / file
            need(p.is_file() and p.stat().st_size > 10, "Original band input missing: " + str(p))
        ends = (c.parse((src / "00" / "POSCAR").read_text()),
                c.parse((src / "06" / "POSCAR").read_text()))
        c.endpoint_check(ends[0], "molecule")
        c.endpoint_check(ends[1], "product")
        for i in range(1, 6):
            folder = src / f"{i:02d}"
            for fname in ("POSCAR", "CONTCAR", "OUTCAR"):
                need((folder / fname).is_file() and (folder / fname).stat().st_size > 100,
                     "Original checkpoint missing: " + str(folder / fname))
            start = c.parse((folder / "POSCAR").read_text())
            end = c.parse((folder / "CONTCAR").read_text())
            need(start["sp"] == end["sp"] == ends[0]["sp"], "Image chemistry changed")
            need(np.array_equal(start["mask"], end["mask"]) and
                 np.array_equal(start["mask"], ends[0]["mask"]),
                 "Image selective-dynamics constraint changed")
            need(np.allclose(start["cell"], end["cell"], atol=1e-7) and
                 np.allclose(ends[0]["cell"], end["cell"], atol=1e-7),
                 "Image cell changed")
            delta = c.mic(end["xyz"] - start["xyz"], end["cell"])
            need(np.max(np.abs(delta[~start["mask"]]), initial=0.0) < 1e-5,
                 "Fixed substrate atoms moved")
            geometry = c.summary_geometry(end)
            need(geometry["H_H_A"] > 0.55, "Hydrogen atoms collide")
            need(min(h["nearest_metals"][0][0] for h in geometry["H_neighbours"]) > 0.65,
                 "H-metal atom collision")
    need((ROOT / "H2" / "clean" / "VERIFIED_FINISHED").is_file(), "Original H2 reference missing")
    return True

def precheck():
    rows = check_slurm()
    check_inputs()
    print("CONTINUATION_PREFLIGHT_OK|" +
          json.dumps({"checked_at": now(), "original_jobs": {k: rows[j][2] for k,j in JOBS.items()},
                      "each_band": "5 validated checkpoint geometries + 2 fixed endpoints",
                      "no_completed_parent_recomputed": True}, ensure_ascii=False), flush=True)

def slurm_script(sysname, phase):
    mode = "coarse" if phase == "coarse" else "ci"
    name = "mg" if sysname == "Mg" else "ni"
    return f"""#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=60
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi {MODULE}
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
python3 resume_H2_NEB_63614948_63614949_20261015.py --run {sysname} {mode}
"""

def copy_band_from_checkpoint(sysname, target):
    source = SOURCE_DIR[sysname]
    need(not target.exists(), "Refuse to replace existing continuation band")
    target.mkdir(parents=True)
    for fname in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(source / fname, target / fname)
    for i in range(7):
        folder = target / f"{i:02d}"
        folder.mkdir()
        source_path = source / f"{i:02d}" / ("CONTCAR" if 1 <= i <= 5 else "POSCAR")
        shutil.copy2(source_path, folder / "POSCAR")
    return {f"{i:02d}/POSCAR": c.sha(target / f"{i:02d}" / "POSCAR") for i in range(7)}

def stage():
    need(not STAGE.exists(), "Continuation directory exists: inspect ledger; do not overwrite or resubmit")
    precheck()
    STAGE.mkdir(mode=0o700)
    (STAGE / "logs").mkdir()
    shutil.copy2(Path(__file__), STAGE / "resume_H2_NEB_63614948_63614949_20261015.py")
    records = {}
    for sysname in JOBS:
        src = SOURCE_DIR[sysname]
        target = STAGE / sysname / "coarse"
        hashes = copy_band_from_checkpoint(sysname, target)
        hashes.update({f: c.sha(target / f) for f in ("INCAR", "KPOINTS", "POTCAR")})
        records[sysname] = {"source_job": JOBS[sysname],
                            "source_checkpoint_sha256": {f"{i:02d}": c.sha(src / f"{i:02d}" / "CONTCAR")
                                                          for i in range(1,6)},
                            "input_hashes": hashes}
        for phase in ("coarse", "ci"):
            fn = ("resume_mg" if sysname == "Mg" else "resume_ni") + "_" + phase + ".slurm"
            (STAGE / fn).write_text(slurm_script(sysname, phase))
            call(["bash", "-n", fn], cwd=STAGE)
    save(STAGE / "manifest.json",
         {"created_at": now(), "original_ids": JOBS,
          "binary_sha256": EXPECTED_BINARY_SHA, "scope": ["Mg", "Mg2Ni_region"],
          "original_preserved": True, "restarts": records})
    print("CONTINUATION_STAGED|" + str(STAGE), flush=True)

def submit():
    if STAGE.exists():
        ledger = STAGE / "submission.json"
        if ledger.exists():
            data = json.loads(ledger.read_text())
            need(data.get("all_submitted") is True,
                 "Partial prior submission: manual review required; no duplicate jobs")
            print("EXISTING_SUBMISSION|" + json.dumps(data), flush=True)
            return
        raise RuntimeError("Continuation directory exists without ledger: manual review required")
    stage()
    manifest = json.loads((STAGE / "manifest.json").read_text())
    for sysname, record in manifest["restarts"].items():
        folder = STAGE / sysname / "coarse"
        for rel, original_sha in record["input_hashes"].items():
            need(c.sha(folder / rel) == original_sha,
                 "Prepared continuation input unexpectedly modified")
    (STAGE / ".submission_attempted").mkdir()
    ledger = {"original": JOBS, "created_at": now(),
              "jobs": {}, "all_submitted": False, "resume_mode": "coarse_then_afterok_CI"}
    save(STAGE / "submission.json", ledger)
    for sysname in JOBS:
        prefix = "mg" if sysname == "Mg" else "ni"
        for phase in ("coarse", "ci"):
            key = sysname + "_" + phase
            args = ["sbatch", "--parsable", "--job-name=h2neb_resume_" + prefix + "_" + phase]
            if phase == "ci":
                dependency = ledger["jobs"][sysname + "_coarse"]
                args.extend(["--dependency=afterok:" + dependency, "--kill-on-invalid-dep=yes"])
            args.append("resume_" + prefix + "_" + phase + ".slurm")
            ledger["attempting"] = key
            save(STAGE / "submission.json", ledger)
            output = call(args, cwd=STAGE).strip()
            ident = output.split(";")[0]
            need(re.fullmatch(r"[0-9]+", ident) is not None, "Ambiguous sbatch response " + output)
            ledger["jobs"][key] = ident
            ledger.pop("attempting", None)
            save(STAGE / "submission.json", ledger)
            print("SUBMITTED|" + key + "|" + ident +
                  (("|afterok=" + dependency) if phase == "ci" else "|independent"), flush=True)
    ledger["all_submitted"] = True
    ledger["submitted_at"] = now()
    save(STAGE / "submission.json", ledger)
    print("H2_NEB_CONTINUATION_SUBMITTED|" + json.dumps(ledger), flush=True)

def binary_check():
    need(int(os.environ.get("SLURM_NTASKS", "0")) == 60,
         "VASP NEB requires precisely 60 allocated MPI processes")
    binary = shutil.which("vasp_std")
    need(binary is not None, "VTST VASP unavailable")
    need(c.sha(Path(binary).resolve()) == EXPECTED_BINARY_SHA,
         "VTST VASP executable differs from original production")

def verify_band(sysname, phase, folder):
    limit = 0.05 if phase == "coarse" else 0.03
    report = []
    for i in range(1,6):
        q = folder / f"{i:02d}"
        inp = c.parse((q / "POSCAR").read_text())
        out = c.parse((q / "CONTCAR").read_text())
        need(inp["sp"] == out["sp"] and np.array_equal(inp["mask"], out["mask"]) and
             np.allclose(inp["cell"], out["cell"], atol=1e-7),
             "Final image chemistry/cell/mask mismatch: " + str(q))
        need(np.max(np.abs(c.mic(out["xyz"]-inp["xyz"], out["cell"])[~inp["mask"]]),
                    initial=0.0) < 1e-5, "Fixed substrate moved: " + str(q))
        result = c.output(q, inp["mask"])
        good = bool(result.get("normal_end") and result.get("last_ediff") and
                    result.get("neb_marker") and result.get("lclimb") == (phase == "ci") and
                    result.get("vtst_force_max") is not None and
                    result["vtst_force_max"] < limit and
                    result.get("F_eV") is not None and math.isfinite(result["F_eV"]))
        result["verified"] = good
        result["H_H_A"] = c.summary_geometry(out)["H_H_A"]
        c.save(q / "validation.json", result)
        report.append(result)
    verified = all(item["verified"] for item in report)
    c.save(folder / "validation.json",
           {"system": sysname, "phase": phase,
            "force_tolerance_eV_A": limit, "verified": verified, "images": report})
    need(verified, "NEB path not converged and verified: " + str(folder))
    (folder / "BAND_VERIFIED").touch()
    print("BAND_VERIFIED|" + sysname + "|" + phase, flush=True)
    return report

def prepare_ci(sysname):
    parent = STAGE / sysname / "coarse"
    need((parent / "BAND_VERIFIED").exists(), "Cannot begin CI before a verified coarse NEB")
    source = json.loads((parent / "validation.json").read_text())
    need(source.get("verified") is True, "Coarse band validation not passed")
    destination = STAGE / sysname / "ci"
    need(not destination.exists(), "Existing CI band would be overwritten")
    destination.mkdir(parents=True)
    for fname in ("KPOINTS", "POTCAR"):
        shutil.copy2(parent / fname, destination / fname)
    start = c.parse((parent / "00" / "POSCAR").read_text())
    nband = 48 if sysname == "Mg" else 96
    (destination / "INCAR").write_text(c.input_text(c.incar(sysname, "neb", nband, False, "ci"), start["sp"]))
    for i in range(7):
        q = destination / f"{i:02d}"
        q.mkdir()
        src = parent / f"{i:02d}" / ("CONTCAR" if 1 <= i <= 5 else "POSCAR")
        shutil.copy2(src, q / "POSCAR")
    print("CI_STAGED|" + sysname, flush=True)

def barrier(sysname, values):
    reference = [json.loads((ROOT/sysname/name/"static"/"validation.json").read_text())
                 for name in ("molecule", "product", "clean")]
    gas = json.loads((ROOT/"H2"/"clean"/"static"/"validation.json").read_text())
    need(all(x.get("verified") for x in reference+[gas]), "Reference/endpoint verification missing")
    F = [reference[0]["F_eV"]] + [v["F_eV"] for v in values] + [reference[1]["F_eV"]]
    E0 = [reference[0]["E0_eV"]] + [v["E0_eV"] for v in values] + [reference[1]["E0_eV"]]
    top = int(np.argmax(F))
    need(0 < top < 6, "No interior activation maximum; do not declare a zero barrier")
    need(values[top-1]["max_free_DFT_force"] < 0.04,
         "Highest image is not stationary in raw DFT force")
    result = {
        "status":"CI_NEB_NUMERICALLY_VERIFIED_TS_VIBRATIONS_NOT_DONE",
        "system":sysname, "highest_image":top, "image_F_eV":F,
        "image_E0_eV":E0, "Ea_from_adsorbed_F_eV":F[top]-F[0],
        "Ea_from_adsorbed_E0_eV":E0[top]-E0[0],
        "Ea_reverse_F_eV":F[top]-F[-1], "reaction_F_eV":F[-1]-F[0],
        "Ea_from_gas_F_eV":F[top]-reference[2]["F_eV"]-gas["F_eV"],
        "Ea_from_gas_E0_eV":E0[top]-reference[2]["E0_eV"]-gas["E0_eV"],
        "limitations":"No TS vibration, ZPE, temperature/pressure or rate constant validation."
    }
    save(STAGE / sysname / "barrier_summary.json", result)
    (STAGE / sysname / "NEB_VERIFIED").touch()
    print("BARRIER_VERIFIED|" + json.dumps(result), flush=True)

def run_band(sysname, phase):
    binary_check()
    ledger = json.loads((STAGE / "submission.json").read_text())
    need(ledger.get("jobs", {}).get(sysname + "_" + phase),
         "No recorded Slurm submission for branch")
    manifest = json.loads((STAGE / "manifest.json").read_text())
    if phase == "ci":
        prepare_ci(sysname)
    folder = STAGE / sysname / phase
    if phase == "coarse":
        for rel, sha in manifest["restarts"][sysname]["input_hashes"].items():
            need(c.sha(folder / rel) == sha, "Coarse restart input changed: " + rel)
    with (folder/"vasp.stdout").open("w") as out, (folder/"vasp.stderr").open("w") as err:
        p = subprocess.run(["mpirun", "-np", "60", "vasp_std"],
                           cwd=folder, stdout=out, stderr=err)
    need(p.returncode == 0, "VASP returned nonzero: " + str(p.returncode))
    values = verify_band(sysname, phase, folder)
    if phase == "ci":
        barrier(sysname, values)

def status():
    print("CHECKED_AT|" + now(), flush=True)
    p = STAGE / "submission.json"
    if not p.is_file():
        print("NOT_SUBMITTED")
        return
    j = json.loads(p.read_text())
    print("SUBMISSION|" + json.dumps(j, ensure_ascii=False))
    ids = ",".join(j.get("jobs", {}).values())
    if ids:
        print("SQUEUE|" + call(["squeue","-r","-h","-j",ids,"-o","%i|%j|%T|%M|%R"]).replace("\n"," || "),flush=True)
        print("SACCT|" + call(["sacct","-X","-n","-P","-j",ids,
                              "-o","JobID,JobName%28,State,ExitCode,Elapsed,Start,End"]).replace("\n"," || "),flush=True)
    for sysname in JOBS:
        for phase in ("coarse","ci"):
            p = STAGE / sysname / phase / "validation.json"
            if p.exists():
                data=json.loads(p.read_text())
                print("VALIDATION|" + sysname + "|" + phase +
                      "|" + json.dumps({"verified":data.get("verified"),
                                        "max_effective_force":max((x.get("vtst_force_max") or 0) for x in data.get("images",[]))}),flush=True)

def main():
    a = argparse.ArgumentParser()
    group = a.add_mutually_exclusive_group(required=True)
    group.add_argument("--precheck", action="store_true")
    group.add_argument("--submit", action="store_true")
    group.add_argument("--status", action="store_true")
    group.add_argument("--run", nargs=2, metavar=("SYSTEM", "PHASE"))
    x = a.parse_args()
    if x.precheck:
        precheck()
    elif x.submit:
        submit()
    elif x.status:
        status()
    else:
        need(x.run[0] in JOBS and x.run[1] in ("coarse", "ci"), "Unknown continuation branch")
        run_band(*x.run)

if __name__ == "__main__":
    main()
