#!/usr/bin/env python3
"""Read-only compact status report for the exact MgNi/H4/B continuation jobs."""
import json
import os
import re
import subprocess
from pathlib import Path
from collections import deque

home=Path.home()
projects=[
 ("MgNi_cycle",home/"MgNi_abs_des_cycle_20261011",
  ["D_S1a_same_Ni","D_S1b_split_Ni","D_Sfar_Mg","A_H2_Ni_parallel","A_H2_Ni_vertical",
   "A_H2_Mg_parallel","A_H2_Mg_vertical","A_2H_Ni","A_2H_Mg","H2_reference",
   "D_S2_refill_one","D_S3_refill_two"]),
 ("B_retry",home/"Mg2Ni_Mg_H_thermo_validation_20261010/retry_20261011_B_v1",
  ["B-H-Mg_Mg_Mg2NiH4","B-H-Ni_Mg_Mg2NiH4"]),
 ("Mg17La2_H4",home/"TaskA_Mg17La2_H4_preference_20261011",
  ["H4_La_partner","H4_Mg_top","H4_Mg_bridge"])
]
ids="63565364,63565365,63565366,63554410,63554411,63554980"
print("READONLY_STATUS_AT|",subprocess.check_output(["date","+%Y-%m-%d_%H:%M:%S_%Z"],text=True).strip(),flush=True)
s=subprocess.run(["sacct","-X","-n","-P","-j",ids,"-o","JobIDRaw,JobName%25,State,ExitCode,Elapsed,Start,End"],capture_output=True,text=True,timeout=30)
print("SACCT_RETURN_CODE|",s.returncode,flush=True)
for line in s.stdout.splitlines():
 if line.strip() and (line.split("|")[0] in ["63565365","63565366","63554410","63554411"] or
                      re.fullmatch(r"63565364_\d+|63554980_\d+",line.split("|")[0])):
  print("JOB|"+line,flush=True)
sq=subprocess.run(["squeue","-u",os.environ.get("USER",""),"-h","-o","%i|%j|%T|%M|%R"],
 capture_output=True,text=True,timeout=20)
for line in sq.stdout.splitlines():
 if re.search(r"6356536[456]|6355441[01]|63554980",line):print("QUEUE|"+line,flush=True)
def read_json(p):
 try:return json.loads(p.read_text())
 except (OSError,ValueError):return None
def limited_tail(p,n=18):
 try:
  with p.open(errors="replace") as f:return [x.strip()[:250] for x in deque(f,n)]
 except OSError:return []
for label,root,names in projects:
 print("PROJECT|"+label+"|"+str(root.relative_to(home)),flush=True)
 for name in names:
  p=root/name
  if not p.is_dir():
   print("MODEL|"+label+"|"+name+"|MISSING",flush=True);continue
  ok=(p/"VERIFIED_FINISHED").exists()
  stages={}
  for phase in ["01_relax","02_static"]:
   q=p/phase; val=read_json(q/"validation.json")
   if val is not None:
    fields=["verified","numerically_verified","F_eV","E0_eV","ionic_steps","max_free_force_eV_A",
            "last_scf_ediff","last_scf_converged","ionic_converged","normal_end","min_H_H_A",
            "minimum_H_H_A","possible_H2_geometry"]
    info={k:val[k] for k in fields if k in val}
    # H4: show the added H80's nearest neighbors, not hundreds of atoms
    hn=val.get("final_H_neighbours")
    if isinstance(hn,list):
     matches=[x for x in hn if x.get("H_atom_1based")==80]
     if matches:info["H80_neighbours"]=matches[0].get("nearest6",[])[:5]
    stages[phase]=info
   elif q.is_dir():
    oszi=q/"OSZICAR"
    outs=q/"OUTCAR"
    stages[phase]={"validation":"MISSING",
      "OUTCAR_bytes":outs.stat().st_size if outs.exists() else 0,
      "OSZICAR_tail":limited_tail(oszi,3)}
  print("MODEL|"+label+"|"+name+"|"+("VERIFIED" if ok else "NOT_VERIFIED")+"|"+
        json.dumps(stages,ensure_ascii=False,separators=(",",":")),flush=True)
 if label=="MgNi_cycle":
  for name in ("D_S1b_split_Ni","A_H2_Mg_parallel","H2_reference"):
   p=root/name
   print("FAILURE_INSPECT|"+name,flush=True)
   for phase in ("01_relax","02_static"):
    d=p/phase
    if not d.is_dir():continue
    inc=d/"INCAR"
    if inc.is_file():
     tags=[]
     for line in inc.read_text(errors="replace").splitlines():
      if re.match(r"^\s*(NSW|POTIM|EDIFFG|IBRION|NELM)\s*=",line):tags.append(line.strip()[:120])
     print("INCAR|"+name+"|"+phase+"|"+";".join(tags),flush=True)
   for log in sorted(p.glob("**/slurm*.err"))[:4]+sorted(p.glob("**/slurm*.out"))[:4]:
    if log.is_file() and log.stat().st_size<1000000:
     print("SLURM_LOG|"+str(log.relative_to(root))+"|"+
           json.dumps(limited_tail(log,16),ensure_ascii=False),flush=True)
   for f in sorted(root.glob("slurm*"))[:12]:
    if f.is_file() and f.stat().st_size<1000000:
     print("ROOT_SLURM_LOG|"+f.name+"|"+json.dumps(limited_tail(f,8),ensure_ascii=False),flush=True)
  job=subprocess.run(["scontrol","show","job","63565366"],capture_output=True,text=True,timeout=15)
  print("REFILL_SCONTROL|"+(job.stdout[:600].replace("\n"," ") if job.returncode==0 else job.stderr[:160]),flush=True)
print("READONLY_DONE",flush=True)
