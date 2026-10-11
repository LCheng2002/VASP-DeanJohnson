#!/usr/bin/env python3
"""Submit only the already-prepared fixed repair roster; never changes old jobs."""
import datetime, hashlib, json, os, re, shutil, subprocess
from pathlib import Path
O=Path.home()/'MgNi_abs_des_cycle_20261011';R=O/'repairs_20261012_v1'
assert (R/'input_manifest.json').is_file(),'Run the input preparer first'
ledger=R/'submission.json'
if ledger.is_file():
 old=json.loads(ledger.read_text())
 if old.get('all_submitted'):
  print('ALREADY_SUBMITTED|'+json.dumps(old));raise SystemExit(0)
 raise RuntimeError('Partial or uncertain submission: inspect ledger, do not blindly retry')
assert not (R/'.submission_lock').exists(),'Ambiguous earlier submission'
for rel,h in json.loads((R/'input_manifest.json').read_text()).items():
 assert hashlib.sha256((R/rel).read_bytes()).hexdigest()==h,'Prepared input changed: '+rel
q=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%j'],capture_output=True,text=True,check=True,timeout=30)
assert not any(n.startswith('mgni_fix_') for n in q.stdout.splitlines()),'A repair is already active'
runner=Path(__file__).resolve().parent/'run_cycle_repairs_20261012.py'
assert runner.is_file();compile(runner.read_text(),str(runner),'exec')
shutil.copy2(runner,R/runner.name)
roster=[('S1b_polish','0',64,'12:00:00','mgni_fix_S1b',None),('MgH2_parallel_polish','1',64,'12:00:00','mgni_fix_MgH2',None),('H2_memory_repair','2',4,'00:30:00','mgni_fix_H2ref','16G'),('refill_recovery','refill',64,'24:00:00','mgni_fix_refill',None)]
for key,mode,ranks,wall,name,mem in roster:
 s='#!/bin/bash\n#SBATCH --partition=64c512g\n#SBATCH --nodes=1\n#SBATCH --ntasks='+str(ranks)+'\n#SBATCH --job-name='+name+'\n#SBATCH --time='+wall+'\n'
 s+=('#SBATCH --mem='+mem+'\n') if mem else '#SBATCH --exclusive\n'
 s+='#SBATCH --output=logs/%x-%j.out\n#SBATCH --error=logs/%x-%j.err\nset -euo pipefail\ncd "$SLURM_SUBMIT_DIR"\nmodule purge\nmodule load oneapi vasp/6.3.2-intel-2021.4.0\nulimit -s unlimited\nulimit -l unlimited\nexport OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\npython3 run_cycle_repairs_20261012.py '+mode+'\n'
 f=R/(key+'.slurm');f.write_text(s);subprocess.run(['bash','-n',str(f)],check=True)
(R/'.submission_lock').mkdir()
record={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'directory':str(R),'jobs':{},'runner_sha256':hashlib.sha256(runner.read_bytes()).hexdigest()}
def save():
 t=ledger.with_suffix('.tmp');t.write_text(json.dumps(record,indent=2));t.replace(ledger)
save()
for key,mode,ranks,wall,name,mem in roster:
 (R/('.attempted_'+key)).touch(exist_ok=False)
 cmd=['sbatch','--parsable']
 if key=='refill_recovery':
  dep=record['jobs']['S1b_polish'];cmd+=['--dependency=afterok:'+dep,'--kill-on-invalid-dep=yes'];record['refill_dependency']='afterok:'+dep
 cmd.append(key+'.slurm')
 p=subprocess.run(cmd,cwd=R,capture_output=True,text=True,timeout=45)
 if p.returncode:raise RuntimeError('sbatch failed '+key+': '+p.stderr[:400])
 jid=p.stdout.strip();assert re.fullmatch(r'[0-9]+',jid),'Unexpected job acknowledgement'
 record['jobs'][key]=jid;save();print('SUBMITTED|'+key+'|'+jid,flush=True)
record['all_submitted']=True;save();print('REPAIR_SUBMISSION_JSON|'+json.dumps(record),flush=True)
