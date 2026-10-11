#!/usr/bin/env python3
"""Prepare the three diagnosed repairs in a new, fixed HPC directory; no sbatch."""
import hashlib, importlib.util, json, shutil, subprocess
from pathlib import Path
O=Path.home()/'MgNi_abs_des_cycle_20261011'
R=O/'repairs_20261012_v1'
PIN='f5696bda76eb2c990e5daab3a09baabf8619a914bf5bed10b37a2a42aa495abc'
SPECS=[('D_S1b_split_Ni','02_static/CONTCAR','9437467a79f82eae1992d71d77739051b644651d3ad2c7acd83acbcd7e98b049'),('A_H2_Mg_parallel','02_static/CONTCAR','760dc2996184545d11706b56ddddd5e85ee613b890d2445ef85374e0c7ae1b83'),('H2_reference','01_relax/POSCAR','07490a433b8716285d1a7213c649132ca1b07cde917160b9e84e4501bff0a06a')]
POT='72b607bd5434dcb71996bb13fd478af3592115762200d3c53e0f47c68b0a84ae'
HPOT='030f79b5d3ab3cf0e668861823c8fb652ff669f3e15e46930bd03bfd63a607b6'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(O/'mgni_cycle_batch.py')==PIN,'Source driver changed'
spec=importlib.util.spec_from_file_location('cycle_reader',O/'mgni_cycle_batch.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert not R.exists(),'Repair directory already exists; do not overwrite or resubmit'
q=subprocess.run(['squeue','-r','-h','-j','63565364_1,63565364_5,63565365,63565366','-o','%i'],capture_output=True,text=True,timeout=30)
assert not q.stdout.strip(),'Original target still active'
changes={'LREAL':'.FALSE.','EDIFF':'1E-7','EDIFFG':'-0.015','NSW':200,'POTIM':'0.10','ISTART':0,'ICHARG':2}
for name,rel,h in SPECS:
 p=O/name
 assert sha(p/rel)==h,'Changed restart geometry: '+name
 assert not (p/'VERIFIED_FINISHED').exists(),'Target already verified'
 for st in ['01_relax','02_static']:assert sha(p/st/'POTCAR')==(HPOT if name=='H2_reference' else POT)
 if name!='H2_reference':
  v=m.read_output(p/'02_static');assert v['last_scf_ediff'] and v['normal_end'] and .02<=v['max_free_force_eV_A']<.022
 else:assert 'oom_kill' in (p/'01_relax/vasp.stderr').read_text(errors='replace'),'Missing OOM evidence'
p=O/'D_S1a_same_Ni';v=m.read_output(p/'02_static')
assert (p/'VERIFIED_FINISHED').is_file() and v['last_scf_ediff'] and v['normal_end'] and v['max_free_force_eV_A']<.02
assert not (O/'D_S2_refill_one').exists() and not (O/'D_S3_refill_two').exists(),'Refill already has outputs'
R.mkdir();(R/'logs').mkdir();(R/'reference').mkdir()
for f in ['parent_D.vasp','energies.json']:shutil.copy2(O/'reference'/f,R/'reference'/f)
t=json.loads((O/'reference/templates.json').read_text())
for key in t:t[key]['01_relax']['incar']=m.patch(t[key]['01_relax']['incar'],changes)
(R/'reference/templates.json').write_text(json.dumps(t,indent=2))
for name,rel,h in SPECS:
 p=O/name;d=R/name;d.mkdir();shutil.copy2(p/'model.json',d/'model.json')
 for st in ['01_relax','02_static']:
  x=d/st;x.mkdir()
  for f in ['KPOINTS','POTCAR']:shutil.copy2(p/st/f,x/f)
  raw=(p/st/'INCAR').read_text()
  if st=='01_relax' and name!='H2_reference':raw=m.patch(raw,changes)
  (x/'INCAR').write_text(raw)
 shutil.copy2(p/rel,d/'01_relax/POSCAR')
 (d/'repair_provenance.json').write_text(json.dumps({'original':str(p),'restart':rel,'sha256':h,'changes':changes if name!='H2_reference' else {'memory':'4G to 16G'}},indent=2))
# Successful S1a is copied as a reference, not submitted again.
d=R/'D_S1a_same_Ni';d.mkdir();p=O/'D_S1a_same_Ni'
for f in ['model.json','VERIFIED_FINISHED']:shutil.copy2(p/f,d/f)
for st in ['01_relax','02_static']:
 (d/st).mkdir()
 for f in ['POSCAR','CONTCAR','INCAR','KPOINTS','POTCAR','OUTCAR','validation.json']:shutil.copy2(p/st/f,d/st/f)
(R/'input_manifest.json').write_text(json.dumps({str(p.relative_to(R)):sha(p) for p in R.rglob('*') if p.is_file()},indent=2))
print('PREPARED_ONLY|'+str(R),flush=True)
