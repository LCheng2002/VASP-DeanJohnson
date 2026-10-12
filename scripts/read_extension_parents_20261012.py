#!/usr/bin/env python3
"""Read only the fixed VASP inputs listed below. Never export a POTCAR dataset."""
import hashlib,json,subprocess
from pathlib import Path
home=Path.home()
paths={
 'C0':'Mg2Ni_Mg_H_thermo_validation_20261010/C0_MgH2_Mg2Ni_clean',
 'A0':'Mg2Ni_Mg_H_thermo_validation_20261010/A0_Mg_Mg2Ni_clean',
 'S1b':'MgNi_abs_des_cycle_20261011/repairs_20261012_v1/D_S1b_split_Ni',
 'S2':'MgNi_abs_des_cycle_20261011/repairs_20261012_v1/D_S2_refill_one',
 'S3':'MgNi_abs_des_cycle_20261011/repairs_20261012_v1/D_S3_refill_two',
 'Hrelay':'Mg_LaH3_Hrelay_endpoints_20261008/S0_base_static_retry_20261010',
 'Mg_bulk':'10_bulk_references_unified/01_Mg',
 'Mg2Ni_bulk':'Mg2Ni_H_transfer_bulk_screen_20261009_supplement_v2/supplement_20261009/02_Mg2Ni'
}
report={'models':{}}
for name,rel in paths.items():
 root=home/rel;data={'path':rel,'files':{},'stages':{}}
 for filename in ('model.json','model_metadata.tsv'):
  p=root/filename
  if p.is_file() and p.stat().st_size<60000:data['files'][filename]=p.read_text(errors='replace')
 for stage in ('01_relax','02_static','01_static'):
  q=root/stage
  if not (q/'POSCAR').is_file():continue
  obj={}
  for fn in ('POSCAR','CONTCAR','INCAR','KPOINTS','validation.json'):
   p=q/fn
   if p.is_file() and p.stat().st_size<80000:
    obj[fn]=p.read_text(errors='replace');obj[fn+'_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
  p=q/'POTCAR'
  if p.is_file():
   obj['potential_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
   with p.open(errors='replace') as f:obj['potential_titles']=[l.strip() for l in f if 'TITEL' in l]
  p=q/'OUTCAR'
  if p.is_file():
   out={}
   with p.open(errors='replace') as f:
    for line in f:
     if 'free  energy   TOTEN' in line:out['F_eV']=float(line.split('=')[1].split()[0])
     if 'reached required accuracy' in line:out['ionic_converged']=True
     if 'General timing and accounting informations' in line:out['normal_end']=True
   obj['output_summary']=out
  data['stages'][stage]=obj
 report['models'][name]=data
report['checked_at']=subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()
print(json.dumps(report,ensure_ascii=False))
