#!/usr/bin/env python3
"""Read-only fixed-path NEB preflight. No VASP execution or pseudopotential export."""
import hashlib, importlib.util, json, os, re, shutil, subprocess
from pathlib import Path
H=Path.home()
paths={
 'Ni_clean':'Mg2Ni_Mg_H_thermo_validation_20261010/A0_Mg_Mg2Ni_clean/02_static',
 'Ni_molecule':'MgNi_abs_des_cycle_20261011/A_H2_Ni_vertical/02_static',
 'Ni_dissociated':'MgNi_abs_des_cycle_20261011/A_2H_Ni/02_static',
 'Ni_dissociated_parallel':'MgNi_abs_des_cycle_20261011/A_H2_Ni_parallel/02_static',
 'Mg_clean':'MgNi_LaH3_extensions_20261012/Mg_parent/02_static',
 'Mg_molecule':'MgNi_LaH3_extensions_20261012/Mg_H2/02_static',
 'Mg_dissociated':'MgNi_LaH3_extensions_20261012/Mg_2H/02_static',
 'H2':'MgNi_abs_des_cycle_20261011/repairs_20261012_v1/H2_reference/02_static',
 'Mg2Ni_bulk':'Mg2Ni_H_transfer_bulk_screen_20261009_supplement_v2/supplement_20261009/02_Mg2Ni/02_static'
}
r={'models':{},'environment':{},'queue':''}
for name,rel in paths.items():
 p=H/rel;d={'path':rel,'exists':p.is_dir(),'files':{}}
 for fn in ('POSCAR','CONTCAR','INCAR','KPOINTS','validation.json'):
  f=p/fn
  if f.is_file() and 0<f.stat().st_size<100000:
   d['files'][fn]={'text':f.read_text(errors='replace'),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
 f=p/'POTCAR'
 if f.is_file():
  d['potcar_sha256']=hashlib.sha256(f.read_bytes()).hexdigest()
  with f.open(errors='replace') as inp:d['potcar_titles']=[l.strip() for l in inp if 'TITEL' in l]
 f=p/'OUTCAR'
 if f.is_file():
  v={};it=None
  with f.open(errors='replace') as inp:
   for l in inp:
    m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',l)
    if m and m.group(1)!=it:it=m.group(1);v['last_ediff']=False
    if 'aborting loop because EDIFF is reached' in l:v['last_ediff']=True
    if 'free  energy   TOTEN' in l:v['F_eV']=float(l.split('=')[1].split()[0])
    if 'energy(sigma->0)' in l:v['E0_eV']=float(l.split('energy(sigma->0)')[1].split('=')[1].split()[0])
    if 'General timing and accounting informations' in l:v['normal_end']=True
  d['outcar_summary']=v
 r['models'][name]=d
exe=shutil.which('vasp_std');r['environment']['vasp_executable']=exe
if exe:
 p=Path(exe).resolve();r['environment']['vasp_resolved']=str(p)
 # Return booleans/hash only, not any licensed binary or source data.
 if p.is_file() and p.stat().st_size<400000000:
  b=p.read_bytes();r['environment']['vasp_sha256']=hashlib.sha256(b).hexdigest()
  r['environment']['feature_strings']={t.decode():t in b for t in (b'LCLIMB',b'ICHAIN',b'IOPT',b'VTST',b'CHAIN:')}
for x in ('numpy','scipy','ase'):
 try:r['environment'][x]=bool(importlib.util.find_spec(x))
 except Exception:r['environment'][x]=False
for key,cmd in [('module_info',['bash','-lc','module -t avail vasp 2>&1 | head -n 35']),('queue',['squeue','-r','-h','-j','63585911,63585912,63585913,63585914,63585915','-o','%i|%j|%T|%M|%R'])]:
 p=subprocess.run(cmd,capture_output=True,text=True,timeout=30);r[key]=p.stdout[:7000];r[key+'_returncode']=p.returncode
r['checked_at']=subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()
print(json.dumps(r,ensure_ascii=False))
