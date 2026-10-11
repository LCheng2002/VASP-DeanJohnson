#!/usr/bin/env python3
"""Read only selected modelling parents; never export POTCAR contents or credentials."""
from pathlib import Path
import json,hashlib,re,subprocess,os,datetime
H=Path.home();R=H/'Mg2Ni_Mg_H_thermo_validation_20261010'
names=['A0_Mg_Mg2Ni_clean','A-H-Mg_Mg_Mg2Ni','A-H-Ni_Mg_Mg2Ni','D0_MgH2_Mg2NiH4_clean','D-VH-Ni_MgH2_Mg2NiH4']
r={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':{},'project_scripts':{}}
for name in names:
 root=R/name;d={'files':{},'potcar':{},'energy_checks':{}}
 for rel in ['metadata.json','01_relax/POSCAR','01_relax/CONTCAR','01_relax/INCAR','01_relax/KPOINTS','02_static/POSCAR','02_static/CONTCAR','02_static/INCAR','02_static/KPOINTS']:
  p=root/rel
  if p.is_file() and p.stat().st_size<50000:
   b=p.read_bytes();d['files'][rel]={'text':b.decode('utf-8'),'sha256':hashlib.sha256(b).hexdigest()}
 for stage in ['01_relax','02_static']:
  p=root/stage/'POTCAR'
  if p.is_file():
   b=p.read_bytes();d['potcar'][stage]={'sha256':hashlib.sha256(b).hexdigest(),'titles':[x.strip() for x in b.decode(errors='replace').splitlines() if 'TITEL' in x]}
  p=root/stage/'OUTCAR';q={'last_toten':None,'last_e0':None,'normal_end':False,'ionic_converged':False,'last_scf_ediff':False};it=None
  if p.is_file():
   for l in p.open(errors='replace'):
    m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',l)
    if m and it!=m.group(1):it=m.group(1);q['last_scf_ediff']=False
    if 'free  energy   TOTEN' in l:q['last_toten']=float(l.split('=')[1].split()[0])
    if 'energy(sigma->0)' in l:q['last_e0']=float(l.split('energy(sigma->0)')[1].split('=')[1].split()[0])
    if 'reached required accuracy' in l:q['ionic_converged']=True
    if 'aborting loop because EDIFF is reached' in l:q['last_scf_ediff']=True
    if 'General timing and accounting informations for this job' in l:q['normal_end']=True
  d['energy_checks'][stage]=q
 r['models'][name]=d
for name in ['run_single.slurm','README.md']:
 p=R/name
 if p.is_file() and p.stat().st_size<80000:r['project_scripts'][name]=p.read_text()
x=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%i|%j|%T|%M|%R'],capture_output=True,text=True,timeout=30,check=True)
r['queue']=x.stdout
print(json.dumps(r,ensure_ascii=False))
