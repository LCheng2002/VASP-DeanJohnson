#!/usr/bin/env python3
"""Read only three failed MgNi calculations and their cancelled dependency."""
import hashlib, importlib.util, json, os, re, subprocess
from pathlib import Path
ROOT=Path.home()/'MgNi_abs_des_cycle_20261011'
NAMES=['D_S1b_split_Ni','A_H2_Mg_parallel','H2_reference']
def digest(p):
 if not p.is_file():return None
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p,limit=16000):
 if not p.is_file():return None
 with p.open('rb') as f:
  size=p.stat().st_size
  if size>limit:f.seek(size-limit)
  return f.read(limit).decode('utf-8',errors='replace')
def tail(p,n=30):
 t=read(p)
 return None if t is None else '\n'.join(t.splitlines()[-n:])
spec=importlib.util.spec_from_file_location('cycle',ROOT/'mgni_cycle_batch.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
report={'checked_at':subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip(),'root':str(ROOT),'driver_sha256':digest(ROOT/'mgni_cycle_batch.py'),'models':{},'logs':{}}
for name in NAMES:
 p=ROOT/name
 item={'model_json':read(p/'model.json'),'verified_finished':(p/'VERIFIED_FINISHED').is_file(),'stages':{}}
 for st in ['01_relax','02_static']:
  d=p/st;v={'files':{},'incar':read(d/'INCAR'),'kpoints':read(d/'KPOINTS'),'validation':None}
  for fn in ['POSCAR','CONTCAR','INCAR','KPOINTS','POTCAR','WAVECAR','CHGCAR','OUTCAR','OSZICAR','vasp.stdout','vasp.stderr']:
   f=d/fn
   if f.is_file():
    v['files'][fn]={'bytes':f.stat().st_size}
    if fn not in ['OUTCAR','OSZICAR','vasp.stdout','vasp.stderr','WAVECAR','CHGCAR']:v['files'][fn]['sha256']=digest(f)
  if (d/'validation.json').is_file():
   r=json.loads((d/'validation.json').read_text());v['validation']={k:x for k,x in r.items() if k not in ['H_nearest_metals','final_H_neighbours']}
  if (d/'OUTCAR').is_file():
   try:v['outcar_recheck']=m.read_output(d)
   except Exception as e:v['outcar_recheck_error']=str(e)
  v['stdout_tail']=tail(d/'vasp.stdout',25);v['stderr_tail']=tail(d/'vasp.stderr',25);v['outcar_tail']=tail(d/'OUTCAR',25);v['oszicar_tail']=tail(d/'OSZICAR',8)
  if name=='H2_reference':v['poscar']=read(d/'POSCAR');v['contcar']=read(d/'CONTCAR')
  item['stages'][st]=v
 report['models'][name]=item
for p in sorted((ROOT/'logs').glob('*')):
 if p.is_file() and any(s in p.name for s in ['63565364_1.','63565364_5.','63565365_','63565366_']):report['logs'][p.name]=tail(p,45)
for cmd,key in [(['squeue','-r','-u',os.environ['USER'],'-h','-o','%i|%j|%T|%M|%R'],'queue'),(['sacct','-X','-n','-P','-j','63565364,63565365,63565366','-o','JobID,JobName%28,State,ExitCode,Elapsed,Start,End'],'accounting')]:
 p=subprocess.run(cmd,capture_output=True,text=True,timeout=30);report[key]={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
report['S1a_validation']=json.loads((ROOT/'D_S1a_same_Ni/02_static/validation.json').read_text());report['S1a_validation'].pop('H_nearest_metals',None)
report['S2_exists']=(ROOT/'D_S2_refill_one').exists();report['S3_exists']=(ROOT/'D_S3_refill_two').exists()
print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
