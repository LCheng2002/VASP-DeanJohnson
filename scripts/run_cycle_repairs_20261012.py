#!/usr/bin/env python3
"""Scheduled execution only; never submits jobs or changes original outputs."""
import argparse, hashlib, importlib.util, itertools, json, math, os, shutil, subprocess
from pathlib import Path
O=Path.home()/'MgNi_abs_des_cycle_20261011';R=O/'repairs_20261012_v1'
PIN='f5696bda76eb2c990e5daab3a09baabf8619a914bf5bed10b37a2a42aa495abc'
NAMES=['D_S1b_split_Ni','A_H2_Mg_parallel','H2_reference']
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def save(path,data):
 tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(data,indent=2));tmp.replace(path)
def engine():
 path=O/'mgni_cycle_batch.py'
 require(hashlib.sha256(path.read_bytes()).hexdigest()==PIN,'Original engine changed')
 spec=importlib.util.spec_from_file_location('verified_cycle_functions',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 # Only this process uses the isolated destination. No original script is edited.
 m.ROOT=R
 return m
def run(m,name):
 p=R/name;require(p.is_dir(),'Model missing')
 if (p/'VERIFIED_FINISHED').exists():return
 require(int(os.environ.get('SLURM_NTASKS','0'))==(4 if name=='H2_reference' else 64),'Wrong or absent Slurm allocation')
 (p/'.execution_lock').mkdir(exist_ok=False)
 expected=m.parse((p/'01_relax/POSCAR').read_text())
 for st in ['01_relax','02_static']:
  d=p/st
  if st=='02_static':shutil.copy2(p/'01_relax/CONTCAR',d/'POSCAR')
  with (d/'vasp.stdout').open('w') as out,(d/'vasp.stderr').open('w') as err:
   proc=subprocess.run(['mpirun','-np',os.environ['SLURM_NTASKS'],'vasp_std'],cwd=d,stdout=out,stderr=err)
  if proc.returncode:
   save(d/'execution_failure.json',{'returncode':proc.returncode,'phase':st})
   raise RuntimeError('VASP_EXECUTION_FAILED|'+name+'/'+st)
  v=m.read_output(d);s=m.parse((d/'CONTCAR').read_text())
  require(s['species']==expected['species'] and s['mask']==expected['mask'],'Composition or constraints changed')
  require(max(abs(x-y) for a,b in zip(s['cell'],expected['cell']) for x,y in zip(a,b))<1e-7,'Cell changed')
  limit=.005 if name=='H2_reference' else (.015 if st=='01_relax' else .02)
  v['force_gate_eV_A']=limit;v['phase']=st
  v['verified']=bool(v['normal_end'] and v['last_scf_ediff'] and all(v[k] is not None and math.isfinite(v[k]) for k in ['F_eV','E0_eV','max_free_force_eV_A']) and v['max_free_force_eV_A']<limit and (st=='02_static' or v['ionic_converged']))
  hs=[i for i,e in enumerate(s['species']) if e=='H']
  v['min_H_H_A']=min(m.distance(s['xyz'][i],s['xyz'][j],s['cell']) for i,j in itertools.combinations(hs,2)) if len(hs)>1 else None
  if name=='H2_reference':v['verified']=v['verified'] and .65<v['min_H_H_A']<.9
  save(d/'validation.json',v)
  save(d/'geometry_summary.json',{'H_nearest_metals':[{'atom':i+1,'Mg':m.nearest(s,s['xyz'][i],'Mg')[:2],'Ni':m.nearest(s,s['xyz'][i],'Ni')[:2]} for i in hs]})
  print('REPAIR_STAGE|'+name+'|'+json.dumps(v),flush=True)
  require(v['verified'],'QUALITY_GATE_FAILED|'+name+'/'+st)
 (p/'VERIFIED_FINISHED').touch()
def refill(m):
 require(int(os.environ.get('SLURM_NTASKS','0'))==64,'No intended Slurm allocation')
 if all((R/n/'VERIFIED_FINISHED').exists() for n in ['D_S2_refill_one','D_S3_refill_two']):return
 (R/'.refill_execution_started').mkdir(exist_ok=False)
 for n in ['D_S1a_same_Ni','D_S1b_split_Ni']:
  require((R/n/'VERIFIED_FINISHED').exists(),'Parent unverified: '+n)
  v=m.read_output(R/n/'02_static')
  require(v['normal_end'] and v['last_scf_ediff'] and v['max_free_force_eV_A']<.02,'Parent force recheck failed: '+n)
 # Keep the original atom mapping, donor criteria, collision checks and S2-before-S3 gate.
 m.run_model=lambda name:run(m,name)
 m.refill()
def status():
 print('CHECKED_AT|'+subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip(),flush=True)
 f=R/'submission.json'
 if not f.is_file():print('NO_SUBMISSION_RECORD');return
 r=json.loads(f.read_text());print('REPAIR_SUBMISSION|'+json.dumps(r),flush=True)
 ids=','.join(r['jobs'].values())
 for cmd in [['squeue','-r','-h','-j',ids,'-o','%i|%j|%T|%M|%R'],['sacct','-X','-n','-P','-j',ids,'-o','JobID,JobName%26,State,ExitCode,Elapsed,Start,End']]:
  q=subprocess.run(cmd,capture_output=True,text=True,timeout=30);print(q.stdout,flush=True)
 for n in NAMES+['D_S2_refill_one','D_S3_refill_two']:
  d=R/n;v={}
  for st in ['01_relax','02_static']:
   f=d/st/'validation.json'
   if f.is_file():v[st]=json.loads(f.read_text())
  print('REPAIR_MODEL|'+n+'|'+('VERIFIED' if (d/'VERIFIED_FINISHED').exists() else ('PRESENT' if d.exists() else 'NOT_CONSTRUCTED'))+'|'+json.dumps(v),flush=True)
 for f in sorted((R/'logs').glob('*.err')):
  if f.stat().st_size and f.stat().st_size<100000:print('ERROR_LOG|'+f.name+'|'+f.read_text(errors='replace')[-2000:].replace('\n',' | '),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['0','1','2','refill','status']);a=p.parse_args()
 if a.mode=='status':status()
 elif a.mode=='refill':refill(engine())
 else:run(engine(),NAMES[int(a.mode)])
