#!/usr/bin/env python3
"""Resume only three unfinished LaH3 interface endpoints after 63636883 TIMEOUT.

- L_vLa continues from its original timed-out 01_relax/CONTCAR in a NEW directory.
- L_MgH2_vMg independently derives from the validated reference.
- L_transfer is created ONLY after the resumed L_vLa passes BOTH relax and static.
- Both previously finished L_vMg models are left untouched.
- No Hx/Mg/NEB or other original project job is submitted.
- Licensed POTCAR data stay in HPC; never published to GitHub.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HOME=Path.home()
SOURCE=HOME/'MgNi_LaH3_extensions_20261012'
STAGE=SOURCE/'La_resume_63636883_20261015'
RUNNER=SOURCE/'runner.py'
MODELS=('L_vLa','L_MgH2_vMg','L_transfer')
# This checksum was verified on the timed-out checkpoint from a read-only audit.
CHECKPOINT_SHA256='e6cde2a098bbe07aa0e18ab3fb267debbe7476eb3abe7aa04b8295a495b45da0'
NAME_TO_JOB={
 'L_vLa':'mechext_La_vLa_resume',
 'L_MgH2_vMg':'mechext_La_MgH2_vMg',
 'L_transfer':'mechext_La_transfer_resume',
}
PER_MODEL_HOURS={'L_vLa':24,'L_MgH2_vMg':16,'L_transfer':24}

def need(ok,msg):
 if not ok:raise RuntimeError('LA_RESUME_STOP|'+msg)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return datetime.datetime.now().astimezone().isoformat()
def save(path,obj):
 path=Path(path);tmp=path.with_name(path.name+'.tmp')
 tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n');tmp.replace(path)
def invoke(cmd,cwd=None,timeout=45):
 p=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,timeout=timeout)
 need(p.returncode==0,cmd[0]+' rc='+str(p.returncode)+': '+p.stderr[-350:])
 return p.stdout
def lib():
 m=json.loads((SOURCE/'manifest.json').read_text())
 need(m.get('preview') is False and (SOURCE/'submission.json').is_file(),'Original project is not verified production')
 need(sha(RUNNER)==m['input_hashes']['runner.py'],'Original production runner.py modified')
 need(sha(SOURCE/'plan.json')==m['input_hashes']['plan.json'],'Original physical parameters changed')
 spec=importlib.util.spec_from_file_location('original_La_extension',RUNNER)
 need(spec is not None and spec.loader is not None,'Unable to load original runner')
 module=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(module)
 need(module.ROOT==SOURCE,'Original runtime root mismatch')
 names=module.PLAN['La_children']
 need(names==['L_vMg_near','L_vMg_far','L_vLa','L_MgH2_vMg','L_transfer'],
      'Original La child roster changed')
 return module
def parent_valid(name,runner,root=SOURCE):
 d=root/name
 need((d/'VERIFIED_FINISHED').is_file(),'Parent or finished child missing validated marker: '+name)
 for stage,limit in (('01_relax',.015),('02_static',.02)):
  p=d/stage/'validation.json'
  need(p.is_file(),name+' missing validation at '+stage)
  info=json.loads(p.read_text())
  need(info.get('verified') and info.get('last_scf_ediff') and info.get('normal_end') and info['max_free_force_eV_A']<limit,
       name+' invalid physical result '+stage)
  if stage=='01_relax':need(info.get('ionic_converged'),name+' did not ionic converge')
 need((d/'02_static/CONTCAR').is_file(),name+' missing valid geometry')
 return runner.finished(name)
def check_slurm():
 output=invoke(['sacct','-X','-n','-P','-j','63636883','-o',
                'JobID,JobName%38,State,ExitCode,Elapsed,Timelimit,Start,End'],timeout=30)
 matches=[x.split('|') for x in output.splitlines() if x.split('|')[0]=='63636883']
 need(len(matches)==1 and matches[0][2]=='TIMEOUT','Original 63636883 is not confirmed TIMEOUT')
 old=json.loads((SOURCE/'recovery_La_Mg_only_20261013/submission.json').read_text())
 need(old.get('jobs',{}).get('La')=='63636883','Original La recovery ledger disagrees')
 return matches[0]
def check_active():
 p=subprocess.run(['squeue','-h','-u',os.environ.get('USER','deanfilan'),'-o','%i|%j|%T'],
                  capture_output=True,text=True,timeout=30)
 need(p.returncode==0,'Unable to verify existing Slurm queue: '+p.stderr[-350:])
 active=[line for line in p.stdout.splitlines()
         if any(('|'+v+'|') in ('|'+line+'|') or '|'+v+'|' in line for v in NAME_TO_JOB.values())]
 need(not active,'Already queued La branch continuation '+str(active[:6]))
def check_source():
 need(SOURCE.is_dir() and RUNNER.is_file(),'Original La study missing')
 state=check_slurm()
 runner=lib()
 for n in ('L_parent','L_MgH2_reference','L_vMg_near','L_vMg_far'):parent_valid(n,runner)
 for name in ('L_MgH2_vMg','L_transfer'):
  need(not (SOURCE/name).exists(),'Child unexpectedly already built in old project: '+name)
 p=SOURCE/'L_vLa'
 need(p.is_dir() and not (p/'VERIFIED_FINISHED').exists(),'La-vacancy already finished or not created')
 need(not (p/'02_static/validation.json').exists(),'Original L_vLa Static already exists')
 start=p/'01_relax/CONTCAR'
 need(start.is_file() and sha(start)==CHECKPOINT_SHA256,'Original timed-out CONTCAR changed, do not resume unknown structure')
 meta=json.loads((p/'model.json').read_text())
 need(meta['name']=='L_vLa' and meta['group']=='L','Wrong vacancy source metadata')
 need(meta['metadata']['deleted_site']['deleted_id']==55,'La vacancy source H55 mapping changed')
 s=runner.parse(start.read_text());s['ids']=meta['ids']
 need(len(s['species'])==len(meta['ids'])==99,'La vacancy resumed structure must contain 99 atoms')
 need({e:s['species'].count(e) for e in set(s['species'])}=={'Mg':12,'La':16,'H':71},
      'La-vacancy checkpoint wrong composition')
 need(len(set(s['ids']))==99 and 55 not in s['ids'],'Deleted-H mapping not preserved')
 for filename in ('INCAR','KPOINTS','POTCAR'):
  d=p/'01_relax'/filename
  need(d.is_file() and d.stat().st_size>({'KPOINTS':10,'INCAR':40,'POTCAR':100}[filename]),'Old checkpoint '+filename+' missing/empty')
 mg,_=runner.finished('L_MgH2_reference')
 need(39 in mg['ids'] and mg['species'][mg['ids'].index(39)]=='H',
      'MgH2 cleavage reference lost original H39')
 need(s['mask'].count([False,False,False])>0,'Unexpected fixed constraints in La checkpoint')
 return runner,s,meta,state
def generator_inputs():
 runner,s,meta,state=check_source()
 # Copy original potentials into new independent directory (licensed HPC-local copies).
 pots_l=runner.copied_pots('L_vLa')
 ref,refmeta=runner.finished('L_MgH2_reference')
 pots_ref=runner.copied_pots('L_MgH2_reference')
 need(set(pots_l)==set(s['species']) and set(pots_ref)==set(ref['species']),
      'PBE potential species mismatch')
 runner.ROOT=STAGE
 runner.write_model('L_vLa',s,'L',meta['metadata'],pots_l)
 old=SOURCE/'L_vLa'
 for phase in ('01_relax','02_static'):
  for fn in ('INCAR','KPOINTS','POTCAR'):
   need(sha(STAGE/'L_vLa'/phase/fn)==sha(old/phase/fn),
        'Changed original La-vacancy physical input '+phase+'/'+fn)
 v=runner.parse((STAGE/'L_vLa/01_relax/POSCAR').read_text())
 need(v['species']==s['species'] and v['mask']==s['mask'],'Restart changed atom roster')
 need(max(runner.distance(x,y,s['cell']) for x,y in zip(s['xyz'],v['xyz']))<1e-6,
      'Restart differs from timed-out CONTCAR')
 idx=ref['ids'].index(39)
 need(ref['species'][idx]=='H','H39 no longer hydrogen')
 vacancy={'deleted_id':39,'xyz':ref['xyz'][idx]}
 new=runner.subset(ref,[i for i in range(len(ref['ids'])) if i!=idx])
 need(new['species'].count('Mg')==12 and new['species'].count('H')==23,
      'MgH2 finite cleavage vacancy not Mg12H23')
 runner.write_model('L_MgH2_vMg',new,refmeta['group'],
                   {'parent':'L_MgH2_reference','deleted_site':vacancy},pots_ref)
 return {'source_checkpoint_SHA256':CHECKPOINT_SHA256,
         'original_job':'63636883','original_state':state,
         'created_initial_models':['L_vLa','L_MgH2_vMg'],
         'late_transfer_requires_verified_L_vLa':True}
def make_job(name):
 mode={'L_vLa':'--run-vLa','L_MgH2_vMg':'--run-MgH2','L_transfer':'--run-transfer'}[name]
 wall=PER_MODEL_HOURS[name]
 contents=f'''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --time={wall:02d}:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
python3 continue_La_63636883_20261015.py {mode}
'''
 p=STAGE/(name+'.slurm');p.write_text(contents)
 invoke(['bash','-n',str(p)])
 return p
def stage():
 need(not STAGE.exists(),'Continuation directory already exists; inspect existing ledger, do not overwrite')
 runner,s,meta,state=check_source()
 check_active()
 STAGE.mkdir(mode=0o700)
 (STAGE/'logs').mkdir()
 shutil.copy2(Path(__file__),STAGE/'continue_La_63636883_20261015.py')
 shutil.copy2(SOURCE/'plan.json',STAGE/'physical_plan.json')
 info=generator_inputs()
 for n in MODELS:make_job(n)
 hashes={str(p.relative_to(STAGE)):sha(p) for p in STAGE.rglob('*') if p.is_file()}
 save(STAGE/'manifest.json',{'created_at':now(),'scope':list(MODELS),
      'hashes':hashes,'source':info,'max_new_models':3,'no_NEB':True,
      'no_Hx_resume':True,'keeps_original_folders_untouched':True})
 print('STAGED|La vacancy from preserved checkpoint; MgH2 vacancy from converged parent; La transfer waits for vacancy',flush=True)
def submit():
 need(shutil.which('sbatch') is not None,'Not a Slurm login node')
 if STAGE.exists():
  p=STAGE/'submission.json'
  if p.is_file():print('EXISTING_SUBMISSION|'+p.read_text(),flush=True);return
  raise RuntimeError('Continuation root exists without submitted ledger; inspect manually, no duplicates')
 stage()
 manifest=json.loads((STAGE/'manifest.json').read_text())
 for name,value in manifest['hashes'].items():
  need(sha(STAGE/name)==value,'Stage input changed before submitting: '+name)
 (STAGE/'.attempted_submission').mkdir(exist_ok=False)
 ledger={'created_at':now(),'prior_job':'63636883','scope':list(MODELS),
         'jobs':{},'completed_submission':False}
 save(STAGE/'submission.json',ledger)
 for model in MODELS:
  args=['sbatch','--parsable','--job-name='+NAME_TO_JOB[model]]
  if model=='L_transfer':
   args.extend(['--dependency=afterok:'+ledger['jobs']['L_vLa'],'--kill-on-invalid-dep=yes'])
  args.append(model+'.slurm')
  raw=invoke(args,cwd=STAGE).strip()
  ident=raw.split(';')[0]
  need(re.fullmatch(r'[0-9]+',ident),'Ambiguous sbatch id for '+model)
  ledger['jobs'][model]=ident
  save(STAGE/'submission.json',ledger)
  print('SUBMITTED|'+model+'|'+ident+('|afterok='+ledger['jobs']['L_vLa'] if model=='L_transfer' else '|independent'),flush=True)
 ledger['completed_submission']=True;ledger['submitted_at']=now()
 save(STAGE/'submission.json',ledger)
 print('LA_CONTINUATION_SUBMITTED|'+json.dumps(ledger,ensure_ascii=False),flush=True)
def checked_stage(name):
 need((STAGE/name/'VERIFIED_FINISHED').exists(),'Parent resumed L_vLa not verified')
 runner=lib();runner.ROOT=STAGE
 s,meta=parent_valid(name,runner,root=STAGE)
 return runner,s,meta
def run_simple(name):
 need(os.environ.get('SLURM_NTASKS')=='64','Job does not have 64 MPI tasks')
 need(name in MODELS[:2],'Unexpected direct model')
 ledger=json.loads((STAGE/'submission.json').read_text())
 need(ledger['jobs'].get(name),'Not submitted by this La-continuation ledger')
 orig=lib();orig.ROOT=STAGE
 manifest=json.loads((STAGE/'manifest.json').read_text())
 need(sha(STAGE/'continue_La_63636883_20261015.py')==manifest['hashes']['continue_La_63636883_20261015.py'],
      'Runner changed after submission')
 for st in ('01_relax','02_static'):
  for fn in ('INCAR','KPOINTS','POTCAR'):
   path=name+'/'+st+'/'+fn
   need(sha(STAGE/path)==manifest['hashes'][path],'Input modified '+path)
 path=name+'/01_relax/POSCAR'
 need(sha(STAGE/path)==manifest['hashes'][path],'Initial checkpoint modified')
 print('RUNNING_MODEL|'+name+'|'+now(),flush=True)
 orig.run_model(name)
 need((STAGE/name/'VERIFIED_FINISHED').exists(),'Model did not pass Relax+Static validation')
 print('MODEL_VERIFIED_FINISHED|'+name+'|'+now(),flush=True)
def run_transfer():
 need(os.environ.get('SLURM_NTASKS')=='64','Job does not have 64 MPI tasks')
 ledger=json.loads((STAGE/'submission.json').read_text())
 need(ledger['jobs'].get('L_transfer') and ledger['jobs'].get('L_vLa'),'Missing transfer job dependency')
 if (STAGE/'L_transfer').exists():
  if (STAGE/'L_transfer/VERIFIED_FINISHED').exists():
   print('ALREADY_FINISHED|L_transfer');return
  raise RuntimeError('L_transfer exists incomplete; refuse duplicate/overwrite')
 runner,s,meta=checked_stage('L_vLa')
 site=meta['metadata']['deleted_site']
 need(site['deleted_id']==55,'Not original La-derived H vacancy')
 anchor=site['anchor_id'];need(anchor in s['ids'],'Original La anchor missing')
 a=s['xyz'][s['ids'].index(anchor)]
 target=[a[k]+site['vector'][k] for k in range(3)]
 cand=[]
 for donor in range(29,53):
  if donor not in s['ids']:continue
  j=s['ids'].index(donor)
  need(s['species'][j]=='H','Candidate donor mapping changed')
  hh=[runner.distance(target,x,s['cell']) for x,e,oid in zip(s['xyz'],s['species'],s['ids']) if e=='H' and oid!=donor]
  met=[runner.distance(target,x,s['cell']) for x,e in zip(s['xyz'],s['species']) if e!='H']
  if hh and met and min(hh)>1.10 and min(met)>1.25 and runner.neighbours(s,s['xyz'][j],'Mg')[0][0]<2.4:
   cand.append((runner.distance(target,s['xyz'][j],s['cell']),donor))
 need(cand,'No physically safe MgH2 donor; stop rather than force invalid transfer')
 travel,donor=min(cand)
 s['xyz'][s['ids'].index(donor)]=target
 need(s['species'].count('Mg')==12 and s['species'].count('La')==16 and s['species'].count('H')==71,
      'Transfer composition changed')
 metadata={'parent':'L_vLa','transferred_H_ids':[donor],
           'acceptor_deleted_H':55,'displacement_A':travel,'no_H_added':True}
 pots=runner.copied_pots('L_vLa')
 runner.write_model('L_transfer',s,'L',metadata,pots)
 print('TRANSFER_INITIALIZED|donor='+str(donor)+'|travel_A='+str(travel)+'|stoichiometry=Mg12La16H71',flush=True)
 runner.run_model('L_transfer')
 need((STAGE/'L_transfer/VERIFIED_FINISHED').exists(),'Transfer did not pass Relax+Static validation')
 print('MODEL_VERIFIED_FINISHED|L_transfer|'+now(),flush=True)
def status():
 print('CHECKED_AT|'+now(),flush=True)
 if not (STAGE/'submission.json').exists():print('NOT_SUBMITTED');return
 rec=json.loads((STAGE/'submission.json').read_text());print('LA_RESUME_LEDGER|'+json.dumps(rec,ensure_ascii=False))
 ids=','.join(rec['jobs'].values())
 for args in (['squeue','-r','-h','-j',ids,'-o','%i|%j|%T|%M|%R'],
              ['sacct','-X','-n','-P','-j',ids,'-o','JobID,JobName%38,State,ExitCode,Elapsed,Timelimit,Start,End']):
  p=subprocess.run(args,capture_output=True,text=True,timeout=30)
  print('SLURM|rc='+str(p.returncode)+'|'+p.stdout.strip().replace('\n',' || '),flush=True)
 for name in MODELS:
  m=STAGE/name
  d={'exists':m.exists(),'VERIFIED_FINISHED':(m/'VERIFIED_FINISHED').is_file()}
  for stage in ('01_relax','02_static'):
   p=m/stage/'validation.json'
   if p.is_file():
    v=json.loads(p.read_text())
    d[stage]={k:v.get(k) for k in ('verified','normal_end','last_scf_ediff','ionic_converged','max_free_force_eV_A','F_eV','E0_eV')}
   else:
    d[stage]={'OUTCAR_exists':(m/stage/'OUTCAR').exists(),'CONTCAR_exists':(m/stage/'CONTCAR').exists()}
  print('MODEL|'+name+'|'+json.dumps(d),flush=True)
def main():
 pa=argparse.ArgumentParser();g=pa.add_mutually_exclusive_group(required=True)
 for x in ('precheck','submit','run-vLa','run-MgH2','run-transfer','status'):
  g.add_argument('--'+x,action='store_true')
 a=pa.parse_args()
 if a.precheck:
  r,s,m,state=check_source()
  print('CONTINUATION_PRECHECK_OK|job=63636883 state='+state[2]+'|checkpoint_SHA256='+CHECKPOINT_SHA256+
        '|L_vLa_atoms='+str(len(s['species']))+'|two_prior_children_finished=True|no_Hx=True',flush=True)
 elif a.submit:submit()
 elif a.run_vLa:run_simple('L_vLa')
 elif a.run_MgH2:run_simple('L_MgH2_vMg')
 elif a.run_transfer:run_transfer()
 else:status()
if __name__=='__main__':
 try:main()
 except Exception as ex:
  print('LA_CONTINUATION_FAILED|'+str(ex),file=sys.stderr,flush=True)
  raise
