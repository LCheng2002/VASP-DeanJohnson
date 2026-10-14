#!/usr/bin/env python3
"""Bounded read-only audit of Ca corrosion models; no VASP or Slurm submission.
Export whitelisted text geometry/parameters and convergence metadata only.
Never export POTCAR datasets, keys, or unrelated projects.
"""
import os,re,json,time,hashlib,subprocess,datetime
from pathlib import Path
HOME=Path.home();start=time.monotonic()
report={'checked_at':datetime.datetime.now().astimezone().isoformat(),'jobs':{},'roots':[],'models':{},'scan':{}}
ids=['63131486','63306398']
p=subprocess.run(['sacct','-X','-n','-P','-j',','.join(ids),'-o','JobID,JobName%50,State,ExitCode,Start,End,WorkDir%250'],capture_output=True,text=True,timeout=25)
report['jobs']={'rc':p.returncode,'records':p.stdout,'stderr':p.stderr[:500]}
roots=[]
for line in p.stdout.splitlines():
 parts=line.split('|')
 if len(parts)>=7:
  q=Path(parts[6])
  if str(q).startswith(str(HOME)+'/') and q.is_dir():roots.append(q)
for q in sorted(HOME.iterdir()):
 if q.is_dir() and not q.is_symlink() and not q.name.startswith('.') and re.search(r'Ca|ca[_-]|63131486|[Cc]orros|OH|H2O',q.name):roots.append(q)
roots=list(dict.fromkeys(roots));report['roots']=[str(q.relative_to(HOME)) for q in roots]
pattern=re.compile(r'^(?:M[01](?:_|-).*(?:OH|H2O)|M[01]_(?:pristine|Ca_surface)|H2O(?:_reference)?|H2(?:_reference)?|water|hydrogen)$',re.I)
found=[];visited=0
for root in roots:
 for b,dirs,files in os.walk(root,followlinks=False):
  q=Path(b);visited+=1;depth=len(q.relative_to(root).parts)
  dirs[:]=[] if depth>=7 or visited>2600 or time.monotonic()-start>40 else [x for x in dirs if not x.startswith('.') and x not in ('__pycache__','node_modules')]
  if pattern.match(q.name):
   if any((q/s/'POSCAR').is_file() for s in ('.','01_relax','02_static','relax','static')):found.append(q)
  if visited>2600 or time.monotonic()-start>40:break
 if visited>2600 or time.monotonic()-start>40:break
report['scan']={'visited':visited,'bounded':visited>2600 or time.monotonic()-start>40,'found':[str(p.relative_to(HOME)) for p in found]}
def text(p,n=80000):
 if not p.is_file() or p.stat().st_size>n:return None
 return p.read_text(errors='replace')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def summary(p):
 if not p.is_file():return None
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-3000000));s=f.read().decode(errors='replace')
 energies=re.findall(r'free\s+energy\s+TOTEN\s*=\s*([-+0-9.Ee]+)',s)
 e0=re.findall(r'energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)',s)
 it=list(re.finditer(r'Iteration\s+(\d+)\(\s*(\d+)\)',s));lastcycle=s
 if it:
  last=it[-1].group(1);idx=next((m.start() for m in reversed(it) if m.group(1)!=last),None)
  if idx is not None:lastcycle=s[idx:]
 lastF=[];parts=s.rsplit('TOTAL-FORCE (eV/Angst)',1)
 if len(parts)==2:
  for l in parts[1].splitlines():
   cols=l.split()
   if len(cols)==6:
    try:lastF.append([float(v) for v in cols[3:6]]);continue
    except ValueError:pass
   if lastF:break
 return {'F_eV':float(energies[-1]) if energies else None,'E0_eV':float(e0[-1]) if e0 else None,'normal_end':'General timing and accounting informations' in s,'ionic_converged':'reached required accuracy' in s,'last_scf_ediff_seen':'aborting loop because EDIFF is reached' in lastcycle,'last_force_vectors':lastF,'last_ionic_iteration':int(it[-1].group(1)) if it else None,'OUTCAR_bytes':p.stat().st_size}
# Prefer OH parents, then existing water endpoints and clean parents. Cap reads.
found=sorted(set(found),key=lambda q:(0 if re.fullmatch(r'M[01]_OH(?:_H2O)?',q.name) else 1,str(q)))
for q in found[:45]:
 d={'path':str(q.relative_to(HOME)),'metadata':{},'stages':{}}
 for fn in ('model.json','metadata.json','model_info.txt','model_metadata.tsv'):
  t=text(q/fn)
  if t is not None:d['metadata'][fn]=t
 for stage in ('01_relax','02_static','relax','static','.'):
  sd=q/stage
  if not (sd/'POSCAR').is_file():continue
  v={}
  for fn in ('POSCAR','CONTCAR','INCAR','KPOINTS','validation.json'):
   t=text(sd/fn)
   if t is not None:v[fn]={'text':t,'sha256':digest(sd/fn)}
  pot=sd/'POTCAR'
  if pot.is_file():
   v['POTCAR_sha256']=digest(pot)
   with pot.open(errors='replace') as f:v['POTCAR_titles']=[l.strip() for l in f if 'TITEL' in l]
  v['outcar']=summary(sd/'OUTCAR')
  oz=sd/'OSZICAR'
  if oz.is_file():
   with oz.open('rb') as f:f.seek(max(0,oz.stat().st_size-1500));v['OSZICAR_tail']=f.read().decode(errors='replace')
  d['stages'][stage]=v
 report['models'][str(q.relative_to(HOME))]=d
# Candidate exact batch existence only, never inspect unrelated output directories.
report['new_batch_exists']=(HOME/'Ca_corrosion_OH_coverage_20261014').exists()
print(json.dumps(report,ensure_ascii=False))
