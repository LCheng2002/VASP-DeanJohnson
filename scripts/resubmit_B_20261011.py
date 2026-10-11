#!/usr/bin/env python3
"""Authorized continuation of two failed VASP geometry relaxations only."""
import os,re,json,hashlib,shutil,subprocess,datetime
from pathlib import Path
H=Path.home();R=H/'Mg2Ni_Mg_H_thermo_validation_20261010';B=R/'retry_20261011_B_v1'
specs=[('B-H-Mg_Mg_Mg2NiH4','8478ae373a9b5df5f0afd916339e7df2583d207fbdb06adac5e02509dc450df9','mgni_Bretry_Mg',63535597),('B-H-Ni_Mg_Mg2NiH4','0d642160d1d57d76249def4a85a4cc64f6a242dc593222c1ee7d8d6f35d9bce0','mgni_Bretry_Ni',63535598)]
POT='72b607bd5434dcb71996bb13fd478af3592115762200d3c53e0f47c68b0a84ae'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def config(p):
 d={}
 for l in p.read_text().splitlines():
  for s in re.split(';',re.split('[!#]',l)[0]):
   if '=' in s:
    k,v=s.split('=',1);d[k.strip().upper()]=v.strip()
 return d
def patch(s,updates):
 lines=[]
 for line in s.splitlines():
  k=re.split('=',re.split('[!#]',line)[0],1)[0].strip().upper()
  if k not in updates:lines.append(line)
 return '\n'.join(lines)+'\n'+'\n'.join(k+' = '+v for k,v in updates.items())+'\n'
validator=r'''import re,json,math,sys
from pathlib import Path
p=Path('.');ls=(p/'POSCAR').read_text().splitlines();n=sum(map(int,ls[6].split()));i=7
select=ls[i].strip().lower().startswith('s')
if select:i+=1
i+=1
masks=[[t.upper().startswith('T') for t in l.split()[3:6]] if select else [True]*3 for l in ls[i:i+n]]
assert len(masks)==n and all(len(m)==3 for m in masks)
F=None;E0=None;ionic=False;footer=False;ediff=False;iteration=0;forces=[];buf=[];reading=False
with open('OUTCAR',errors='replace') as f:
 for line in f:
  m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
  if m and int(m.group(1))!=iteration:iteration=int(m.group(1));ediff=False
  if 'aborting loop because EDIFF is reached' in line:ediff=True
  if 'reached required accuracy' in line:ionic=True
  if 'General timing and accounting informations for this job' in line:footer=True
  if 'free  energy   TOTEN' in line:F=float(line.split('=')[1].split()[0])
  if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
  if 'TOTAL-FORCE (eV/Angst)' in line:buf=[];reading=True;continue
  if reading:
   z=line.split()
   if len(z)==6:
    try:buf.append([float(t) for t in z[3:6]]);continue
    except ValueError:pass
   if buf:forces=buf;reading=False
assert len(forces)==n
maxf=max(math.sqrt(sum(v*v for v,free in zip(vec,mask) if free)) for vec,mask in zip(forces,masks))
valid=footer and ediff and F is not None and maxf<0.0200001 and (ionic or sys.argv[1]=='static')
r={'phase':sys.argv[1],'F_eV':F,'E0_eV':E0,'last_scf_converged':ediff,'ionic_converged':ionic,'normal_end':footer,'max_free_force_eV_A':maxf,'verified':valid}
Path('validation.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
if not valid:raise SystemExit(40)
'''
job='''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# Refuse to execute this job again in an already-started directory.
mkdir .execution_started
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
(cd 01_relax && mpirun -np "$SLURM_NTASKS" vasp_std && python3 ../check_stage.py relax)
cp 01_relax/CONTCAR 02_static/POSCAR
(cd 02_static && mpirun -np "$SLURM_NTASKS" vasp_std && python3 ../check_stage.py static)
touch VERIFIED_FINISHED
'''
q=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%j'],capture_output=True,text=True,check=True,timeout=30).stdout.splitlines()
# Preflight both sources before any submission; no edits to the original models.
for name,chash,jobname,oldid in specs:
 src=R/name/'01_relax'
 assert sha(src/'CONTCAR')==chash,'Source CONTCAR changed: '+name
 assert sha(src/'POTCAR')==POT,'POTCAR changed: '+name
 ls=(src/'CONTCAR').read_text().splitlines();assert ls[5].split()==['Mg','Ni','H'] and list(map(int,ls[6].split()))==[44,8,33]
 cfg=config(src/'INCAR');assert cfg['NSW']=='250' and cfg['EDIFFG']=='-0.02' and cfg['POTIM']=='0.20' and cfg['IBRION']=='2'
 if ('mgni_'+name in q or jobname in q) and not (B/name/'submission.json').is_file():raise RuntimeError('Related job is already active: '+name)
 assert (R/'B0_Mg_Mg2NiH4_clean/02_static/INCAR').is_file()
B.mkdir(exist_ok=True)
results=[]
for name,chash,jobname,oldid in specs:
 dst=B/name;src=R/name/'01_relax'
 if (dst/'submission.json').is_file():
  results.append(json.loads((dst/'submission.json').read_text()));continue
 # A directory from an uncertain earlier attempt is never blindly resubmitted.
 dst.mkdir(exist_ok=False)
 for stage in ('01_relax','02_static'):
  d=dst/stage;d.mkdir()
  for fn in ('KPOINTS','POTCAR'):shutil.copy2(src/fn,d/fn)
 shutil.copy2(src/'CONTCAR',dst/'01_relax/POSCAR')
 updates={'NSW':'500','POTIM':'0.10','ISTART':'0','ICHARG':'2'}
 (dst/'01_relax/INCAR').write_text(patch((src/'INCAR').read_text(),updates))
 shutil.copy2(R/'B0_Mg_Mg2NiH4_clean/02_static/INCAR',dst/'02_static/INCAR')
 (dst/'check_stage.py').write_text(validator)
 (dst/'job.slurm').write_text(job)
 subprocess.run(['bash','-n',str(dst/'job.slurm')],check=True)
 compile(validator,str(dst/'check_stage.py'),'exec')
 provenance={'source_job':oldid,'source_model':name,'source_contcar_sha256':chash,'potcar_sha256':POT,'changes':updates,'n_tasks':64,'time_limit':'24:00:00','preserve_original':True,'directory':str(dst.relative_to(H))}
 (dst/'provenance.json').write_text(json.dumps(provenance,indent=2))
 # Persist before calling sbatch to protect against lost acknowledgements.
 (dst/'.submission_attempted').touch(exist_ok=False)
 proc=subprocess.run(['sbatch','--parsable','--job-name='+jobname,'job.slurm'],cwd=dst,capture_output=True,text=True,timeout=45)
 if proc.returncode!=0:
  (dst/'submission_error.txt').write_text(proc.stderr);raise RuntimeError('sbatch rejected '+name+': '+proc.stderr[:400])
 jid=proc.stdout.strip();assert re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?',jid),'Unexpected sbatch response'
 provenance['job_id']=jid;provenance['submitted_utc']=datetime.datetime.utcnow().isoformat()+'Z'
 (dst/'submission.json').write_text(json.dumps(provenance,indent=2));(dst/'.slurm_jobid').write_text(jid+'\n')
 print('SUBMITTED|'+name+'|'+jid,flush=True);results.append(provenance)
# Small, read-only TaskA provenance supplement for interpretation of the third H.
T=H/'TaskA_Mg17La2_H_enrichment_20261010'
meta={}
for model in ('H3_addLa','H3_addMg'):
 p=T/model/'model_metadata.tsv'
 if p.is_file():meta[model]=p.read_text()[:12000]
p=T/'parent_source/S+2H_La_static_OUTCAR';parent={}
if p.is_file():
 for line in p.open(errors='replace'):
  if 'free  energy   TOTEN' in line:parent['F_eV']=float(line.split('=')[1].split()[0])
  if 'energy(sigma->0)' in line:parent['E0_eV']=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
print('SUBMISSION_SUMMARY_JSON|'+json.dumps({'jobs':results,'TaskA_metadata':meta,'TaskA_parent_2H':parent},ensure_ascii=False),flush=True)
