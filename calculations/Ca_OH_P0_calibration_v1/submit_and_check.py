#!/usr/bin/env python3
"""Ca OH P0: four static dipole recalibrations from archived relaxed geometries.

Source geometry: versioned, audited nonlicensed files already in this repository.
Potential data: exact old PAW datasets, acquired ONLY inside user's HPC account.
No continuation to P1/P2, no relaxed parent re-runs, no cancelled CaFresh restarts.
"""
import argparse, datetime, hashlib, json, math, os, re, shutil, subprocess, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
PLAN=json.loads((HERE/'plan.json').read_text())
ROOT=Path.home()/PLAN['production_directory']
NAMES=list(PLAN['models'])
ARCH=STAGE/PLAN['archive_source']

def need(ok,msg):
    if not ok:raise RuntimeError('PRECHECK_FAILED|'+msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now().astimezone().isoformat()
def save(p,value):
    p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n');tmp.replace(p)
def call(args,**kwargs):
    p=subprocess.run(args,text=True,capture_output=True,timeout=45,**kwargs)
    need(p.returncode==0,'command failed '+str(args[0])+': '+p.stderr[-280:])
    return p.stdout
def incar_dict(text):
    v={}
    for l in text.splitlines():
        l=re.split('[!#]',l)[0]
        if '=' in l:
            a,b=l.split('=',1);v[a.strip().upper()]=b.strip()
    return v
def poscar_info(text):
    ls=text.splitlines();need(len(ls)>9,'Empty POSCAR')
    scale=float(ls[1]);need(scale>0,'Bad scale')
    cell=[[float(t)*scale for t in x.split()[:3]] for x in ls[2:5]]
    els=ls[5].split();cnt=list(map(int,ls[6].split()))
    need(len(els)==len(cnt),'Elements/counts mismatch')
    sp=sum(([e]*n for e,n in zip(els,cnt)),[])
    i=7;sel=ls[i].lower().startswith('s');i+=int(sel)
    direct=ls[i].lower().startswith('d');i+=1
    rows=[x.split() for x in ls[i:i+len(sp)]]
    need(len(rows)==len(sp),'Missing atomic coordinate rows')
    masks=[[(x.upper()=='T') for x in r[3:6]] if sel else [True]*3 for r in rows]
    need(all(len(x)==3 for x in masks),'Bad selective dynamics flags')
    comp={e:sp.count(e) for e in els}
    return {'lines':ls,'cell':cell,'species':sp,'masks':masks,'coord_end':i+len(sp),'composition':comp,'direct':direct,'selective':sel}
def build_inputs():
    need(set(NAMES)=={'M0_OH','M1_OH','M0_OH_H2O','M1_OH_H2O'},'Unexpected case roster')
    out={}
    for name,entry in PLAN['models'].items():
        d=ARCH/name
        con=d/'results/01_relax/CONTCAR'
        inc=d/'calculation_inputs/INCAR.static'
        k=d/'calculation_inputs/KPOINTS'
        for f in (con,inc,k):need(f.is_file(),'Missing archived file '+str(f))
        need(sha(con)==entry['parent_con_sha256'],'Relax CONTCAR changed: '+name)
        need(sha(inc)==entry['parent_incar_sha256'],'Original INCAR changed: '+name)
        need(sha(k)==PLAN['common_KPOINTS_sha256'],'KPOINTS changed: '+name)
        s=poscar_info(con.read_text())
        need(s['composition']==entry['composition'],'Composition incorrect '+name)
        need(s['selective'] and sum(not any(x) for x in s['masks'])==32,'Changed bottom constraints '+name)
        cell=s['cell'];need(30<math.sqrt(sum(z*z for z in cell[2]))<40,'Incorrect slab vacuum/cell')
        raw=inc.read_text()
        tags=incar_dict(raw)
        need(tags.get('ENCUT')=='520' and tags.get('ISMEAR')=='1' and tags.get('SIGMA')=='0.2','Old physical settings mismatch '+name)
        # Keep ALL original INCAR tags apart from explicitly recalibrated electrostatics
        # and explicit defaults, to isolate the dipole variable.
        changes={'SYSTEM':name+'_P0_dipole_calibration',
                 'EDIFF':'1E-7','LREAL':'.FALSE.',
                 'LDIPOL':'.TRUE.','IDIPOL':'3',
                 'DIPOL':'0.500000000000 0.500000000000 %.12f'%entry['DIPOL_fractional_z'],
                 'ISPIN':'1','ISTART':'0','ICHARG':'2','IBRION':'-1','NSW':'0'}
        edited=[line for line in raw.splitlines() if re.split('[!#]',line)[0].split('=',1)[0].strip().upper() not in changes]
        edited.extend(k+' = '+v for k,v in changes.items())
        # Isolate original positions/constraints: VASP only reads the first atom-coordinate block.
        geom='\n'.join(s['lines'][:s['coord_end']])+'\n'
        out[name]={'POSCAR':geom,'INCAR':'\n'.join(edited)+'\n','KPOINTS':k.read_text(),
                   'composition':s['composition'],'old_F_eV':entry['old_F'],'old_E0_eV':entry['old_E0']}
    return out

def split_potcar(raw):
    parts=[];buf=[]
    for line in raw.splitlines(keepends=True):
        buf.append(line)
        if line.strip()==b'End of Dataset':parts.append(b''.join(buf));buf=[]
    need(parts and not b''.join(buf).strip(),'Not recognized POTCAR')
    out={}
    for chunk in parts:
        names=re.findall(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',chunk)
        need(len(names)==1,'Missing PAW label')
        title=names[0].decode()
        need(title not in out,'Duplicate PAW dataset')
        out[title]=chunk
    return out

def potentials():
    want=PLAN['required_potential_hashes']
    found={}
    seen_paths=set()
    def load(path):
        if not path.is_file() or str(path) in seen_paths or path.stat().st_size>7000000:return
        seen_paths.add(str(path))
        for label,b in split_potcar(path.read_bytes()).items():
            if label in want and sha_bytes(b)==want[label]:
                found[label]=b
    def sha_bytes(b):return hashlib.sha256(b).hexdigest()
    # First prefer the verified, locally retained Ca license cache.
    cache=Path.home()/'VASP/potential_sets/Ca_corrosion_fresh_20261014'
    for label in want:load(cache/label/'POTCAR')
    # If exact old datasets are different, audit only bounded prior Ca job roots.
    if set(found)!=set(want):
        roots=[Path.home()/'next_stage_63131486_20261002']
        r=subprocess.run(['sacct','-X','-n','-P','-j','63131486','-o','WorkDir%300'],
                         capture_output=True,text=True,timeout=20)
        if r.returncode==0:
            for line in r.stdout.splitlines():
                t=line.split('|')[0].strip()
                try:p=Path(t).resolve();p.relative_to(Path.home().resolve())
                except (ValueError,OSError):continue
                if p.is_dir():roots.append(p);roots.append(p.parent)
        visited=0;files=0;start=time.monotonic()
        for root in dict.fromkeys(roots):
            if not root.is_dir():continue
            for base,dirs,names in os.walk(root,followlinks=False):
                visited+=1;q=Path(base);depth=len(q.relative_to(root).parts)
                dirs[:]=[] if depth>=6 else [x for x in dirs if not x.startswith('.') and x!='__pycache__']
                if 'POTCAR' in names:load(q/'POTCAR');files+=1
                if set(found)==set(want) or visited>=800 or files>=100 or time.monotonic()-start>20:break
            if set(found)==set(want) or visited>=800 or files>=100 or time.monotonic()-start>20:break
    print('POTENTIAL_CHECK|available_exact_labels='+','.join(sorted(found))+
          '|matching_all='+str(set(found)==set(want)),flush=True)
    need(set(found)==set(want),'Missing the EXACT original PAW blocks; no submission and no data substituted')
    results={}
    for name,entry in PLAN['models'].items():
        b=b''.join(found[t] for t in entry['PAW_labels'])
        need(sha_bytes(b)==entry['combined_POTCAR_sha256'],'Combined POTCAR mismatch '+name)
        results[name]=b
    return results

def prepare_and_submit():
    # Always verify archived, nonlicensed inputs BEFORE touching cluster files.
    inputs=build_inputs()
    if (ROOT/'submission.json').is_file():
        print('EXISTING_SUBMISSION|'+(ROOT/'submission.json').read_text(),flush=True);return
    need(not ROOT.exists(),'Production directory already exists; refuse unsafe overwrite')
    pots=potentials()  # No POTCAR sent to GitHub or in logs.
    ROOT.mkdir(mode=0o700);(ROOT/'logs').mkdir();(ROOT/'models').mkdir()
    try:
        for name,data in inputs.items():
            d=ROOT/'models'/name;d.mkdir()
            for f in ('POSCAR','INCAR','KPOINTS'):(d/f).write_text(data[f])
            (d/'POTCAR').write_bytes(pots[name]);(d/'POTCAR').chmod(0o600)
            save(d/'metadata.json',{'composition':data['composition'],
                 'old_F_eV':data['old_F_eV'],'old_E0_eV':data['old_E0_eV'],
                 'parent_con_sha256':PLAN['models'][name]['parent_con_sha256'],
                 'nominal_dipol':PLAN['models'][name]['DIPOL_fractional_z']})
        shutil.copy2(Path(__file__),ROOT/'submit_and_check.py')
        shutil.copy2(HERE/'plan.json',ROOT/'plan.json')
        slurm='''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
ulimit -s unlimited
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 submit_and_check.py --run-index "$SLURM_ARRAY_TASK_ID"
'''
        (ROOT/'p0.slurm').write_text(slurm)
        call(['bash','-n',str(ROOT/'p0.slurm')])
        pinned={str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file()}
        save(ROOT/'production_manifest.json',{'input_sha256':pinned,'models':NAMES,'scope':'FOUR_STATIC_DIPOL_ONLY'})
        (ROOT/'.submission_attempted').mkdir()
        receipt=call(['sbatch','--parsable','--job-name=CaOH_P0_dipole','p0.slurm'],cwd=ROOT)
        job_id=receipt.strip().split(';')[0]
        need(re.fullmatch(r'\d+',job_id),'Ambiguous sbatch receipt; inspect before any retry')
        save(ROOT/'submission.json',{'job_id':job_id,'array':'0-3','models':NAMES,
             'submitted_at':now(),'CaFresh_resubmitted':False,'P1_P2_submitted':False})
        print('SUBMITTED|CaOH_P0|'+job_id+'|array=0-3',flush=True)
    except Exception:
        print('STOP|Production directory or attempt lock may exist. Inspect before retry. No automatic cleanup.',flush=True)
        raise

def read_results(name):
    d=ROOT/'models'/name
    geom=poscar_info((d/'POSCAR').read_text());N=len(geom['species'])
    fmax=None;F=E0=None;normal=scf=False;cyc=None;fs=[];reading=False;buffer=[]
    dipflag=False;idflag=False;spinflag=False
    with (d/'OUTCAR').open(errors='replace') as r:
        for line in r:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and m.group(1)!=cyc:cyc=m.group(1);scf=False
            if 'aborting loop because EDIFF is reached' in line:scf=True
            if 'General timing and accounting informations' in line:normal=True
            if 'free  energy   TOTEN' in line:F=float(line.split('=')[1].split()[0])
            if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if re.search(r'\bLDIPOL\s*=\s*T\b',line):dipflag=True
            if re.search(r'\bIDIPOL\s*=\s*3\b',line):idflag=True
            if re.search(r'\bISPIN\s*=\s*1\b',line):spinflag=True
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buffer=[];continue
            if reading:
                q=line.split()
                if len(q)==6:
                    try:buffer.append(list(map(float,q)));continue
                    except ValueError:pass
                if buffer:fs=buffer;reading=False
    need(len(fs)==N,'Missing complete VASP final force table '+name)
    fmax=max(math.sqrt(sum(v*v for v,m in zip(row[3:],mask) if m))
             for row,mask in zip(fs,geom['masks']))
    valid=normal and scf and F is not None and E0 is not None and dipflag and idflag and spinflag
    result={'name':name,'normal_end':normal,'last_scf_ediff':scf,'LDIPOL':dipflag,'IDIPOL_3':idflag,
            'ISPIN_1':spinflag,'F_eV':F,'E0_eV':E0,'max_free_force_eV_A':fmax,
            'electronic_verified':bool(valid),'force_below_0p02':fmax<.02,
            'geometry_ready_for_followup':bool(valid and fmax<.02)}
    return result

def run_case(idx):
    need(0<=idx<4,'Unexpected array index')
    need(os.environ.get('SLURM_NTASKS')=='64','Not in intended 64-task allocation')
    name=NAMES[idx];d=ROOT/'models'/name
    p=json.loads((ROOT/'production_manifest.json').read_text())
    for f in ('INCAR','POSCAR','KPOINTS','POTCAR'):
        rel='models/'+name+'/'+f
        need(sha(d/f)==p['input_sha256'][rel],'Production input changed '+rel)
    (d/'.execution_started').mkdir()
    with (d/'vasp.stdout').open('w') as out,(d/'vasp.stderr').open('w') as err:
        z=subprocess.run(['mpirun','-np','64','vasp_std'],cwd=d,stdout=out,stderr=err)
    need(z.returncode==0,'VASP exited nonzero '+name)
    record=read_results(name);save(d/'validation.json',record)
    print('VALIDATION|'+json.dumps(record,ensure_ascii=False),flush=True)
    need(record['electronic_verified'],'Static electronic verification failed '+name)
    (d/'STATIC_VERIFIED').touch()
    if record['force_below_0p02']:(d/'GEOMETRY_READY').touch()
    else:(d/'NEEDS_SHORT_REOPTIMIZATION').touch()

def status():
    print('CHECKED_AT|'+now())
    if not (ROOT/'submission.json').is_file():print('NOT_SUBMITTED');return
    rec=json.loads((ROOT/'submission.json').read_text())
    print('SUBMISSION|'+json.dumps(rec))
    jid=rec['job_id']
    for cmd in (['squeue','-r','-h','-j',jid,'-o','%i|%j|%T|%M|%R'],
                ['sacct','-X','-n','-P','-j',jid,'-o','JobID,State,ExitCode,Elapsed,Start,End']):
        z=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
        print('SLURM|'+z.stdout.strip()+'|RC='+str(z.returncode))
    old=PLAN['models'];new={}
    for n in NAMES:
        p=ROOT/'models'/n/'validation.json'
        if p.is_file():
            v=json.loads(p.read_text());new[n]=v
            print('MODEL|'+n+'|'+json.dumps(v))
        else:print('MODEL|'+n+'|NO_FINAL_VALIDATION')
    if len(new)==4 and all(r['electronic_verified'] for r in new.values()):
        def dd(arr,key):
            return (arr['M1_OH_H2O'][key]-arr['M1_OH'][key])-(arr['M0_OH_H2O'][key]-arr['M0_OH'][key])
        result={'old_E0_double_diff_eV':dd({n:old[n] for n in NAMES},'old_E0'),
                'new_E0_double_diff_eV':dd(new,'E0_eV'),
                'old_F_double_diff_eV':dd({n:old[n] for n in NAMES},'old_F'),
                'new_F_double_diff_eV':dd(new,'F_eV'),
                'needs_short_relax':[n for n in NAMES if not new[n]['force_below_0p02']]}
        save(ROOT/'P0_REPORT.json',result)
        print('P0_REPORT|'+json.dumps(result))
        print('NEXT_ACTION|Review this calibrated result before starting P1/P2')

def main():
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--check-inputs',action='store_true')
    g.add_argument('--submit',action='store_true')
    g.add_argument('--status',action='store_true')
    g.add_argument('--run-index',type=int)
    a=p.parse_args()
    if a.check_inputs:
        inp=build_inputs()
        for name,v in inp.items():
            print('INPUT_VERIFIED|'+name+'|'+json.dumps({'composition':v['composition'],
                 'POSCAR_sha256':hashlib.sha256(v['POSCAR'].encode()).hexdigest(),
                 'DIPOL':PLAN['models'][name]['DIPOL_fractional_z']}))
    elif a.submit:prepare_and_submit()
    elif a.status:status()
    else:run_case(a.run_index)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print('STOP|'+str(exc),file=sys.stderr,flush=True)
        raise
