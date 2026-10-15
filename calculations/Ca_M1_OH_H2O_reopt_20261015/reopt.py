#!/usr/bin/env python3
"""One and only one M1_OH_H2O ionic re-relaxation + matched Static.

Starts from completed/canonical P0 job 63686231_3. All licensed PAW data remain
inside the user's existing HPC account. Do not touch old models or jobs.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())
HOME=Path.home()
SOURCE=HOME/PLAN['source']['production_directory']
ROOT=HOME/PLAN['batch']
CASE=PLAN['model']
FLOAT=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'

def need(ok,msg):
    if not ok:raise RuntimeError('PRECHECK_OR_VALIDATION_FAILED|'+msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now().astimezone().isoformat()
def save(p,obj):
    p=Path(p);tmp=p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n');tmp.replace(p)
def load_parent_parser():
    path=SOURCE/'calibrate.py'
    need(path.is_file(),'Original P0 parser unavailable')
    spec=importlib.util.spec_from_file_location('verified_ca_p0_parser',path)
    need(spec is not None and spec.loader is not None,'Cannot load P0 parser')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m
def incar_map(raw):
    d={}
    for line in raw.splitlines():
        line=re.split(r'[!#]',line)[0]
        for z in line.split(';'):
            if '=' not in z:continue
            k,v=z.split('=',1);k=k.strip().upper();v=v.strip()
            need(not(k in d and d[k]!=v),'Conflicting source INCAR '+k)
            d[k]=v
    return d
def patch_incar(raw,mode):
    changes={
      'SYSTEM':'Ca_M1_OH_H2O_20261015_'+mode,
      'ISIF':'2', 'IBRION':'2' if mode=='relax' else '-1',
      'NSW':'160' if mode=='relax' else '0',
      'POTIM':'0.10', 'EDIFFG':'-0.012',
      'LREAL':'.FALSE.', 'EDIFF':'1E-7',
      'ISTART':'0','ICHARG':'2',
      'LWAVE':'.FALSE.','LCHARG':'.FALSE.','LVTOT':'.FALSE.'
    }
    original=incar_map(raw)
    need(original.get('LDIPOL','').upper() in ('T','.TRUE.') and original.get('IDIPOL')=='3','Calibrated dipole correction not active')
    need(original.get('DIPOL') is not None,'No calibrated fractional DIPOL')
    for key,value in {'ENCUT':'520','ISPIN':'1','ISMEAR':'1','SIGMA':'0.2','NBANDS':'160','ISYM':'0'}.items():
        need(original.get(key)==value,'Changed canonical numerical setting '+key)
    lines=[l for l in raw.splitlines() if re.split(r'[!#]',l)[0].split('=',1)[0].strip().upper() not in changes]
    lines.extend(k+' = '+v for k,v in changes.items())
    updated='\n'.join(lines)+'\n'
    vals=incar_map(updated)
    need(len(vals)==len(set(vals)),'Duplicate resulting INCAR keys')
    return updated
def validate_parent():
    need(SOURCE.is_dir(),'Canonical P0 directory missing')
    submit=SOURCE/'submission.json'
    need(submit.is_file(),'P0 submission record missing')
    receipt=json.loads(submit.read_text())
    need(receipt.get('job_id')=='63686231' and receipt.get('four_statics_only') and CASE in receipt['models'],'Not expected canonical completed P0')
    original=SOURCE/CASE
    vfile=original/'validation.json';marker=original/'CALIBRATION_COMPLETE'
    need(marker.is_file() and vfile.is_file(),'M1_OH_H2O P0 was not completed')
    val=json.loads(vfile.read_text())
    need(val.get('electronic_calibration_verified') and val.get('last_scf_ediff') and val.get('normal_end'),'P0 electronic gate failed')
    need(not val.get('geometry_force_gate_passed') and val.get('needs_consistent_relaxation'),'P0 model already satisfied geometry; refuse duplicate work')
    ref=PLAN['source'];need(abs(val['F_eV']-ref['parent_calibration_F_eV'])<1e-7,'Unexpected P0 F energy')
    need(abs(val['E0_eV']-ref['parent_calibration_E0_eV'])<1e-7,'Unexpected P0 E0 energy')
    need(abs(val['max_free_force_eV_A']-ref['parent_calibration_max_free_force_eV_A'])<1e-7,'Unexpected P0 force')
    ev=val['effective_parameters'];settings=PLAN['identical_numerical_settings']
    for key,num in [('ENCUT',520),('ISPIN',1),('ISMEAR',1),('SIGMA',.2),('NBANDS',160),('NELECT',215),('IDIPOL',3)]:
        need(ev.get(key) is not None and abs(ev[key]-num)<1e-5,'Unexpected actual P0 VASP '+key)
    need(ev.get('LDIPOL')=='T','P0 LDIPOL disabled')
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    for file in ('POSCAR','INCAR','KPOINTS','POTCAR','model.json'):
        key=CASE+'/'+file
        need(manifest['input_hashes'].get(key)==sha(original/file),'P0 input modified '+file)
    need(sha(original/'POTCAR')==settings['original_exact_POTCAR_SHA256'],'P0 licensed POTCAR hash differs')
    incar=(original/'INCAR').read_text();params=incar_map(incar)
    need(params.get('DIPOL'),'No calibrated DIPOL')
    from_plan=[float(x) for x in params['DIPOL'].split()]
    need(len(from_plan)==3 and all(0<=x<=1 for x in from_plan),'Not a fractionally specified P0 DIPOL')
    m=load_parent_parser()
    poscar=(original/'POSCAR').read_text()
    geom=m.parse_structure(poscar)
    need(geom['composition']==PLAN['stoichiometry'] and len(geom['species'])==PLAN['expected_atoms'],'Changed P0 composition')
    need(sum(not any(x) for x in geom['masks'])==32,'P0 frozen layers altered')
    need((original/'KPOINTS').is_file() and (original/'POTCAR').stat().st_size>100,'Missing P0 authorized inputs')
    need((original/'OUTCAR').is_file(),'P0 OUTCAR unavailable')
    need(geom['species'][95]=='Ca' and geom['species'][96:]==['O','O','H','H','H'],
        'Unexpected M1_OH_H2O atom ID map: Ca96, O97/O98, H99–101')
    return original,val,geom
def query_duplicate():
    names={'Ca_M1_OH_H2O_reopt'}
    proc=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%j'],
                        capture_output=True,text=True,timeout=30)
    need(proc.returncode==0,'Cannot check existing Slurm jobs')
    need(not any(n.strip() in names for n in proc.stdout.splitlines()),'An M1_OH_H2O reopt Slurm job is already active')
def prepare():
    original,reference,geom=validate_parent()
    if ROOT.exists():
        if (ROOT/'submission.json').is_file():
            print('ALREADY_SUBMITTED|'+(ROOT/'submission.json').read_text(),flush=True);return False
        raise RuntimeError('Existing unconfirmed production root: inspect; never overwrite')
    query_duplicate()
    relaxed_incar=patch_incar((original/'INCAR').read_text(),'relax')
    static_incar=patch_incar((original/'INCAR').read_text(),'static')
    ROOT.mkdir(mode=0o700)
    (ROOT/'logs').mkdir()
    for name,text in [('plan.json',(HERE/'plan.json').read_text()),('reopt.py',Path(__file__).read_text())]:
        (ROOT/name).write_text(text)
    for mode,source_incar in [('01_relax',relaxed_incar),('02_static',static_incar)]:
        dest=ROOT/mode;dest.mkdir()
        (dest/'INCAR').write_text(source_incar)
        shutil.copyfile(original/'KPOINTS',dest/'KPOINTS')
        shutil.copyfile(original/'POTCAR',dest/'POTCAR')
        (dest/'POTCAR').chmod(0o600)
        if mode=='01_relax':
            shutil.copyfile(original/'POSCAR',dest/'POSCAR')
    job='''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
python3 reopt.py --run
'''
    (ROOT/'job.slurm').write_text(job)
    p=subprocess.run(['bash','-n',str(ROOT/'job.slurm')],capture_output=True,text=True)
    need(p.returncode==0,'Invalid Slurm job script')
    checks={str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file()}
    save(ROOT/'manifest.json',{'source_job_id':ref_job(),'source_validation_F_eV':reference['F_eV'],
        'initial_P0_force_eV_A':reference['max_free_force_eV_A'],'models':[CASE],
        'file_hashes':checks,'authorized_scope':'ONLY_M1_OH_H2O_REOPT_STATIC',
        'input_unchanged_except_geometry_relax_and_output_flags':True,'created_at':now()})
    print('PREPARED|'+str(ROOT)+'|one M1_OH_H2O Relax+Static|no old task restarted',flush=True)
    return True
def ref_job():return PLAN['source']['parent_job_id']
def submit():
    need(sys.platform.startswith('linux'),'Submit only on HPC login host')
    need(shutil.which('sbatch'),'No Slurm environment')
    if (ROOT/'submission.json').is_file():
        print('EXISTING_SUBMISSION|'+(ROOT/'submission.json').read_text(),flush=True);return
    was=prepare()
    need(was,'Existing production root; stop')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for key,value in manifest['file_hashes'].items():
        need(sha(ROOT/key)==value,'Production file hash mismatch '+key)
    (ROOT/'.submission_attempted').mkdir(exist_ok=False)
    p=subprocess.run(['sbatch','--parsable','--job-name=Ca_M1_OH_H2O_reopt','job.slurm'],
                     cwd=ROOT,text=True,capture_output=True,timeout=45)
    need(p.returncode==0,'sbatch failed '+p.stderr[-350:])
    jid=p.stdout.strip().split(';')[0]
    need(re.fullmatch(r'[0-9]+',jid),'Noncanonical sbatch receipt, inspect lock')
    save(ROOT/'submission.json',{'job_id':jid,'name':CASE,'submitted_at':now(),
        'reference_job':ref_job(),'reopt_then_static':True,
        'other_Ca_jobs_submitted':False,'cancelled_CaFresh_resubmitted':False})
    print('SUBMITTED|Ca_M1_OH_H2O_reopt|'+jid+'|one Relax->Static task',flush=True)

def read_outcar(folder,expected,geometry):
    maxfree=F=E0=None
    end=False;ediff=False;ionic=False;step=None
    forces=[];rows=[];reading=False;counter=0
    final_params={}
    for line in (folder/'OUTCAR').open(errors='replace'):
        match=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
        if match and match.group(1)!=step:step=match.group(1);ediff=False
        if 'aborting loop because EDIFF is reached' in line:ediff=True
        if 'reached required accuracy' in line:ionic=True
        if 'General timing and accounting informations' in line:end=True
        if 'free  energy   TOTEN' in line:
            F=float(line.split('=')[1].split()[0]);counter+=1
        if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
        for key in ('ISPIN','NBANDS','NELECT','IDIPOL','ENCUT','ISMEAR','SIGMA'):
            m=re.search(r'\b'+key+r'\s*=\s*('+FLOAT+')',line)
            if m:final_params[key]=float(m.group(1).replace('D','E'))
        m=re.search(r'\bLDIPOL\s*=\s*([TF])',line)
        if m:final_params['LDIPOL']=m.group(1)
        if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;rows=[];continue
        if reading:
            vals=line.split()
            if len(vals)==6:
                try:rows.append(list(map(float,vals)));continue
                except ValueError:pass
            if rows:forces=rows;reading=False
    need(len(forces)==len(expected['species']),'Final force vector count incomplete')
    maxfree=max(math.sqrt(sum(value*value for value,free in zip(row[3:],mask) if free))
                for row,mask in zip(forces,expected['masks']))
    need(math.isfinite(maxfree) and F is not None and E0 is not None,'Missing VASP final energies or forces')
    for key,value in [('ISPIN',1),('NBANDS',160),('NELECT',215),('IDIPOL',3),('ENCUT',520),('ISMEAR',1),('SIGMA',.2)]:
        need(key in final_params and abs(final_params[key]-value)<1e-4,'Changed effective P0 physical settings '+key)
    need(final_params.get('LDIPOL')=='T','LDIPOL not active')
    val={'phase':folder.name,'F_eV':F,'E0_eV':E0,'max_free_force_eV_A':maxfree,
         'ionic_converged':ionic,'last_scf_ediff':ediff,'normal_end':end,
         'full_force_table_atoms':len(forces),'ionic_steps':counter,
         'effective_parameters':final_params}
    lim=PLAN['relaxation']['max_free_force_eV_A_lt'] if folder.name=='01_relax' else PLAN['static']['max_free_force_eV_A_lt']
    val['force_limit_eV_A']=lim
    val['verified']=bool(end and ediff and maxfree<lim and (folder.name=='02_static' or ionic))
    return val
def run():
    need(os.environ.get('SLURM_NTASKS')=='64','Must run with 64 Slurm tasks')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    rec=json.loads((ROOT/'submission.json').read_text())
    need(rec.get('name')==CASE and rec.get('reopt_then_static'),'Unexpected submission record')
    for rel,checksum in manifest['file_hashes'].items():
        need(sha(ROOT/rel)==checksum,'Pinned input changed '+rel)
    (ROOT/'.execution_started').mkdir(exist_ok=False)
    m=load_parent_parser()
    inp=m.parse_structure((ROOT/'01_relax/POSCAR').read_text())
    geom=inp
    for mode in ('01_relax','02_static'):
        d=ROOT/mode
        if mode=='02_static':
            shutil.copyfile(ROOT/'01_relax/CONTCAR',d/'POSCAR')
            inp_static=m.parse_structure((d/'POSCAR').read_text())
            need(inp_static['composition']==geom['composition'] and inp_static['masks']==geom['masks'],'Relax altered atom roster or frozen layer constraints')
            need(max(abs(a-b) for x,y in zip(inp_static['cell'],geom['cell']) for a,b in zip(x,y))<1e-6,'Relax changed fixed slab cell')
        with (d/'vasp.stdout').open('w') as out,(d/'vasp.stderr').open('w') as err:
            p=subprocess.run(['mpirun','-np','64','vasp_std'],cwd=d,stdout=out,stderr=err)
        need(p.returncode==0,'VASP returned nonzero in '+mode)
        info=read_outcar(d,m.parse_structure((d/'POSCAR').read_text()),geom)
        con=m.parse_structure((d/'CONTCAR').read_text())
        need(con['composition']==geom['composition'] and con['masks']==geom['masks'],'CONTCAR changed stoichiometry/constraints '+mode)
        need(max(abs(a-b) for x,y in zip(con['cell'],geom['cell']) for a,b in zip(x,y))<1e-6,'CONTCAR changed cell '+mode)
        info['final_O_H_distances_A']={
            'original_O97_to_H99':m.mic_length(con['xyz'][96],con['xyz'][98],con['cell']),
            'water_O98_to_H100':m.mic_length(con['xyz'][97],con['xyz'][99],con['cell']),
            'water_O98_to_H101':m.mic_length(con['xyz'][97],con['xyz'][100],con['cell']),
        }
        info['water_remains_molecular_initial_bond_pattern']=all(
            0.7<info['final_O_H_distances_A'][key]<1.25
            for key in ('water_O98_to_H100','water_O98_to_H101'))
        save(d/'validation.json',info)
        print('VALIDATION|'+mode+'|'+json.dumps(info,ensure_ascii=False),flush=True)
        need(info['verified'],'Relax/Static verification failed '+mode+'; no false COMPLETED claim')
        (d/'VERIFIED_FINISHED').touch()
    result=json.loads((ROOT/'02_static/validation.json').read_text())
    sourceval=json.loads((SOURCE/CASE/'validation.json').read_text())
    final={'new_F_eV':result['F_eV'],'new_E0_eV':result['E0_eV'],
           'delta_F_eV':result['F_eV']-sourceval['F_eV'],
           'delta_E0_eV':result['E0_eV']-sourceval['E0_eV'],
           'max_free_force_eV_A':result['max_free_force_eV_A'],
           'water_remains_molecular_initial_bond_pattern':result['water_remains_molecular_initial_bond_pattern'],
           'warning':'Cannot recompute final Ca-vs-Mg adsorption double difference without consistently relaxed M1_OH; only requested M1_OH_H2O was run.'}
    save(ROOT/'final_result.json',final)
    (ROOT/'VERIFIED_FINISHED').touch()
    print('REOPT_COMPLETE|'+json.dumps(final,ensure_ascii=False),flush=True)
def status():
    print('CHECKED_AT|'+now())
    if not (ROOT/'submission.json').is_file():print('NOT_SUBMITTED');return
    rec=json.loads((ROOT/'submission.json').read_text())
    print('SUBMISSION|'+json.dumps(rec))
    jid=rec['job_id']
    for args in (['squeue','-j',jid,'-h','-o','%i|%j|%T|%M|%R'],
                 ['sacct','-X','-n','-P','-j',jid,'-o','JobID,JobName%28,State,ExitCode,Elapsed,Start,End']):
        r=subprocess.run(args,text=True,capture_output=True,timeout=30)
        print('QUEUE_OR_ACCOUNTING|'+str(r.returncode)+'|'+r.stdout.strip())
    for mode in ('01_relax','02_static'):
        p=ROOT/mode/'validation.json'
        print('MODEL|'+mode+'|'+(p.read_text() if p.is_file() else 'NO_FINAL_VALIDATION'))
    p=ROOT/'final_result.json'
    if p.is_file():print('RESULT|'+p.read_text())
def main():
    parser=argparse.ArgumentParser()
    opt=parser.add_mutually_exclusive_group(required=True)
    opt.add_argument('--submit',action='store_true')
    opt.add_argument('--run',action='store_true')
    opt.add_argument('--status',action='store_true')
    opt.add_argument('--precheck',action='store_true')
    args=parser.parse_args()
    if args.precheck:
        orig,v,g=validate_parent()
        print('PARENT_VALID|'+CASE+'|P0F='+str(v['F_eV'])+'|P0_E0='+str(v['E0_eV'])+
              '|P0_max_force='+str(v['max_free_force_eV_A'])+'|source='+str(orig))
    elif args.submit:submit()
    elif args.run:run()
    else:status()
if __name__=='__main__':
    try:main()
    except Exception as error:
        print('STOP|'+str(error),file=sys.stderr,flush=True)
        raise
