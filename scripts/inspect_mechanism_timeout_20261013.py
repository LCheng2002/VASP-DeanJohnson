#!/usr/bin/env python3
"""Read-only structured diagnostic for the fixed 63585911 VASP batch; no job submission."""
import json,os,re,subprocess,datetime
from pathlib import Path
ROOT=Path.home()/'MgNi_LaH3_extensions_20261012'
NAMES=['L_parent','L_MgH2_reference','Mg_parent','X025_A','X025_B','H2_reference']
def tail(p,cap=1500000):
    if not p.is_file():return ''
    with p.open('rb') as f:
        f.seek(max(0,p.stat().st_size-cap));b=f.read()
    return b.decode(errors='replace')
def parse_outcar(p):
    s=tail(p,1750000)
    if not s:return {'exists':False}
    en=re.findall(r'free\s+energy\s+TOTEN\s*=\s*([-+\d.]+)',s)
    e0=re.findall(r'energy\(sigma->0\)\s*=\s*([-+\d.]+)',s)
    ionic=bool('reached required accuracy' in s)
    scf='aborting loop because EDIFF is reached' in s
    normal='General timing and accounting informations for this job' in s
    raw=list(re.finditer(r'TOTAL-FORCE \(eV/Angst\)',s))
    max_f=None;n_force=0
    if raw:
        frame=s[raw[-1].end():].splitlines()[2:]
        forces=[]
        for line in frame:
            row=line.split()
            if len(row)==6:
                try:forces.append([float(z) for z in row[-3:]])
                except ValueError:pass
            elif forces:break
        n_force=len(forces)
        if forces:max_f=max(sum(v*v for v in xyz)**.5 for xyz in forces)
    return {'exists':True,'bytes':p.stat().st_size,'energy_tail':list(map(float,en[-6:])), 'E0_tail':list(map(float,e0[-3:])), 'ionic_converged_in_tail':ionic,'last_scf_ediff_in_tail':scf,'normal_end':normal,'last_force_count':n_force,'last_force_max_all_eV_A':max_f,'tail_has_error':bool(re.search(r'ERROR|internal error|out of memory|ZBRENT|BRMIX|segmentation fault',s,re.I))}
def peek_validation(p):
    try:
        v=json.loads(p.read_text())
        return {k:v.get(k) for k in ('F_eV','E0_eV','verified','normal_end','ionic_converged','max_free_force_eV_A','ionic_steps','min_H_H_A')}
    except FileNotFoundError:return None
def snapshot():
    report={'at':datetime.datetime.now().astimezone().isoformat(),'jobs':{},'models':{},'submission':None}
    f=ROOT/'submission.json'
    if f.is_file():report['submission']=json.loads(f.read_text())
    ids='63585911,63585912,63585913,63585914,63585915'
    for label,cmd in [('sacct',['sacct','-n','-P','-j',ids,'-o','JobID,JobName%28,State,ExitCode,Elapsed,Timelimit,Start,End','-X']),('squeue',['squeue','-r','-h','-j',ids,'-o','%i|%j|%T|%M|%R'])]:
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=25)
        report['jobs'][label]={'rc':p.returncode,'stdout':p.stdout[-18000:],'stderr':p.stderr[-600:]}
    for name in NAMES:
        d=ROOT/name;v={'model':name,'exists':d.is_dir(),'complete_marker':(d/'VERIFIED_FINISHED').exists()}
        if not d.is_dir():report['models'][name]=v;continue
        for stage in ('01_relax','02_static'):
            p=d/stage
            files={}
            for fn in ('INCAR','KPOINTS','POSCAR','CONTCAR','OSZICAR','OUTCAR','vasp.stdout','vasp.stderr','validation.json'):
                f=p/fn
                if f.is_file():files[fn]={'bytes':f.stat().st_size,'updated':datetime.datetime.fromtimestamp(f.stat().st_mtime,datetime.timezone.utc).isoformat()}
            o={'files':files,'validation':peek_validation(p/'validation.json'),'outcar':parse_outcar(p/'OUTCAR')}
            if (p/'INCAR').is_file():
                incar=(p/'INCAR').read_text()
                o['INCAR']={key:re.findall('(?mi)^\\s*'+key+r'\\s*=\\s*([^!#\\n]+)',incar)[-1].strip() if re.search('(?mi)^\\s*'+key+r'\\s*=',incar) else None for key in ('NSW','EDIFF','EDIFFG','NELM','POTIM','ISPIN')}
            oz=tail(p/'OSZICAR',18000);ol=oz.splitlines()
            o['oszicar_tail']=ol[-9:]
            o['ionic_summary']=[s for s in ol if re.match(r'^\s*\d+\s+F=',s)][-8:]
            for fn in ('vasp.stdout','vasp.stderr'):
                s=tail(p/fn,18000).splitlines()
                o[fn+'_tail']=s[-14:]
            v[stage]=o
        report['models'][name]=v
    logs=ROOT/'logs'
    if logs.is_dir():
        report['slurm_errors']={}
        for p in sorted(logs.glob('*.err')):
            if '63585911' in p.name or '6358591' in p.name:report['slurm_errors'][p.name]=tail(p,1800)
    return report
if __name__=='__main__':
    result=snapshot()
    print('AT|'+result['at'])
    print('JOBS|'+json.dumps(result['jobs'],ensure_ascii=False))
    for name,v in result['models'].items():print('MODEL_AUDIT|'+name+'|'+json.dumps(v,ensure_ascii=False))
    print('SLURM_ERRORS|'+json.dumps(result.get('slurm_errors',{}),ensure_ascii=False))
