#!/usr/bin/env python3
"""Build 3 fourth-H starting configurations on a pinned H3 CONTCAR.
Run on Siyuan: python3 prepare_submit_H4.py --prepare / --submit / --status.
No POTCAR or authentication material leaves the HPC account.
"""
import argparse, datetime, hashlib, itertools, json, math, os, re, shutil, subprocess, sys
from pathlib import Path

PARENT_REL = 'TaskA_Mg17La2_H_enrichment_20261010/H3_addLa'
ROOT_REL = 'TaskA_Mg17La2_H4_preference_20261011'
PARENT_SHA = '7fdf86b6153683dda3ec22231ddcb039db5a48fb6bda0caf5a63308ce66e3b24'
POTCAR_SHA = '0d78810830835216364407b2112c9faf7d684083398095c1c9519fb18ca9de44'
MODELS = ('H4_La_partner', 'H4_Mg_top', 'H4_Mg_bridge')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(65536), b''): h.update(block)
    return h.hexdigest()
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def cart(v,c): return [sum(v[k]*c[k][j] for k in range(3)) for j in range(3)]
def frac(v,c):
    d=dot(c[0],cross(c[1],c[2])); assert abs(d)>1e-8
    return [dot(v,cross(c[1],c[2]))/d,dot(v,cross(c[2],c[0]))/d,dot(v,cross(c[0],c[1]))/d]
def norm(v): return math.sqrt(dot(v,v))
def mic(delta,c):
    f=frac(delta,c);n=[round(x) for x in f]; candidates=[]
    for t in itertools.product((-1,0,1),repeat=3):
        trans=cart([a+b for a,b in zip(n,t)],c)
        candidates.append([x-y for x,y in zip(delta,trans)])
    return min(candidates,key=lambda v:dot(v,v))
def distance(a,b,c): return norm(mic([x-y for x,y in zip(a,b)],c))
def parse_structure(content):
    lines=content.splitlines(); assert len(lines)>=8
    scale=float(lines[1]); assert scale>0
    cell=[[float(v)*scale for v in l.split()[:3]] for l in lines[2:5]]
    elements=lines[5].split();counts=list(map(int,lines[6].split()));n=sum(counts)
    i=7;select=lines[i].strip().lower().startswith('s')
    if select:i+=1
    mode=lines[i].strip().lower();assert mode.startswith(('d','c','k'));i+=1
    rows=lines[i:i+n];assert len(rows)==n
    coords=[];masks=[]
    for row in rows:
        parts=row.split();v=list(map(float,parts[:3]));assert len(v)==3 and all(math.isfinite(x) for x in v)
        coords.append(cart(v,cell) if mode.startswith('d') else [x*scale for x in v])
        m=[x.upper().startswith('T') for x in parts[3:6]] if select else [True]*3
        assert len(m)==3;masks.append(m)
    return {'lines':lines,'cell':cell,'elements':elements,'counts':counts,'xyz':coords,'masks':masks,'start':i,'select':select,'direct':mode.startswith('d'),'scale':scale,'species':sum(([e]*n for e,n in zip(elements,counts)),[])}
def config(p):
    result={}
    for line in p.read_text().splitlines():
        for seg in re.split(';',re.split('[!#]',line)[0]):
            if '=' in seg:
                k,v=seg.split('=',1);result[k.strip().upper()]=v.strip()
    return result

def outcar_summary(p,s):
    state={'last_scf_ediff':False,'ionic_converged':False,'normal_end':False};ionic_index=None
    forces=[];buf=[];reading=False
    with p.open(errors='replace') as f:
        for line in f:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and int(m.group(1))!=ionic_index:
                ionic_index=int(m.group(1));state['last_scf_ediff']=False
            if 'aborting loop because EDIFF is reached' in line:state['last_scf_ediff']=True
            if 'reached required accuracy' in line:state['ionic_converged']=True
            if 'General timing and accounting informations for this job' in line:state['normal_end']=True
            if 'free  energy   TOTEN' in line:state['F_eV']=float(line.split('=')[1].split()[0])
            if 'energy(sigma->0)' in line:state['E0_eV']=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buf=[];continue
            if reading:
                parts=line.split()
                if len(parts)==6:
                    try:buf.append(list(map(float,parts[3:6])));continue
                    except ValueError:pass
                if buf:forces=buf;reading=False
    assert len(forces)==len(s['xyz']), 'Missing/incomplete final force block'
    state['max_free_force_eV_A']=max(norm([x if move else 0 for x,move in zip(vec,mask)]) for vec,mask in zip(forces,s['masks']))
    return state

def candidate_summary(p,s):
    rows=sorted((distance(p,q,s['cell']),e,i+1) for i,(q,e) in enumerate(zip(s['xyz'],s['species'])))
    distances={e:sorted(r[0] for r in rows if r[1]==e) for e in s['elements']}
    return {'new_H_cart_A':p,'nearest8':rows[:8],'min_H_H_A':distances['H'][0], 'min_H_Mg_A':distances['Mg'][0], 'second_H_Mg_A':distances['Mg'][1], 'min_H_La_A':distances['La'][0]}

def build_models(content):
    s=parse_structure(content)
    assert s['elements']==['Mg','La','H'] and s['counts']==[68,8,3]
    assert all(all(v) for v in s['masks']), 'Expected all-free H3 parent; do not change constraints'
    xyz=s['xyz'];cell=s['cell']
    assert abs(cell[2][0])+abs(cell[2][1])<1e-8 and abs(cell[0][1])+abs(cell[0][2])+abs(cell[1][2])<1e-8
    # Atom IDs are 1-based from the audited H3 parent. Not a claim of exact symmetry.
    la=xyz[74];h3=xyz[78];delta=mic([x-y for x,y in zip(h3,la)],cell)
    la_partner=[la[0]-delta[0],la[1]+delta[1],la[2]+delta[2]]
    mg_top=[xyz[60][0],xyz[60][1],xyz[60][2]+1.80]
    span=mic([x-y for x,y in zip(xyz[55],xyz[53])],cell)
    mg_bridge=[xyz[53][j]+span[j]/2+(1.20 if j==2 else 0.0) for j in range(3)]
    positions=(la_partner,mg_top,mg_bridge)
    results={}
    for name,raw in zip(MODELS,positions):
        fr=frac(raw,cell);fr[0]%=1;fr[1]%=1;p=cart(fr,cell)
        info=candidate_summary(p,s)
        assert info['min_H_H_A']>=2.0, name+': too close to existing H'
        assert 1.60<=info['min_H_Mg_A']<=2.45, name+': Mg-H geometry outside screening bounds'
        assert p[2]>max(x[2] for x in xyz[:76])-3.0, name+': not an upper-surface-near site'
        assert cell[2][2]-p[2]>5.0, name+': insufficient upper-cell margin'
        if name=='H4_La_partner':
            assert 2.00<=info['min_H_La_A']<=2.80 and info['second_H_Mg_A']<=2.70
        else:
            assert info['min_H_La_A']>=3.10, name+': not separated from La-associated shell'
            if name=='H4_Mg_bridge': assert info['second_H_Mg_A']<=2.55
        lines=s['lines'][:s['start']+79]
        lines[0]=name+' from verified H3_addLa; H80 added, original 79 input atoms unchanged'
        lines[6]='68 8 4'
        coord=fr if s['direct'] else [x/s['scale'] for x in p]
        row='  '+'  '.join(f'{x:.16f}' for x in coord)+('   T   T   T' if s['select'] else '')
        output='\n'.join(lines+[row])+'\n'
        child=parse_structure(output)
        assert child['counts']==[68,8,4] and len(child['xyz'])==80
        assert max(abs(a-b) for r1,r2 in zip(child['xyz'][:79],xyz) for a,b in zip(r1,r2))<1e-12
        assert child['cell']==cell
        info.update({'model':name,'new_H_atom_1based':80,'composition':'Mg68La8H4','initial_coordinates_sha256':hashlib.sha256(output.encode()).hexdigest(),'screening_only_not_energy_ranking':True})
        results[name]={'POSCAR':output,'metadata':info}
    return results

JOB_SCRIPT='''#!/bin/bash
#SBATCH --job-name=Mg17La2-H4
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive
#SBATCH --array=0-2%2
#SBATCH --time=12:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
models=(H4_La_partner H4_Mg_top H4_Mg_bridge)
model="${models[$SLURM_ARRAY_TASK_ID]}"
cd "$model"
mkdir .execution_started
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
(cd 01_relax && mpirun -np "$SLURM_NTASKS" vasp_std && python3 ../../prepare_submit_H4.py --verify relax)
cp 01_relax/CONTCAR 02_static/POSCAR
(cd 02_static && mpirun -np "$SLURM_NTASKS" vasp_std && python3 ../../prepare_submit_H4.py --verify static)
touch VERIFIED_FINISHED
'''

def prepare(root,parent):
    p=parent/'01_relax/CONTCAR';assert sha(p)==PARENT_SHA, 'H3 source CONTCAR changed'
    c=p.read_text();s=parse_structure(c)
    a=outcar_summary(parent/'01_relax/OUTCAR',parse_structure((parent/'01_relax/POSCAR').read_text()))
    b=outcar_summary(parent/'02_static/OUTCAR',parse_structure((parent/'02_static/POSCAR').read_text()))
    assert a['ionic_converged'] and a['last_scf_ediff'] and a['normal_end'] and a['max_free_force_eV_A']<0.02
    assert b['last_scf_ediff'] and b['normal_end'] and b['max_free_force_eV_A']<0.02
    assert abs(b['F_eV']-(-146.97883350))<1e-7
    assert sha(parent/'02_static/POSCAR')==PARENT_SHA
    assert sha(parent/'01_relax/POTCAR')==POTCAR_SHA and sha(parent/'02_static/POTCAR')==POTCAR_SHA
    relcfg=config(parent/'01_relax/INCAR');stcfg=config(parent/'02_static/INCAR')
    assert relcfg['ISPIN']=='1' and relcfg['ISMEAR']=='1' and float(relcfg['SIGMA'])==0.1
    assert relcfg['IBRION']=='2' and relcfg['ISIF']=='2' and relcfg['NSW']=='300' and float(relcfg['EDIFFG'])==-0.02
    assert 'NELECT' not in relcfg and 'NELECT' not in stcfg, 'Fixed electron count unsuitable when adding H'
    assert sha(parent/'01_relax/KPOINTS')==sha(parent/'02_static/KPOINTS')
    new=build_models(c)
    if root.exists():
        manifest=json.loads((root/'manifest.json').read_text())
        assert manifest['parent_contcar_sha256']==PARENT_SHA
        assert set(manifest['models'])==set(MODELS)
        for name in MODELS:
            for fn,h in manifest['input_hashes'][name].items():assert sha(root/name/fn)==h,'Prepared input changed: '+name+'/'+fn
        print('ALREADY_PREPARED|'+str(root));return manifest
    root.mkdir();(root/'logs').mkdir()
    shutil.copy2(Path(__file__).resolve(),root/'prepare_submit_H4.py')
    (root/'job_array.slurm').write_text(JOB_SCRIPT)
    (root/'parent_H3_CONTCAR').write_text(c)
    manifest={'prepared_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'parent_model':PARENT_REL,'parent_contcar_sha256':PARENT_SHA,'parent_static':b,'parent_potcar_sha256':POTCAR_SHA,'models':{},'input_hashes':{},'job_type':'three 80-atom Relax->Static trials; array 0-2%2; no NEB','method_note':'Same H3 parent; only fourth-H starting position varies; NOT exhaustive global configuration search.'}
    for name,data in new.items():
        d=root/name;d.mkdir()
        for stage in ('01_relax','02_static'):
            (d/stage).mkdir()
            for fn in ('INCAR','KPOINTS','POTCAR'):shutil.copy2(parent/stage/fn,d/stage/fn)
        (d/'01_relax/POSCAR').write_text(data['POSCAR'])
        (d/'model_metadata.json').write_text(json.dumps(data['metadata'],ensure_ascii=False,indent=2))
        manifest['models'][name]=data['metadata'];manifest['input_hashes'][name]={}
        for stage in ('01_relax','02_static'):
            for fn in ('INCAR','KPOINTS','POTCAR'):
                manifest['input_hashes'][name][stage+'/'+fn]=sha(d/stage/fn)
        manifest['input_hashes'][name]['01_relax/POSCAR']=sha(d/'01_relax/POSCAR')
    manifest['job_script_sha256']=sha(root/'job_array.slurm')
    manifest['validator_script_sha256']=sha(root/'prepare_submit_H4.py')
    subprocess.run(['bash','-n',str(root/'job_array.slurm')],check=True)
    (root/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'PREPARED':str(root),'models':manifest['models']},ensure_ascii=False,indent=2))
    return manifest

def verify(stage):
    p=Path.cwd();s=parse_structure((p/'POSCAR').read_text());assert s['counts']==[68,8,4]
    r=outcar_summary(p/'OUTCAR',s);r['stage']=stage
    valid=r['last_scf_ediff'] and r['normal_end'] and r['max_free_force_eV_A']<0.0200001 and (r['ionic_converged'] or stage=='static')
    if (p/'CONTCAR').is_file():
        end=parse_structure((p/'CONTCAR').read_text());assert end['counts']==s['counts']
        r['final_H_neighbours']=[]
        for i in range(76,80):
            neigh=sorted((distance(end['xyz'][i],v,end['cell']),end['species'][j],j+1) for j,v in enumerate(end['xyz']) if j!=i)
            r['final_H_neighbours'].append({'H_atom_1based':i+1,'nearest6':neigh[:6]})
        hh=[distance(end['xyz'][i],end['xyz'][j],end['cell']) for i,j in itertools.combinations(range(76,80),2)]
        r['minimum_H_H_A']=min(hh);r['possible_H2_geometry']=min(hh)<1.0
    r['numerically_verified']=bool(valid)
    (p/'validation.json').write_text(json.dumps(r,ensure_ascii=False,indent=2))
    print(json.dumps(r,ensure_ascii=False));return 0 if valid else 40

def main():
    a=argparse.ArgumentParser(description=__doc__);g=a.add_mutually_exclusive_group(required=True)
    g.add_argument('--prepare',action='store_true');g.add_argument('--submit',action='store_true');g.add_argument('--status',action='store_true');g.add_argument('--verify',choices=['relax','static'])
    args=a.parse_args()
    if args.verify:return verify(args.verify)
    root=Path.home()/ROOT_REL;parent=Path.home()/PARENT_REL
    if args.status:
        p=root/'submission.json'
        if not p.is_file():print('NO_CONFIRMED_JOB_ID');return 0
        rec=json.loads(p.read_text());print(json.dumps(rec,indent=2));jid=rec['job_id'].split(';')[0]
        subprocess.run(['squeue','-r','-j',jid,'-o','%.20i %.20j %.12T %.10M %R'],timeout=30)
        subprocess.run(['sacct','-n','-P','-j',jid,'-o','JobID,JobName%24,State,ExitCode,Elapsed,Start,End'],timeout=30)
        for name in MODELS:
            p=root/name/'02_static/validation.json'
            if p.is_file():print(name,p.read_text())
        return 0
    manifest=prepare(root,parent)
    if not args.submit:return 0
    receipt=root/'submission.json'
    if receipt.is_file():print('ALREADY_SUBMITTED|'+receipt.read_text());return 0
    assert sha(root/'job_array.slurm')==manifest['job_script_sha256']
    assert sha(root/'prepare_submit_H4.py')==manifest['validator_script_sha256']
    # An ambiguous previous attempt keeps this lock. Never delete it to retry blindly.
    marker=root/'.submission_attempted';fd=os.open(str(marker),os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600);os.close(fd)
    res=subprocess.run(['sbatch','--parsable','job_array.slurm'],cwd=root,capture_output=True,text=True,timeout=45)
    if res.returncode:
        (root/'submission_error.txt').write_text(res.stderr);raise RuntimeError('sbatch failed; inspect submission_error.txt before retry')
    jid=res.stdout.strip();assert re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?',jid),'Ambiguous sbatch output; do not retry'
    rec={'job_id':jid,'array_indices':dict(enumerate(MODELS)),'max_concurrent_models':2,'tasks_per_model':64,'walltime_per_model':'12:00:00','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'directory':str(root)}
    receipt.write_text(json.dumps(rec,ensure_ascii=False,indent=2));print('SUBMITTED_ARRAY_JOB_ID|'+jid);print(json.dumps(rec,ensure_ascii=False))
    return 0

if __name__=='__main__':
    try:sys.exit(main())
    except (AssertionError,ValueError,OSError,RuntimeError,subprocess.SubprocessError) as exc:
        print('STOPPED|'+str(exc),file=sys.stderr);sys.exit(1)
