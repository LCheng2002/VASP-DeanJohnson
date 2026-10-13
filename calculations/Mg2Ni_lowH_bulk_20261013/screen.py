#!/usr/bin/env python3
"""Five small bulk Mg2NiHx structures: prepare, Slurm submit, run, verify, status.
Only initial chemical models and thermodynamic endpoints; no NEB or 147-atom restarts.
Licensed POTCARs are read and assembled locally on the HPC; never exported.
"""
import argparse,datetime,hashlib,itertools,json,math,os,re,shutil,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())
HOME=Path.home()
ROOT=HOME/PLAN['batch']
NAMES=[v['name'] for v in PLAN['models']]
SHIFTS=list(itertools.product((-1,0,1),repeat=3))
def require(ok,msg):
    if not ok: raise RuntimeError('CHECK_FAILED|'+msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p=Path(p);t=p.with_name(p.name+'.tmp')
    t.write_text(json.dumps(obj,indent=2,ensure_ascii=False));t.replace(p)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def norm(a):return math.sqrt(dot(a,a))
def cart(f,c):return [sum(f[k]*c[k][j] for k in range(3)) for j in range(3)]
def mic_distance(f,g,c):
    dif=[f[i]-g[i] for i in range(3)]
    return min(norm(cart([dif[k]-s[k] for k in range(3)],c)) for s in SHIFTS)
def parse(s):
    lines=s.splitlines();require(len(lines)>8,'Invalid POSCAR')
    scale=float(lines[1]);require(scale>0,'POSCAR scale must be positive')
    cell=[[float(z)*scale for z in row.split()[:3]] for row in lines[2:5]]
    els=lines[5].split();nums=[int(x) for x in lines[6].split()]
    sp=sum(([e]*n for e,n in zip(els,nums)),[])
    require(len(els)==len(nums) and len(sp)>0,'Invalid species list')
    k=7;select=lines[k].strip().lower().startswith('s');k+=int(select)
    direct=lines[k].strip().lower().startswith('d');k+=1
    xyz=[]
    for line in lines[k:k+len(sp)]:
        q=line.split();v=[float(x) for x in q[:3]];require(len(v)==3,'Missing positions')
        if direct:xyz.append([x%1 for x in v])
        else:
            # Inverse cell via scalar triple product; normalize fractional coordinates.
            a,b,c=cell
            cross=lambda x,y:[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
            det=dot(a,cross(b,c));require(abs(det)>1e-6,'Singular cell')
            v=[z*scale for z in v]
            xyz.append([(dot(v,cross(b,c))/det)%1,(dot(v,cross(c,a))/det)%1,(dot(v,cross(a,b))/det)%1])
    require(len(xyz)==len(sp),'Incomplete positions')
    return {'cell':cell,'species':sp,'frac':xyz}
def poscar(s,name):
    els=list(dict.fromkeys(s['species']))
    require(sum(([e]*s['species'].count(e) for e in els),[])==s['species'],'Elements ungrouped')
    lines=[name,'1.0']+[' '.join('%.14f'%v for v in row) for row in s['cell']]
    lines+=[' '.join(els),' '.join(str(s['species'].count(e)) for e in els),'Direct']
    lines+=[' '.join('%.14f'%(x%1) for x in f) for f in s['frac']]
    return '\n'.join(lines)+'\n'
def nearest_metal(s,f):
    return sorted([(mic_distance(f,q,s['cell']),e,i+1) for i,(q,e) in enumerate(zip(s['frac'],s['species'])) if e!='H'])
def read_outcar(path):
    F=E0=None;normal=ionic=scf=False;it=None;forces=[];buffer=[];reading=False;steps=0
    with Path(path).open(errors='replace') as q:
        for line in q:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and it!=m.group(1):it=m.group(1);scf=False
            if 'aborting loop because EDIFF is reached' in line:scf=True
            if 'reached required accuracy' in line:ionic=True
            if 'General timing and accounting informations for this job' in line:normal=True
            if 'free  energy   TOTEN' in line:F=float(line.split('=')[1].split()[0]);steps+=1
            if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buffer=[];continue
            if reading:
                cols=line.split()
                if len(cols)==6:
                    try:buffer.append([float(z) for z in cols[3:6]]);continue
                    except ValueError:pass
                if buffer:forces=buffer;reading=False
    return {'F_eV':F,'E0_eV':E0,'normal_end':normal,'last_scf_ediff':scf,
            'ionic_converged':ionic,'ionic_steps':steps,
            'force_count':len(forces),'max_force_eV_A':max((norm(x) for x in forces),default=None)}
def split_potcar(p):
    chunks=[];buffer=[]
    for line in Path(p).read_bytes().splitlines(keepends=True):
        buffer.append(line)
        if line.strip()==b'End of Dataset':chunks.append(b''.join(buffer));buffer=[]
    require(len(chunks)>0 and not buffer,'Could not split licensed POTCAR')
    blocks={}
    for data in chunks:
        title=re.search(rb'TITEL\s*=\s*PAW_PBE\s+(\w+)',data)
        require(title is not None,'Only labeled PBE potentials are permitted')
        e=title.group(1).decode().split('_')[0]
        require(e not in blocks,'Duplicate potential block')
        blocks[e]=data
    return blocks
def model_structure(base,meta):
    require(base['species']==['Mg']*12+['Ni']*6,'Primitive reference Mg12Ni6 mismatch')
    cell=[list(v) for v in base['cell']];cell[2]=[2*x for x in cell[2]]
    species=[];coords=[]
    for element in ('Mg','Ni'):
        for rep in (0,1):
            for e,f in zip(base['species'],base['frac']):
                if e!=element:continue
                species.append(e);coords.append([f[0],f[1],(f[2]+rep)/2])
    require(species==['Mg']*24+['Ni']*12,'Supercell order mismatch')
    base_super={'cell':cell,'species':species,'frac':coords}
    info=[];pos=[]
    for label,rep in meta['h_sites']:
        require(rep in (0,1) and label in PLAN['fractional_sites_primitive_parent'],'Unknown void label')
        site=PLAN['fractional_sites_primitive_parent'][label]
        f=site['f'];q=[f[0],f[1],(f[2]+rep)/2]
        env=nearest_metal(base_super,q)
        actual=sorted(e for _,e,_ in env[:4]);expected=sorted(site['expected_four_nn'])
        require(actual==expected,'Incorrect tetrahedral local coordination for '+label+str(rep))
        require(1.55<env[0][0]<2.15 and env[3][0]<2.20,'Bad H-metal tetrahedron distance')
        require(env[4][0]>2.45,'Not a distinct metal tetrahedron')
        info.append({'site':label,'half_cell':rep,'fractional':q,'nearest_metal_4':env[:4],
                     'fifth_neighbour':env[4]})
        pos.append(q)
    hh=[mic_distance(a,b,cell) for a,b in itertools.combinations(pos,2)]
    require(not hh or min(hh)>3.0,'Initial H-H collision')
    require(len(pos)==len(meta['h_sites']),'H placement mismatch')
    s={'cell':cell,'species':species+['H']*len(pos),'frac':coords+pos}
    require(len(s['species'])==meta['atom_count'],'Model atom count mismatch')
    return s,{'positions':info,'min_initial_HH_A':min(hh) if hh else None,
              'nominal_x':len(pos)/12,'coordination_only_not_refined_X3_X7':True}
def validate_external():
    src=HOME/PLAN['source']['path_on_HPC']
    require(src.is_file() and sha(src)==PLAN['source']['sha256'],'Pinned relaxed Mg2Ni source POSCAR mismatch')
    base=parse(src.read_text())
    require(base['species']==['Mg']*12+['Ni']*6,'Wrong source chemical formula')
    static=read_outcar(src.parent/'OUTCAR')
    require(static['normal_end'] and static['last_scf_ediff'] and static['max_force_eV_A']<.020001,'Mg2Ni source static not validated')
    parent=src.parent.parent/'01_relax/OUTCAR'
    pcheck=read_outcar(parent)
    require(pcheck['ionic_converged'] and pcheck['normal_end'],'Mg2Ni source relaxation not converged')
    h2=HOME/PLAN['h2_reference']['root_on_HPC']
    v=json.loads((h2/'02_static/validation.json').read_text())
    require((h2/'VERIFIED_FINISHED').is_file() and v.get('verified') and v.get('normal_end') and v.get('last_scf_ediff'),'H2 reference invalid')
    require(abs(v['F_eV']-PLAN['h2_reference']['expected_F_eV'])<.001,'H2 reference energy changed')
    pot=HOME/PLAN['potential']['source_on_HPC']
    require(pot.is_file() and sha(pot)==PLAN['potential']['sha256'],'Pinned licensed Mg-Ni-H POTCAR changed')
    blocks=split_potcar(pot)
    require(set(blocks)==set(PLAN['potential']['dataset_labels']),'Potential element order/identities mismatch')
    hfile=h2/'02_static/POTCAR'
    require(hfile.is_file() and hfile.read_bytes()==blocks['H'],'H2 potential differs from reference MgNiH potential')
    return base,blocks,{'source_sha256':sha(src),'source_force_eV_A':static['max_force_eV_A'],
                        'H2_reference_F_eV':v['F_eV'],'POTCAR_SHA256':sha(pot)}
def incar(name,n,stage):
    relax=stage=='01_relax'
    rows=[
        'SYSTEM = '+name,
        'PREC = Accurate','ENCUT = 520','ISPIN = 2',
        'MAGMOM = 24*0.0 12*1.0'+(' '+str(n)+'*0.0' if n else ''),
        'LDAU = .FALSE.','ALGO = Normal','ISMEAR = 0','SIGMA = 0.05',
        'ISYM = 0','LREAL = .FALSE.','LASPH = .TRUE.','ADDGRID = .TRUE.',
        'EDIFF = 1E-7','NELM = 200','NCORE = 4','ISTART = 0','ICHARG = 2',
        'ISIF = 3' if relax else 'ISIF = 2',
        'IBRION = 2' if relax else 'IBRION = -1',
        'NSW = 240' if relax else 'NSW = 0',
        'POTIM = 0.12','EDIFFG = -0.015',
        'LWAVE = .FALSE.','LCHARG = .FALSE.',
        'LORBIT = 11' if not relax else 'LORBIT = 0',
    ]
    return '\n'.join(rows)+'\n'
def kpoints():
    return 'Gamma centered bulk\n0\nGamma\n5 5 2\n0 0 0\n'
def write_inputs(base,blocks,preview=False):
    ROOT.mkdir(parents=True,exist_ok=False);(ROOT/'logs').mkdir()
    for meta in PLAN['models']:
        name=meta['name'];struct,siteinfo=model_structure(base,meta)
        folder=ROOT/name;folder.mkdir()
        els=list(dict.fromkeys(struct['species']))
        pbytes=b''.join(blocks[e] for e in els) if blocks else None
        for stage in ('01_relax','02_static'):
            d=folder/stage;d.mkdir()
            (d/'INCAR').write_text(incar(name,len(meta['h_sites']),stage))
            (d/'KPOINTS').write_text(kpoints())
            if pbytes is not None:(d/'POTCAR').write_bytes(pbytes)
        (folder/'01_relax/POSCAR').write_text(poscar(struct,name+' | candidate T31/T22 tetrahedral sites'))
        save(folder/'model.json',{'model':name,'composition':{e:struct['species'].count(e) for e in els},
                                  'initial':siteinfo,'source_18_atom_sha':PLAN['source']['sha256'],
                                  'not_proven_global_ground_state':True})
        print('PREPARED_MODEL|'+name+'|atoms='+str(len(struct['species']))+'|x='+str(siteinfo['nominal_x'])+
              '|min_HH='+str(siteinfo['min_initial_HH_A']),flush=True)
    shutil.copy2(HERE/'plan.json',ROOT/'plan.json')
    shutil.copy2(Path(__file__),ROOT/'screen.py')
    job='''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
ulimit -s unlimited
ulimit -l unlimited
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 screen.py --run-index "$SLURM_ARRAY_TASK_ID"
'''
    (ROOT/'lowH.slurm').write_text(job)
    subprocess.run(['bash','-n',str(ROOT/'lowH.slurm')],check=True)
    hashes={str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file()}
    save(ROOT/'manifest.json',{'prepared_at':datetime.datetime.now().astimezone().isoformat(),
                               'preview_without_pseudopotentials':preview,'names':NAMES,'input_hashes':hashes,
                               'number_of_models':5,'source_H2_eV':PLAN['h2_reference']['expected_F_eV']})
    print('PREPARED_BATCH|'+str(ROOT)+'|five configurations|no NEB and no large model',flush=True)
def prepare(preview=False,source=None):
    if preview:
        require(source is not None,'Preview requires --source-poscar')
        base=parse(Path(source).read_text())
        write_inputs(base,blocks=None,preview=True)
    else:
        base,blocks,checks=validate_external()
        if ROOT.exists():
            require((ROOT/'manifest.json').is_file(),'Existing incomplete production directory; do not overwrite')
            print('ALREADY_PREPARED|'+str(ROOT),flush=True);return
        write_inputs(base,blocks,preview=False)
        save(ROOT/'source_verification.json',checks)
        print('ACTUAL_SOURCE_CHECK|'+json.dumps(checks),flush=True)
def run_model(idx):
    require(0<=idx<len(NAMES),'Invalid array index')
    require(int(os.environ.get('SLURM_NTASKS','0'))==32,'VASP only in prescribed Slurm allocation')
    name=NAMES[idx];d=ROOT/name
    if (d/'VERIFIED_FINISHED').exists():print('ALREADY_VERIFIED|'+name);return
    (d/'.execution_started').mkdir(exist_ok=False)
    expected=parse((d/'01_relax/POSCAR').read_text());n=len(expected['species'])
    for stage in ('01_relax','02_static'):
        p=d/stage
        if stage=='02_static':shutil.copy2(d/'01_relax/CONTCAR',p/'POSCAR')
        with (p/'vasp.stdout').open('w') as out,(p/'vasp.stderr').open('w') as err:
            proc=subprocess.run(['mpirun','-np','32','vasp_std'],cwd=p,stdout=out,stderr=err)
        require(proc.returncode==0,'VASP returned nonzero '+name+' '+stage)
        v=read_outcar(p/'OUTCAR');s=parse((p/'CONTCAR').read_text())
        require(s['species']==expected['species'] and len(s['frac'])==n,'Atom composition or order changed '+name)
        require(v['force_count']==n,'Not all final forces found '+name+' '+stage)
        cv=lambda c: abs(dot(c[0],[c[1][1]*c[2][2]-c[1][2]*c[2][1],c[1][2]*c[2][0]-c[1][0]*c[2][2],c[1][0]*c[2][1]-c[1][1]*c[2][0]]))
        volume_ratio=cv(s['cell'])/cv(expected['cell'])
        require(.8<volume_ratio<1.2,'Unexpected bulk volume shift '+name)
        limit=.015 if stage=='01_relax' else .020
        v['force_gate_eV_A']=limit;v['cell_volume_ratio_to_input']=volume_ratio
        v['verified']=bool(v['normal_end'] and v['last_scf_ediff'] and
                           v['F_eV'] is not None and v['E0_eV'] is not None and
                           v['max_force_eV_A'] is not None and v['max_force_eV_A']<limit and
                           (stage=='02_static' or v['ionic_converged']))
        hs=[j for j,e in enumerate(s['species']) if e=='H']
        v['H_nearest_metals']=[{'index':j+1,'nearest4':nearest_metal(s,s['frac'][j])[:4]}
                               for j in hs]
        v['min_HH_A']=min((mic_distance(s['frac'][i],s['frac'][j],s['cell'])
                          for i,j in itertools.combinations(hs,2)),default=None)
        save(p/'validation.json',v)
        print('STAGE_RESULT|'+name+'|'+stage+'|'+json.dumps(v),flush=True)
        require(v['verified'],'VASP numerical gate failed '+name+' '+stage)
    (d/'VERIFIED_FINISHED').touch()
def submit():
    require(sys.platform.startswith('linux') and shutil.which('sbatch'),'Must submit inside the Slurm cluster')
    if (ROOT/'submission.json').is_file():
        print('ALREADY_SUBMITTED|'+(ROOT/'submission.json').read_text());return
    if not ROOT.exists():prepare(False)
    manifest=json.loads((ROOT/'manifest.json').read_text())
    require(not manifest['preview_without_pseudopotentials'],'Refuse to submit preview without licensed POTCAR')
    for rel,checksum in manifest['input_hashes'].items():
        require(sha(ROOT/rel)==checksum,'Pinned production input changed '+rel)
    require(not (ROOT/'.submission_attempted').exists(),'Uncertain earlier submission, inspect before retry')
    # No duplicate older 147-atom interface work; only this five-model array.
    require(manifest['number_of_models']==5 and manifest['names']==NAMES,'Model roster changed')
    (ROOT/'.submission_attempted').mkdir()
    p=subprocess.run(['sbatch','--parsable','--job-name=mg2ni_lowH_bulk_small','lowH.slurm'],
                     cwd=ROOT,capture_output=True,text=True,timeout=45)
    require(p.returncode==0,'sbatch failed: '+p.stderr[:500])
    jobid=p.stdout.strip().split(';')[0]
    require(re.fullmatch(r'\d+',jobid),'Unexpected sbatch response '+p.stdout[:100])
    record={'job_id':jobid,'models':NAMES,'array':'0-4',
            'submitted_at':datetime.datetime.now().astimezone().isoformat(),
            'Hx_147_atom_interface_resubmitted':False,'no_NEB':True}
    save(ROOT/'submission.json',record)
    print('SUBMITTED|lowH_bulk_screen|'+jobid+'|array 0-4',flush=True)
def status():
    print('CHECKED_AT|'+datetime.datetime.now().astimezone().isoformat())
    if not (ROOT/'submission.json').is_file():print('NOT_SUBMITTED');return
    rec=json.loads((ROOT/'submission.json').read_text())
    print('SUBMISSION|'+json.dumps(rec))
    jobid=rec['job_id']
    for cmd in [['squeue','-r','-h','-j',jobid,'-o','%i|%j|%T|%M|%R'],
                ['sacct','-X','-n','-P','-j',jobid,'-o','JobID,State,ExitCode,Elapsed,Start,End']]:
        z=subprocess.run(cmd,text=True,capture_output=True,timeout=30)
        print('QUEUE_OR_ACCOUNTING|'+z.stdout.strip())
    for name in NAMES:
        p=ROOT/name/'02_static/validation.json'
        label='VERIFIED_FINISHED' if (ROOT/name/'VERIFIED_FINISHED').exists() else 'NOT_VERIFIED'
        record=json.loads(p.read_text()) if p.exists() else None
        print('MODEL|'+name+'|'+label+'|'+json.dumps(record))
def main():
    global ROOT
    q=argparse.ArgumentParser()
    group=q.add_mutually_exclusive_group(required=True)
    group.add_argument('--preview',action='store_true')
    group.add_argument('--submit',action='store_true')
    group.add_argument('--status',action='store_true')
    group.add_argument('--run-index',type=int)
    q.add_argument('--root',type=Path)
    q.add_argument('--source-poscar',type=Path)
    args=q.parse_args()
    if args.root:
        require(args.preview,'--root only allowed for preview')
        ROOT=args.root.resolve()
    if args.preview:prepare(True,args.source_poscar)
    elif args.submit:submit()
    elif args.status:status()
    else:run_model(args.run_index)
if __name__=='__main__':
    try:main()
    except Exception as e:
        print('STOP|'+str(e),file=sys.stderr,flush=True)
        raise
