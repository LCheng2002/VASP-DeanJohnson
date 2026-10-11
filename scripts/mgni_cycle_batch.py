#!/usr/bin/env python3
"""Bounded Mg2Ni/Mg2NiH4 endpoint screen; no NEB and no original-file edits.
Prepare 9 solids + H2. Submit one solids array, one small H2 job and
one afterok continuation that constructs exactly S2 then S3 from relaxed S1.
"""
import argparse, copy, datetime, hashlib, itertools, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
HOME=Path.home()
SOURCE=HOME/'Mg2Ni_Mg_H_thermo_validation_20261010'
ROOT=HOME/'MgNi_abs_des_cycle_20261011'
NAMES=['D_S1a_same_Ni','D_S1b_split_Ni','D_Sfar_Mg','A_H2_Ni_parallel','A_H2_Ni_vertical','A_H2_Mg_parallel','A_H2_Mg_vertical','A_2H_Ni','A_2H_Mg']
PARENTS={'A':('A0_Mg_Mg2Ni_clean','90c33aa81b094fe82a157de942ed87b71c7b874010ab0a334b30b9289a6c71f8'), 'D':('D0_MgH2_Mg2NiH4_clean','a385f495e97cc8a08d3d72bfbbc5c0c9ac9f594f6e843079f4ac949b5999ec7f')}
POT_REL='A-H-Ni_Mg_Mg2Ni/01_relax/POTCAR'
POT_SHA='72b607bd5434dcb71996bb13fd478af3592115762200d3c53e0f47c68b0a84ae'
SHIFTS=list(itertools.product((-1,0,1),repeat=3))
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def norm(v):return math.sqrt(dot(v,v))
def mic(a,b,c):
    d=[x-y for x,y in zip(a,b)];det=dot(c[0],cross(c[1],c[2]));require(abs(det)>1e-8,'Singular cell')
    f=[dot(d,cross(c[1],c[2]))/det,dot(d,cross(c[2],c[0]))/det,dot(d,cross(c[0],c[1]))/det]
    r=[round(x) for x in f]
    vs=[[d[j]-sum((r[k]+s[k])*c[k][j] for k in range(3)) for j in range(3)] for s in SHIFTS]
    return min(vs,key=dotself)
def dotself(v):return dot(v,v)
def distance(a,b,c):return norm(mic(a,b,c))
def parse(s):
    l=s.splitlines();require(len(l)>8,'Incomplete POSCAR');scale=float(l[1]);require(scale>0,'Positive single scale required')
    c=[[float(v)*scale for v in x.split()[:3]] for x in l[2:5]]
    els=l[5].split();counts=list(map(int,l[6].split()));sp=sum(([e]*n for e,n in zip(els,counts)),[])
    i=7;select=l[i].strip().lower().startswith('s');i+=int(select);direct=l[i].strip().lower().startswith('d');i+=1
    xyz=[];mask=[]
    for row in l[i:i+len(sp)]:
        q=row.split();v=list(map(float,q[:3]));xyz.append([sum(v[k]*c[k][j] for k in range(3)) for j in range(3)] if direct else [x*scale for x in v]);mask.append([x.upper().startswith('T') for x in q[3:6]] if select else [True]*3)
    require(len(xyz)==len(sp) and all(len(x)==3 for x in mask),'Bad atom or mask count')
    return {'cell':c,'species':sp,'xyz':xyz,'mask':mask,'ids':list(range(1,len(sp)+1))}
def poscar(s,title):
    els=list(dict.fromkeys(s['species']));require(sum(([e]*s['species'].count(e) for e in els),[])==s['species'],'Elements not grouped')
    lines=[title,'1.0']+[' '.join('%.14f'%v for v in row) for row in s['cell']]+[' '.join(els),' '.join(str(s['species'].count(e)) for e in els),'Selective dynamics','Cartesian']
    lines+=[' '.join('%.14f'%x for x in v)+' '+' '.join('T' if x else 'F' for x in m) for v,m in zip(s['xyz'],s['mask'])]
    return '\n'.join(lines)+'\n'
def nearest(s,p,element=None):
    return sorted((distance(p,x,s['cell']),e,i) for x,e,i in zip(s['xyz'],s['species'],s['ids']) if element is None or e==element)
def patch(raw,changes):
    lines=[]
    for line in raw.splitlines():
        part=re.split('[!#]',line)[0];k=part.split('=',1)[0].strip().upper()
        if k not in changes:lines.append(line)
    return '\n'.join(lines)+'\n'+'\n'.join(k+' = '+str(v) for k,v in changes.items())+'\n'
def config(raw):
    d={}
    for l in raw.splitlines():
        for x in re.split(';',re.split('[!#]',l)[0]):
            if '=' in x:
                k,v=x.split('=',1);d[k.strip().upper()]=v.strip()
    return d

def make_specs(a,d):
    require(a['species'].count('Mg')==36 and a['species'].count('Ni')==6 and len(a['ids'])==42,'A0 composition changed')
    require(d['species']==['Mg']*20+['Ni']*8+['H']*40,'D0 composition changed')
    out={}
    for name,deleted in [('D_S1a_same_Ni',[43,45]),('D_S1b_split_Ni',[43,42]),('D_Sfar_Mg',[29,31])]:
        s=copy.deepcopy(d);keep=[i for i,x in enumerate(s['ids']) if x not in deleted]
        for k in ['species','xyz','mask','ids']:s[k]=[s[k][i] for i in keep]
        vacancies=[]
        for oldid in deleted:
            j=d['ids'].index(oldid);p=d['xyz'][j];dist,_,ni=nearest(d,p,'Ni')[0]
            vacancies.append({'parent_H_id':oldid,'parent_xyz':p,'anchor_Ni_id':ni,'anchor_distance_A':dist,'Ni_H_vector':mic(p,d['xyz'][d['ids'].index(ni)],d['cell'])})
        if name=='D_S1a_same_Ni':require([v['anchor_Ni_id'] for v in vacancies]==[22,22],'Same-Ni sites invalid')
        if name=='D_S1b_split_Ni':require([v['anchor_Ni_id'] for v in vacancies]==[22,21],'Split-Ni sites invalid')
        out[name]={'structure':s,'source':'D','kind':'double_H_vacancy_endpoint','deleted_H_parent_ids':deleted,'vacancies':vacancies,'gas_release_path_not_modelled':True}
    for element,anchor in [('Ni',41),('Mg',32)]:
        require(a['species'][anchor-1]==element,'A anchor identity changed')
        for orient in ['parallel','vertical','2H']:
            sep=2.0 if orient=='2H' else .75; height=1.65 if orient=='2H' else 2.30
            vec=[0,sep/2,0] if orient!='vertical' else [0,0,sep/2]
            center=[a['xyz'][anchor-1][0],a['xyz'][anchor-1][1],a['xyz'][anchor-1][2]+height]
            pts=[[center[j]+sign*vec[j] for j in range(3)] for sign in (1,-1)]
            name='A_2H_'+element if orient=='2H' else 'A_H2_'+element+'_'+orient
            for p in pts:require(nearest(a,p)[0][0]>1.4,'Initial H-metal collision '+name)
            require(abs(distance(pts[0],pts[1],a['cell'])-sep)<1e-6,'H-H geometry invalid')
            s=copy.deepcopy(a);s['xyz']+=pts;s['species']+=['H','H'];s['mask']+=[[True]*3,[True]*3];s['ids']+=[43,44]
            out[name]={'structure':s,'source':'A','kind':'2H_endpoint' if orient=='2H' else 'molecular_H2_start','anchor_parent_id':anchor,'initial_H_H_A':sep,'initial_new_H_neighbours':[nearest(a,p)[:5] for p in pts],'surface_comparison':'Ni versus Mg-rich sites of same exposed interface slab; not isolated pure Mg reference'}
    require(set(out)==set(NAMES),'Model roster mismatch')
    for name,v in out.items():
        s=v['structure'];n=66 if name.startswith('D') else 44
        require(len(s['ids'])==n and len(set(s['ids']))==n,'Atom mapping invalid '+name)
        require(s['species'].count('H')==(38 if name.startswith('D') else 2),'H balance invalid')
    return out

def read_output(folder):
    p=Path(folder);s=parse((p/'POSCAR').read_text());fval=e0=None;footer=ionic=ediff=False;it=None;reading=False;buf=[];forces=[];steps=0
    with (p/'OUTCAR').open(errors='replace') as f:
        for line in f:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and it!=m.group(1):it=m.group(1);ediff=False
            if 'aborting loop because EDIFF is reached' in line:ediff=True
            if 'reached required accuracy' in line:ionic=True
            if 'General timing and accounting informations for this job' in line:footer=True
            if 'free  energy   TOTEN' in line:fval=float(line.split('=')[1].split()[0]);steps+=1
            if 'energy(sigma->0)' in line:e0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buf=[];continue
            if reading:
                q=line.split()
                if len(q)==6:
                    try:buf.append(list(map(float,q[3:6])));continue
                    except ValueError:pass
                if buf:forces=buf;reading=False
    require(len(forces)==len(s['ids']),'Missing final forces '+str(p))
    fm=max(math.sqrt(sum(x*x for x,m in zip(v,mask) if m)) for v,mask in zip(forces,s['mask']))
    return {'F_eV':fval,'E0_eV':e0,'last_scf_ediff':ediff,'normal_end':footer,'ionic_converged':ionic,'max_free_force_eV_A':fm,'ionic_steps':steps}

def write_model(name,spec,templates,pot):
    path=ROOT/name;path.mkdir(exist_ok=False)
    source=spec['source'];s=spec['structure'];inp=poscar(s,name+' | endpoint screen, not a reaction trajectory')
    for st in ['01_relax','02_static']:
        p=path/st;p.mkdir();(p/'POTCAR').write_bytes(pot)
        inc=templates[source][st]['incar'];kp=templates[source][st]['kpoints']
        if st=='01_relax':inc=patch(inc,{'NSW':100 if source=='H' else 400,'POTIM':'0.10','ISTART':0,'ICHARG':2})
        (p/'INCAR').write_text(inc);(p/'KPOINTS').write_text(kp)
    (path/'01_relax/POSCAR').write_text(inp)
    meta={k:v for k,v in spec.items() if k!='structure'};meta.update({'model':name,'parent_atom_ids':s['ids'],'composition':{e:s['species'].count(e) for e in dict.fromkeys(s['species'])},'cell_A':s['cell'],'new_structure_sha256':sha(path/'01_relax/POSCAR')})
    (path/'model.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False))

def job_text(mode,ranks=64,wall='24:00:00',array=False):
    text='#!/bin/bash\n#SBATCH --partition=64c512g\n#SBATCH --nodes=1\n#SBATCH --ntasks='+str(ranks)+'\n#SBATCH --time='+wall+'\n'
    text+= '#SBATCH --exclusive\n' if ranks==64 else '#SBATCH --mem=4G\n'
    if array:text+='#SBATCH --array=0-8%2\n'
    text+='#SBATCH --output=logs/%x-%A_%a.out\n#SBATCH --error=logs/%x-%A_%a.err\nset -euo pipefail\ncd "$SLURM_SUBMIT_DIR"\nmodule purge\nmodule load oneapi vasp/6.3.2-intel-2021.4.0\nulimit -s unlimited\nulimit -l unlimited\nexport OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\n'
    text+='python3 mgni_cycle_batch.py '+mode+'\n';return text

def prepare():
    if (ROOT/'manifest.json').is_file():
        verify_inputs();return json.loads((ROOT/'manifest.json').read_text())
    require(not ROOT.exists(),'An incomplete target directory already exists; inspect before retry')
    data={};templates={};refs={}
    for key,(name,h) in PARENTS.items():
        p=SOURCE/name;require(sha(p/'01_relax/CONTCAR')==h,'Parent geometry changed '+key)
        data[key]=parse((p/'01_relax/CONTCAR').read_text());templates[key]={}
        for st in ['01_relax','02_static']:
            q=read_output(p/st);require(q['last_scf_ediff'] and q['normal_end'] and q['max_free_force_eV_A']<.0200001,'Parent numerical check failed '+key+'/'+st)
            if st=='01_relax':require(q['ionic_converged'],'Parent relaxation unconverged')
            templates[key][st]={'incar':(p/st/'INCAR').read_text(),'kpoints':(p/st/'KPOINTS').read_text()}
            if st=='02_static':refs[key]=q
    require(sha(SOURCE/POT_REL)==POT_SHA,'Mg Ni H potential changed')
    require(sha(SOURCE/PARENTS['D'][0]/'02_static/POTCAR')==POT_SHA,'A D combined potential mismatch')
    pot=(SOURCE/POT_REL).read_bytes();title=[l for l in pot.decode(errors='replace').splitlines() if 'TITEL' in l]
    require(len(title)==3 and all(re.search(r'PAW_PBE\s+'+e+r'\s',l) for e,l in zip(['Mg','Ni','H'],title)),'Potential order invalid')
    blocks=[];current=[]
    for l in pot.splitlines(keepends=True):
        current.append(l)
        if l.strip()==b'End of Dataset':blocks.append(b''.join(current));current=[]
    require(len(blocks)==3,'Cannot split licensed potential locally')
    hp=blocks[2];require(re.search(rb'TITEL\s*=\s*PAW_PBE\s+H\s',hp),'H potential block invalid')
    specs=make_specs(data['A'],data['D'])
    ROOT.mkdir();(ROOT/'logs').mkdir();(ROOT/'reference').mkdir()
    for key in data:(ROOT/'reference'/('parent_'+key+'.vasp')).write_text(poscar(data[key],key+' optimized reference'))
    (ROOT/'reference/templates.json').write_text(json.dumps(templates,indent=2))
    (ROOT/'reference/energies.json').write_text(json.dumps(refs,indent=2))
    for name in NAMES:write_model(name,specs[name],templates,pot)
    mol={'cell':[[20.,0,0],[0,20.,0],[0,0,20.]],'species':['H','H'],'xyz':[[10,10,9.625],[10,10,10.375]],'mask':[[True]*3,[True]*3],'ids':[1,2]}
    htemp={'H':{}}
    for st in ['01_relax','02_static']:
        changes={'NCORE':1,'ISPIN':1,'LREAL':'.FALSE.','EDIFF':'1E-8','ISMEAR':0,'SIGMA':'0.05','LORBIT':0,'EDIFFG':'-0.005'}
        htemp['H'][st]={'incar':patch(templates['A'][st]['incar'],changes),'kpoints':'H2 box reference\n0\nGamma\n1 1 1\n0 0 0\n'}
    write_model('H2_reference',{'structure':mol,'source':'H','kind':'singlet_H2_reference','box_size_convergence_not_claimed':True},htemp,hp)
    script=Path(__file__).resolve();shutil.copy2(script,ROOT/'mgni_cycle_batch.py')
    for fname,mode,n,wall,arr in [('solids.slurm','--run "$SLURM_ARRAY_TASK_ID"',64,'24:00:00',True),('h2.slurm','--run 9',4,'00:30:00',False),('refill.slurm','--refill',64,'24:00:00',False)]:
        p=ROOT/fname;p.write_text(job_text(mode,n,wall,arr));subprocess.run(['bash','-n',str(p)],check=True)
    manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_root':str(SOURCE),'parent_sha256':{k:v[1] for k,v in PARENTS.items()},'potential_sha256':POT_SHA,'H_potential_sha256':hashlib.sha256(hp).hexdigest(),'array_indices':dict(enumerate(NAMES)),'H2_job_separate':True,'deferred_refill':'afterok S1a and S1b; exactly S2 then S3','screen_scope':'endpoint thermodynamics only, no kinetic acceleration claim','gas_exit_caveat':'D missing Ni-H sites are buried near interface. Removing 2 H references gas reservoir; it does not model H2 escaping through solid. D far H29 H31 are bottom-surface controls, not bulk.','input_hashes':{}}
    for p in ROOT.rglob('*'):
        if p.is_file():manifest['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
    (ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2));print('PREPARED|9 solid configurations and H2; S2/S3 deferred until S1 relaxation',flush=True)
    return manifest

def verify_inputs():
    m=json.loads((ROOT/'manifest.json').read_text())
    for rel,h in m['input_hashes'].items():require(sha(ROOT/rel)==h,'Pinned input changed: '+rel)
    return m

def run_model(name):
    path=ROOT/name;require(path.is_dir(),'Missing model '+name)
    if (path/'VERIFIED_FINISHED').is_file():return
    (path/'.execution_lock').mkdir(exist_ok=False)
    rank=int(os.environ.get('SLURM_NTASKS','0'));require(rank in [4,64],'Run only within scheduled Slurm allocation')
    for st in ['01_relax','02_static']:
        p=path/st
        if st=='02_static':shutil.copy2(path/'01_relax/CONTCAR',p/'POSCAR')
        with (p/'vasp.stdout').open('w') as out,(p/'vasp.stderr').open('w') as err:
            proc=subprocess.run(['mpirun','-np',str(rank),'vasp_std'],cwd=p,stdout=out,stderr=err)
        require(proc.returncode==0,'VASP execution failed '+name+'/'+st)
        r=read_output(p);r['phase']=st
        r['verified']=r['normal_end'] and r['last_scf_ediff'] and r['F_eV'] is not None and r['max_free_force_eV_A']<.0200001 and (st=='02_static' or r['ionic_converged'])
        if (p/'CONTCAR').is_file():
            s=parse((p/'CONTCAR').read_text());hs=[i for i,e in enumerate(s['species']) if e=='H']
            r['min_H_H_A']=min(distance(s['xyz'][i],s['xyz'][j],s['cell']) for i,j in itertools.combinations(hs,2)) if len(hs)>1 else None
            r['H_nearest_metals']=[{'atom':i+1,'Mg':nearest(s,s['xyz'][i],'Mg')[:2],'Ni':nearest(s,s['xyz'][i],'Ni')[:2]} for i in hs]
        (p/'validation.json').write_text(json.dumps(r,indent=2));print('STAGE_RESULT|'+name+'|'+json.dumps(r),flush=True)
        require(r['verified'],'Convergence check failed '+name+'/'+st)
    (path/'VERIFIED_FINISHED').touch()

def build_refill(source_name,new_name,already_filled):
    source=ROOT/source_name;require((source/'VERIFIED_FINISHED').is_file(),'Refill parent must pass validation')
    base=parse((ROOT/'reference/parent_D.vasp').read_text());meta=json.loads((source/'model.json').read_text())
    s=parse((source/'02_static/CONTCAR').read_text());s['ids']=meta['parent_atom_ids'];require(len(s['ids'])==66 and s['species']==['Mg']*20+['Ni']*8+['H']*38,'Refill composition changed')
    vacancies=meta['vacancies'];moved=meta.get('transferred_parent_H_ids',[])
    remaining=[v for v in vacancies if v['parent_H_id'] not in already_filled];candidates=[]
    for vacancy in remaining:
        ni=vacancy['anchor_Ni_id'];j=s['ids'].index(ni)
        displacement=distance(s['xyz'][j],base['xyz'][ni-1],s['cell'])
        require(displacement<1.2,'Ni reconstruction too large for translated-site refill; manual reassessment required')
        target=[s['xyz'][j][k]+vacancy['Ni_H_vector'][k] for k in range(3)]
        for donor in range(29,37):
            if donor not in s['ids'] or donor in moved:continue
            dj=s['ids'].index(donor)
            if nearest(s,s['xyz'][dj],'Ni')[0][0]<2.10:continue
            metals=[distance(target,p,s['cell']) for p,e in zip(s['xyz'],s['species']) if e!='H']
            hh=[distance(target,p,s['cell']) for k,(p,e) in enumerate(zip(s['xyz'],s['species'])) if e=='H' and k!=dj]
            if min(metals)<1.10 or min(hh)<1.05:continue
            candidates.append((distance(target,s['xyz'][dj],s['cell']),donor,vacancy,target))
    require(candidates,'No safe actual-MgH2 donor/vacancy pair; no forced refill will run')
    travel,donor,vacancy,target=min(candidates,key=lambda x:x[0]);s['xyz'][s['ids'].index(donor)]=target
    require(s['species'].count('H')==38,'Refill must move H, never add H')
    templates=json.loads((ROOT/'reference/templates.json').read_text());pot=(ROOT/'D_S1a_same_Ni/01_relax/POTCAR').read_bytes()
    spec={'structure':s,'source':'D','kind':'one_H_cross_interface_refill_endpoint','vacancies':vacancies,'filled_vacancy_parent_ids':already_filled+[vacancy['parent_H_id']],'transferred_parent_H_ids':moved+[donor],'refill_from_model':source_name,'last_transfer':{'donor_parent_H_id':donor,'acceptor_vacancy_parent_H_id':vacancy['parent_H_id'],'endpoint_displacement_A':travel,'new_xyz':target},'not_a_diffusion_barrier':True}
    write_model(new_name,spec,templates,pot)
    print('REFILL_PREPARED|'+new_name+'|'+json.dumps(spec['last_transfer']),flush=True)
    return spec['filled_vacancy_parent_ids']

def refill():
    require((ROOT/'D_S1a_same_Ni/VERIFIED_FINISHED').exists() and (ROOT/'D_S1b_split_Ni/VERIFIED_FINISHED').exists(),'Both S1 candidates must converge')
    r={n:json.loads((ROOT/n/'02_static/validation.json').read_text()) for n in NAMES[:2]}
    chosen=min(r,key=lambda x:r[x]['F_eV']);chosen0=min(r,key=lambda x:r[x]['E0_eV'])
    require(chosen==chosen0,'F/E0 select different S1 parents; do not build an automatic mechanistic sequence')
    report={'selected_S1':chosen,'candidate_energies':r,'selection':'lowest of two sampled converged endpoints, not global ground state'}
    (ROOT/'refill_selection.json').write_text(json.dumps(report,indent=2))
    filled=build_refill(chosen,'D_S2_refill_one',[]);run_model('D_S2_refill_one')
    build_refill('D_S2_refill_one','D_S3_refill_two',filled);run_model('D_S3_refill_two')

def submit():
    prepare()
    ledger=ROOT/'submission.json'
    if ledger.is_file():
        old=json.loads(ledger.read_text())
        if old.get('all_submitted'):print('ALREADY_SUBMITTED|'+json.dumps(old));return
    require(not (ROOT/'.submission_lock').exists(),'Submission result uncertain; inspect existing queue and ledger, do not retry blindly')
    q=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%j'],capture_output=True,text=True,check=True,timeout=30).stdout
    require(not any(x.startswith('mgni_cycle_') for x in q.splitlines()),'Matching cycle jobs already active')
    verify_inputs();(ROOT/'.submission_lock').mkdir();records={}
    def send(key,args):
        proc=subprocess.run(['sbatch','--parsable']+args,cwd=ROOT,capture_output=True,text=True,timeout=45)
        require(proc.returncode==0,'sbatch failed '+key+': '+proc.stderr[:300]);jid=proc.stdout.strip()
        require(re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?',jid),'Unexpected sbatch response')
        records[key]=jid;ledger.write_text(json.dumps(records,indent=2));print('SUBMITTED|'+key+'|'+jid,flush=True);return jid.split(';')[0]
    array=send('solids_array',['--job-name=mgni_cycle_screen','solids.slurm'])
    send('H2_reference',['--job-name=mgni_cycle_H2','h2.slurm'])
    send('refill_chain',['--job-name=mgni_cycle_refill','--dependency=afterok:'+array+'_0:'+array+'_1','--kill-on-invalid-dep=yes','refill.slurm'])
    records.update({'all_submitted':True,'submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'array_indices':dict(enumerate(NAMES)),'root':str(ROOT),'max_simultaneous_solid_models':3,'refill_expected_models':['D_S2_refill_one','D_S3_refill_two']});ledger.write_text(json.dumps(records,indent=2));print('SUBMISSION_COMPLETE_JSON|'+json.dumps(records),flush=True)

def status():
    r=json.loads((ROOT/'submission.json').read_text());print('SUBMISSION_RECORD|'+json.dumps(r),flush=True)
    ids=[r[k].split(';')[0] for k in ['solids_array','H2_reference','refill_chain'] if k in r]
    if ids:
        subprocess.run(['squeue','-r','-j',','.join(ids),'-h','-o','%i|%j|%T|%M|%R'],check=False)
        subprocess.run(['sacct','-X','-n','-P','-j',','.join(ids),'-o','JobID,State,ExitCode,Elapsed,Start,End'],check=False)
    for name in NAMES+['H2_reference','D_S2_refill_one','D_S3_refill_two']:
        p=ROOT/name
        if not p.is_dir():print('MODEL|'+name+'|NOT_YET_CONSTRUCTED');continue
        phase=[]
        for st in ['01_relax','02_static']:
            if (p/st/'validation.json').is_file():phase.append(st+':'+(p/st/'validation.json').read_text().replace('\n',''))
        print('MODEL|'+name+'|'+('VERIFIED_FINISHED' if (p/'VERIFIED_FINISHED').exists() else 'NOT_FINISHED')+'|'+str(phase))

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--submit',action='store_true');g.add_argument('--status',action='store_true');g.add_argument('--run',type=int);g.add_argument('--refill',action='store_true');a=p.parse_args()
    if a.prepare:prepare()
    elif a.submit:submit()
    elif a.status:status()
    elif a.refill:refill()
    elif a.run is not None:
        require(0<=a.run<=9,'Invalid model index');run_model((NAMES+['H2_reference'])[a.run])
if __name__=='__main__':
    try:main()
    except Exception as exc:print('STOP|'+str(exc),file=sys.stderr,flush=True);raise
