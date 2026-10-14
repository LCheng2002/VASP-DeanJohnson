#!/usr/bin/env python3
"""Fresh Ca/Mg water-OH screening; no old geometry or energy import.

Preview: generate all 23 *unrelaxed* inputs from an ideal lattice.
Production: three new references -> two clean slabs -> own-parent dependent
adsorbate jobs. Only licensed local potentials may be reused, never structures.
"""
import argparse, copy, datetime, hashlib, itertools, json, math, os, re
import shlex, shutil, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())
HOME=Path.home()
ROOT=HOME/PLAN['batch']
SHIFTS=list(itertools.product((-1,0,1),repeat=3))
ORDER=['Mg','Ca','O','H']
FLOAT=r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?'

def need(test,message):
    if not test: raise RuntimeError(message)
def timestamp():return datetime.datetime.now().astimezone().isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    p=Path(p);t=p.with_name(p.name+'.tmp')
    t.write_text(json.dumps(d,indent=2,ensure_ascii=False));t.replace(p)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def norm(a):return math.sqrt(dot(a,a))
def plus(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,t):return [x*t for x in a]
def cart(f,c):return [sum(f[k]*c[k][j] for k in range(3)) for j in range(3)]
def frac(x,c):
    det=dot(c[0],cross(c[1],c[2]));need(abs(det)>1e-8,'Singular cell')
    return [dot(x,cross(c[1],c[2]))/det,dot(x,cross(c[2],c[0]))/det,dot(x,cross(c[0],c[1]))/det]
def mic(a,b,c):
    f=frac([x-y for x,y in zip(a,b)],c);r=[round(v) for v in f]
    return min((cart([f[k]-r[k]-d[k] for k in range(3)],c) for d in SHIFTS),key=lambda x:dot(x,x))
def dist(a,b,c):return norm(mic(a,b,c))
def volume(c):return abs(dot(c[0],cross(c[1],c[2])))
def atom(s,label):return s['xyz'][s['labels'].index(label)]
def labels_metal(s):return [l for l,e in zip(s['labels'],s['species']) if e in ('Mg','Ca')]
def grouped(s):
    ix=sorted(range(len(s['species'])),key=lambda i:ORDER.index(s['species'][i]))
    return {**copy.deepcopy(s),**{k:[copy.deepcopy(s[k][i]) for i in ix] for k in ('species','xyz','mask','labels')}}
def parse(text):
    ls=text.splitlines();need(len(ls)>=9,'Incomplete POSCAR')
    scale=float(ls[1]);need(scale>0,'Positive scale required')
    c=[[float(t)*scale for t in l.split()[:3]] for l in ls[2:5]]
    sp=sum(([e]*int(n) for e,n in zip(ls[5].split(),ls[6].split())),[])
    i=7;select=ls[i].strip().lower().startswith('s');i+=int(select)
    direct=ls[i].strip().lower().startswith('d');i+=1
    xyz=[];masks=[]
    for l in ls[i:i+len(sp)]:
        q=l.split();v=list(map(float,q[:3]));need(len(v)==3,'Incomplete atom row')
        xyz.append(cart(v,c) if direct else mul(v,scale))
        masks.append([v.upper().startswith('T') for v in q[3:6]] if select else [True]*3)
    need(len(xyz)==len(sp) and all(len(m)==3 for m in masks),'Atom/mask count')
    return {'cell':c,'species':sp,'xyz':xyz,'mask':masks,'labels':list(range(len(sp)))}
def poscar(s,name):
    els=list(dict.fromkeys(s['species']))
    need(sum(([e]*s['species'].count(e) for e in els),[])==s['species'],'Group species first')
    rows=[name,'1.0']+[' '.join('%.14f'%v for v in r) for r in s['cell']]
    rows+=[' '.join(els),' '.join(str(s['species'].count(e)) for e in els),'Selective dynamics','Cartesian']
    rows+=[' '.join('%.14f'%v for v in p)+' '+' '.join('T' if a else 'F' for a in m) for p,m in zip(s['xyz'],s['mask'])]
    return '\n'.join(rows)+'\n'
def append(s,e,p,label):
    need(label not in s['labels'],'Duplicate label')
    s['species'].append(e);s['xyz'].append(p);s['labels'].append(label);s['mask'].append([True]*3)

def graph():
    g={n:[] for n in PLAN['reference_models']}
    for side in ('M0','M1'):
        g[side+'_clean']=['Mg_bulk']
        for kind in PLAN['per_surface_models'][1:]:
            parent='clean' if kind in ('water','OH') else 'water' if kind=='OH_H' else 'OH_water' if kind=='2OH_H' else 'OH'
            g[side+'_'+kind]=[side+'_'+parent]
    need(len(g)==PLAN['maximum_models']==23,'Bounded model roster')
    return g
DAG=graph()

def reference(name):
    if name=='Mg_bulk':
        a,c=PLAN['initial_hcp_a_A'],PLAN['initial_hcp_c_A']
        cell=[[a,0,0],[-a/2,a*math.sqrt(3)/2,0],[0,0,c]]
        return {'cell':cell,'species':['Mg','Mg'],'xyz':[cart([0,0,0],cell),cart([2/3,1/3,.5],cell)],'mask':[[True]*3]*2,'labels':['b0','b1']}
    s={'cell':[[20.,0,0],[0,20.,0],[0,0,20.]],'species':[],'xyz':[],'mask':[],'labels':[]}
    if name=='H2':
        append(s,'H',[10,10,9.625],'gH0');append(s,'H',[10,10,10.375],'gH1')
    elif name=='H2O':
        append(s,'O',[10,10,10],'gO')
        half=math.radians(104.5/2);r=.97
        for sign,lab in ((1,'gH0'),(-1,'gH1')):append(s,'H',[10+sign*r*math.sin(half),10,10+r*math.cos(half)],lab)
    else:raise ValueError(name)
    return s

def bulk_parameters(s):
    c=s['cell'];a=(norm(c[0])+norm(c[1]))/2;h=norm(c[2])
    angle=math.degrees(math.acos(max(-1,min(1,dot(c[0],c[1])/norm(c[0])/norm(c[1])))))
    need(2.9<a<3.6 and 4.7<h<5.8,'Mg bulk lattice out of screening bounds')
    need(abs(norm(c[0])-norm(c[1]))<.003*a and abs(angle-120)<.2,'Bulk no longer hexagonal')
    need(abs(dot(c[0],c[2]))<1e-4 and abs(dot(c[1],c[2]))<1e-4,'Bulk axes not basal/normal')
    return a,h

def slab(side,bulk):
    a,c=bulk_parameters(bulk);u=[a,0,0];v=[-a/2,a*math.sqrt(3)/2,0];zvac=20.
    s={'cell':[mul(u,4),mul(v,4),[0,0,2.5*c+zvac]],'species':[],'xyz':[],'mask':[],'labels':[]}
    for layer in range(6):
        off=(0.,0.) if layer%2==0 else (2/3,1/3)
        for i,j in itertools.product(range(4),repeat=2):
            lab='m_%d_%d_%d'%(layer,i,j)
            p=plus(mul(u,i+off[0]),mul(v,j+off[1]));p[2]=zvac/2+layer*c/2
            append(s,'Ca' if side=='M1' and lab=='m_5_1_1' else 'Mg',p,lab)
            s['mask'][-1]=[layer>=2]*3
    return grouped(s)

def top(s,i,j):return atom(s,'m_5_%d_%d'%(i%4,j%4))
def bridge(s,i,j):
    p=top(s,i,j);q=top(s,i+1,j)
    return plus(p,mul(mic(q,p,s['cell']),.5))
def empty_h_site(s,old_label,o_label):
    options=[]
    for i,j in itertools.product(range(4),repeat=2):
        a=top(s,i,j);b=plus(a,mic(top(s,i+1,j),a,s['cell']));c=plus(a,mic(top(s,i+1,j+1),a,s['cell']))
        p=[(a[k]+b[k]+c[k])/3 for k in range(3)];p[2]+=1.10
        valid=True
        for q,e,l in zip(s['xyz'],s['species'],s['labels']):
            if l==old_label:continue
            limit=1.25 if e=='H' else 1.45 if e=='O' else 1.30
            if dist(p,q,s['cell'])<limit:valid=False;break
        length=dist(p,atom(s,o_label),s['cell'])
        if valid and 1.5<length<4.5:options.append((length,i,j,p))
    need(options,'No safe nearby H adsorption candidate; no forced dissociation')
    return min(options,key=lambda v:v[:3])[-1]
def water_ok(s):
    return all(.70<dist(atom(s,'wO'),atom(s,h),s['cell'])<1.25 for h in ('wH0','wH1'))
def OH_ok(s):return .70<dist(atom(s,'qO0'),atom(s,'qH0'),s['cell'])<1.25

def build(name,parent=None):
    if name in PLAN['reference_models']:return reference(name)
    side,kind=name.split('_',1)
    if kind=='clean':return slab(side,parent)
    s=copy.deepcopy(parent)
    if kind=='water':
        p=plus(top(s,1,1),[0,0,2.30]);append(s,'O',p,'wO')
        r=.97;t=math.radians(104.5/2)
        for sign,label in ((1,'wH0'),(-1,'wH1')):append(s,'H',plus(p,[sign*r*math.sin(t),0,r*math.cos(t)]),label)
    elif kind=='OH':
        p=plus(bridge(s,1,1),[0,0,1.65]);append(s,'O',p,'qO0');append(s,'H',plus(p,[0,0,.98]),'qH0')
    elif kind=='OH_water':
        need(OH_ok(s),'Parent hydroxyl dissociated; do not relabel as OH-covered')
        o=atom(s,'qO0');d=[0,2.65,.60];p=plus(o,d);u=mul(d,-1/norm(d))
        z=[0,0,1];perp=[z[k]-dot(z,u)*u[k] for k in range(3)];perp=mul(perp,1/norm(perp))
        t=math.radians(104.5);second=plus(mul(u,math.cos(t)),mul(perp,math.sin(t)))
        append(s,'O',p,'wO');append(s,'H',plus(p,mul(u,.97)),'wH0');append(s,'H',plus(p,mul(second,.97)),'wH1')
    elif kind in ('OH_H','2OH_H'):
        need(water_ok(s),'Selected molecular parent did not retain intact water; inspect actual product first')
        if kind=='2OH_H':need(OH_ok(s),'Original OH did not survive in second-water parent')
        s['xyz'][s['labels'].index('wH1')]=empty_h_site(s,'wH1','wO')
    elif kind in PLAN['coverage_patterns']:
        need(OH_ok(s),'Cannot extend a non-OH parent as OH coverage')
        p0=bridge(s,1,1);shift=mic(atom(s,'qO0'),p0,s['cell']);oh=mic(atom(s,'qH0'),atom(s,'qO0'),s['cell'])
        for n,(i,j) in enumerate(PLAN['coverage_patterns'][kind][1:],1):
            p=plus(bridge(s,i,j),shift);append(s,'O',p,'qO'+str(n));append(s,'H',plus(p,oh),'qH'+str(n))
    else:raise ValueError(kind)
    return grouped(s)

def geometry_check(s,name):
    n=len(s['species']);need(n==len(s['xyz'])==len(s['labels'])==len(s['mask']) and len(set(s['labels']))==n,'Atom identities')
    minima={}
    for i in range(n):
        for j in range(i):
            pair='-'.join(sorted([s['species'][i],s['species'][j]]));r=dist(s['xyz'][i],s['xyz'][j],s['cell'])
            minima[pair]=min(minima.get(pair,1e9),r)
            e={s['species'][i],s['species'][j]}
            lim=.65 if e=={'H'} else .70 if e=={'O','H'} else 1.60 if e=={'O'} else 1.10 if 'H' in e else 1.30 if 'O' in e else 2.30
            need(r>lim,'Collision '+name+' '+pair+' %.4f'%r)
    if name.startswith('M') and name!='Mg_bulk':
        side,kind=name.split('_',1)
        need(len(labels_metal(s))==96 and sum(not any(m) for m in s['mask'])==32,'Slab metal/free counts')
        need(s['species'].count('Ca')==(1 if side=='M1' else 0),'Surface Ca count')
    return minima

def conditions(name,phase):
    gas=name in ('H2','H2O');bulk=name=='Mg_bulk';relax=phase=='relax'
    nums=PLAN['numerics'];s={'SYSTEM':name+' '+phase,'PREC':'Accurate','ENCUT':520,'ISPIN':2,'LDAU':'.FALSE.','ISMEAR':0,'SIGMA':.05,'ALGO':'Normal','ISYM':2 if bulk else 0,'LREAL':'.FALSE.','LASPH':'.TRUE.','ADDGRID':'.TRUE.','EDIFF':'1E-8' if gas else '1E-7','NELM':200,'NCORE':1 if gas else 4,'ISTART':0,'ICHARG':2,'ISIF':3 if bulk and relax else 2,'IBRION':2 if relax else -1,'NSW':(120 if bulk or gas else 300) if relax else 0,'POTIM':.10,'EDIFFG':-.005 if gas or bulk else -.015,'LWAVE':'.TRUE.' if relax else '.FALSE.','LCHARG':'.FALSE.'}
    if not bulk:s.update({'LDIPOL':'.TRUE.','IDIPOL':4 if gas else 3,'DIPOL':'0.5 0.5 0.5'})
    return s

def write_model(name,s,preview=False):
    d=ROOT/name;need(not d.exists(),'Existing model directory: '+name);d.mkdir()
    minima=geometry_check(s,name);els=list(dict.fromkeys(s['species']))
    for phase in ('relax','static'):
        p=d/phase;p.mkdir();v=conditions(name,phase)
        v['MAGMOM']=' '.join(str(.2 if e=='O' and name not in ('H2O',) else 0.) for e in s['species'])
        (p/'INCAR').write_text('\n'.join(k+' = '+str(x) for k,x in v.items())+'\n')
        key='molecule_kmesh' if name in ('H2','H2O') else 'bulk_kmesh' if name=='Mg_bulk' else 'slab_kmesh'
        (p/'KPOINTS').write_text('Fresh common settings\n0\nGamma\n'+' '.join(map(str,PLAN['numerics'][key]))+'\n0 0 0\n')
        if not preview:(p/'POTCAR').write_bytes(b''.join((ROOT/'potentials'/e/'POTCAR').read_bytes() for e in els))
    (d/'relax/POSCAR').write_text(poscar(s,name+' | rebuilt; unrelaxed starting structure'))
    save(d/'model.json',{'name':name,'labels':s['labels'],'composition':{e:s['species'].count(e) for e in els},'parent':DAG[name],'preview':preview,'initial_minimum_distances_A':minima,'source':'new ideal hcp lattice and analytical adsorbates; no old geometry','kind':'local endpoint, not a passivating film or activation barrier'})

def local_potentials():
    """Read only configured licensed libraries. No old calculation geometry access."""
    labels=PLAN['chemistry']['potential_labels'];roots=[]
    for key in ('CA_PBE_POTENTIAL_ROOT','VASP_PP_PATH','VASP_PSP_DIR','PBE_PATH'):
        if os.environ.get(key):roots.append(Path(os.path.expandvars(os.path.expanduser(os.environ[key]))))
    cfg=HOME/'.vaspkit'
    if cfg.is_file():
        for line in cfg.read_text(errors='replace').splitlines():
            q=re.split('[#!]',line)[0].strip()
            if re.match(r'^PBE_PATH\s*=',q):
                pieces=shlex.split(q.split('=',1)[1]);
                if pieces:roots.append(Path(os.path.expandvars(os.path.expanduser(pieces[0]))))
    roots.append(HOME/'VASP/potcar_profiles')
    variants=['','PBE','potpaw_PBE','potpaw_PBE.52','potpaw_PBE.54','potpaw_PBE.64','potpaw_PBE_54']
    found={};info={}
    for e,label in labels.items():
        for root in roots:
            for variant in variants:
                for folder in (label,e):
                    p=root/variant/folder/'POTCAR'
                    if not p.is_file() or p.stat().st_size>4000000:continue
                    raw=p.read_bytes();titles=re.findall(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',raw)
                    if titles!=[label.encode()] or raw.count(b'End of Dataset')!=1:continue
                    enmax=re.findall(rb'ENMAX\s*=\s*([\d.]+)',raw)
                    need(enmax and 520+1e-5>=1.3*float(enmax[0]),'ENCUT inadequate for '+label+'; revise whole batch before submission')
                    found[e]=raw;info[e]={'label':label,'sha256':hashlib.sha256(raw).hexdigest(),'ENMAX_eV':float(enmax[0])};break
                if e in found:break
            if e in found:break
    need(set(found)==set(labels),'Missing configured licensed PBE potentials: '+','.join(sorted(set(labels)-set(found)))+'; no Slurm jobs submitted. Set CA_PBE_POTENTIAL_ROOT to the existing licensed PBE directory.')
    return found,info

def output(p):
    F=E0=None;ionic=end=scf=False;last_it=None;reading=False;buf=[];forces=[];stress=None;steps=0
    with (p/'OUTCAR').open(errors='replace') as f:
        for line in f:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and m.group(1)!=last_it:last_it=m.group(1);scf=False
            if 'aborting loop because EDIFF is reached' in line:scf=True
            if 'reached required accuracy' in line:ionic=True
            if 'General timing and accounting informations' in line:end=True
            if 'free  energy   TOTEN' in line:F=float(line.split('=')[1].split()[0]);steps+=1
            if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if re.match(r'^\s*in kB\s+',line):stress=list(map(float,line.split()[2:8]))
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buf=[];continue
            if reading:
                vals=line.split()
                if len(vals)==6:
                    try:buf.append(list(map(float,vals[3:])));continue
                    except ValueError:pass
                if buf:forces=buf;reading=False
    return {'F_eV':F,'E0_eV':E0,'normal_end':end,'last_scf_ediff':scf,'ionic_converged':ionic,'ionic_steps':steps,'forces':forces,'stress_kbar':stress}

def completed(name):
    d=ROOT/name;need((d/'VERIFIED_FINISHED').is_file(),'Parent not verified: '+name)
    meta=json.loads((d/'model.json').read_text());v=json.loads((d/'static/validation.json').read_text());need(v['verified'],'Parent final validation false')
    s=parse((d/'static/CONTCAR').read_text());s['labels']=meta['labels'];return s

def run(name):
    need(name in DAG,'Unapproved model')
    resource='gas' if name in ('H2','H2O') else 'bulk' if name=='Mg_bulk' else 'slab'
    mpi=PLAN['resources'][resource]['mpi'];need(int(os.environ.get('SLURM_NTASKS',0))==mpi,'Run only in prescribed Slurm allocation')
    root_meta=json.loads((ROOT/'manifest.json').read_text());need(not root_meta['preview'],'Preview cannot be executed as production')
    for e,v in root_meta['potentials'].items():need(sha(ROOT/'potentials'/e/'POTCAR')==v['sha256'],'Local potential changed')
    d=ROOT/name
    if (d/'VERIFIED_FINISHED').exists():print('ALREADY_COMPLETE|'+name);return
    if not d.exists():
        parent=completed(DAG[name][0]) if DAG[name] else None
        write_model(name,build(name,parent))
    (d/'.execution_started').mkdir(exist_ok=False)
    initial=parse((d/'relax/POSCAR').read_text());meta=json.loads((d/'model.json').read_text());initial['labels']=meta['labels']
    for phase in ('relax','static'):
        p=d/phase
        if phase=='static':shutil.copy2(d/'relax/CONTCAR',p/'POSCAR')
        with (p/'vasp.stdout').open('w') as out,(p/'vasp.stderr').open('w') as err:
            z=subprocess.run(['mpirun','-np',str(mpi),'vasp_std'],cwd=p,stdout=out,stderr=err)
        need(z.returncode==0,'VASP nonzero exit '+name+'/'+phase)
        v=output(p);s=parse((p/'CONTCAR').read_text());s['labels']=meta['labels']
        need(s['species']==initial['species'] and s['mask']==initial['mask'] and len(v['forces'])==len(s['species']),'Composition/constraints/forces mismatch')
        if name!='Mg_bulk':
            need(max(abs(x-y) for r,t in zip(s['cell'],initial['cell']) for x,y in zip(r,t))<1e-7,'Fixed cell changed')
            for a,b,m in zip(s['xyz'],initial['xyz'],s['mask']):
                if not any(m):need(dist(a,b,s['cell'])<1e-6,'Fixed atom moved')
        else:
            bulk_parameters(s);need(.8<volume(s['cell'])/volume(initial['cell'])<1.2,'Bulk volume outside initial tolerance')
        maxf=max(math.sqrt(sum(x*x for x,f in zip(vec,m) if f)) for vec,m in zip(v.pop('forces'),s['mask']))
        lim=.005 if name in PLAN['reference_models'] else .015 if phase=='relax' else .02
        valid=v['normal_end'] and v['last_scf_ediff'] and (phase=='static' or v['ionic_converged']) and maxf<lim and all(v[k] is not None and math.isfinite(v[k]) for k in ('F_eV','E0_eV'))
        if name=='Mg_bulk':valid=valid and v['stress_kbar'] is not None and max(abs(x) for x in v['stress_kbar'])<2.0
        v.update({'max_free_force_eV_A':maxf,'force_gate':lim,'verified':bool(valid),'volume_A3':volume(s['cell'])})
        os_=[i for i,e in enumerate(s['species']) if e=='O'];hs=[i for i,e in enumerate(s['species']) if e=='H']
        v['O_H_coordination']=[{'O_label':s['labels'][i],'H_bonds':[s['labels'][j] for j in hs if dist(s['xyz'][i],s['xyz'][j],s['cell'])<1.25]} for i in os_]
        v['minimum_HH_A']=min((dist(s['xyz'][i],s['xyz'][j],s['cell']) for i,j in itertools.combinations(hs,2)),default=None)
        save(p/'validation.json',v);print('VALIDATION|'+name+'|'+phase+'|'+json.dumps(v),flush=True)
        need(valid,'Convergence/force/stress gate failed: '+name+'/'+phase)
    (d/'VERIFIED_FINISHED').touch()

def shell(name):
    group='gas' if name in ('H2','H2O') else 'bulk' if name=='Mg_bulk' else 'slab';r=PLAN['resources'][group]
    return ('#!/bin/bash\n#SBATCH --partition=64c512g\n#SBATCH --nodes=1\n#SBATCH --ntasks='+str(r['mpi'])+'\n#SBATCH --mem='+str(r['memory_GB'])+'G\n#SBATCH --time='+r['walltime']+'\n#SBATCH --output=logs/%x-%j.out\n#SBATCH --error=logs/%x-%j.err\n'
            'set -euo pipefail\ncd "$SLURM_SUBMIT_DIR"\nmodule purge\nmodule load oneapi vasp/6.3.2-intel-2021.4.0\nexport OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\nulimit -s unlimited\nulimit -l unlimited\npython3 fresh.py --run '+name+'\n')

def setup(preview=False):
    need(not ROOT.exists(),'Destination exists; do not overwrite')
    potentials,info=({}, {}) if preview else local_potentials()
    ROOT.mkdir(parents=True);(ROOT/'logs').mkdir();(ROOT/'jobs').mkdir()
    shutil.copy2(HERE/'plan.json',ROOT/'plan.json');shutil.copy2(Path(__file__),ROOT/'fresh.py')
    if not preview:
        for e,raw in potentials.items():
            p=ROOT/'potentials'/e;p.mkdir(parents=True);(p/'POTCAR').write_bytes(raw)
    seeds={}
    for name in DAG:
        if preview or name in PLAN['reference_models']:
            parent=seeds[DAG[name][0]] if DAG[name] else None;s=build(name,parent);seeds[name]=s
            write_model(name,s,preview)
        f=ROOT/'jobs'/(name+'.slurm');f.write_text(shell(name));subprocess.run(['bash','-n',str(f)],check=True)
    save(ROOT/'manifest.json',{'created_at':timestamp(),'preview':preview,'potentials':info,'DAG':DAG,'maximum_models':23,'generated_initial_models':list(seeds),'source':'Fresh analytical geometry only','code_sha256':sha(ROOT/'fresh.py'),'plan_sha256':sha(ROOT/'plan.json')})
    print('PREPARED|'+str(ROOT)+'|initial_models='+str(len(seeds))+'|preview='+str(preview),flush=True)

def submit():
    need(shutil.which('sbatch') is not None,'Slurm not available here')
    if not ROOT.exists():setup(False)
    m=json.loads((ROOT/'manifest.json').read_text());need(not m['preview'],'Cannot submit preview')
    need(sha(ROOT/'fresh.py')==m['code_sha256'] and sha(ROOT/'plan.json')==m['plan_sha256'],'Production driver/plan changed')
    ledger=ROOT/'submission.json'
    if ledger.exists():
        print('EXISTING_SUBMISSION|'+ledger.read_text());return
    (ROOT/'.submission_attempted').mkdir(exist_ok=False);record={'created_at':timestamp(),'jobs':{},'complete':False,'max_models':23}
    save(ledger,record)
    for name,parents in DAG.items():
        cmd=['sbatch','--parsable','--job-name=CaFresh_'+name]
        if parents:cmd+=['--dependency=afterok:'+':'.join(record['jobs'][p] for p in parents),'--kill-on-invalid-dep=yes']
        cmd+=['jobs/'+name+'.slurm']
        p=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=45)
        need(p.returncode==0,'Submission rejected '+name+': '+p.stderr[:300])
        j=p.stdout.strip().split(';')[0];need(re.fullmatch(r'\d+',j),'Uncertain sbatch response; inspect ledger, never blindly retry')
        record['jobs'][name]=j;save(ledger,record);print('SUBMITTED|'+name+'|'+j,flush=True)
    record['complete']=True;record['submitted_at']=timestamp();save(ledger,record)

def status():
    print('CHECKED_AT|'+timestamp());p=ROOT/'submission.json'
    if not p.exists():print('NOT_SUBMITTED');return
    d=json.loads(p.read_text());print('SUBMISSION|'+json.dumps(d))
    if d['jobs']:
        ids=','.join(d['jobs'].values())
        for cmd in (['squeue','-h','-j',ids,'-o','%i|%j|%T|%M|%R'],['sacct','-X','-n','-P','-j',ids,'-o','JobID,JobName%35,State,ExitCode,Elapsed,Start,End']):
            q=subprocess.run(cmd,capture_output=True,text=True,timeout=30);print('SLURM_RC|'+str(q.returncode));print(q.stdout)
    for name in DAG:
        q=ROOT/name/'static/validation.json';print('MODEL|'+name+'|'+(q.read_text() if q.exists() else 'NO_FINAL_VALIDATION'))

def main():
    global ROOT
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--preview',action='store_true');g.add_argument('--submit',action='store_true');g.add_argument('--status',action='store_true');g.add_argument('--run',choices=list(DAG));p.add_argument('--root',type=Path)
    a=p.parse_args()
    if a.root:need(a.preview,'Custom root only permitted for preview');ROOT=a.root.resolve()
    if a.preview:setup(True)
    elif a.submit:submit()
    elif a.status:status()
    else:run(a.run)
if __name__=='__main__':main()
