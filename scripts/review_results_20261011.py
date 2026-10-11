#!/usr/bin/env python3
"""Read-only, bounded audit of the explicitly requested VASP calculations."""
import os, re, json, math, hashlib, itertools, subprocess
from pathlib import Path
from collections import deque
H=Path.home()
R=H/'Mg2Ni_Mg_H_thermo_validation_20261010'
T=H/'TaskA_Mg17La2_H_enrichment_20261010'
def text(p,limit=120000):
    try:
        if p.stat().st_size>limit:return '[OMITTED: too large]'
        return p.read_text(errors='replace')
    except OSError:return ''
def sha(p):
    if not p.is_file():return None
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def incar(p):
    d={}
    for line in text(p).splitlines():
        for s in re.split('[;]',re.split('[!#]',line)[0]):
            if '=' in s:
                k,v=s.split('=',1);d[k.strip().upper()]=v.strip()
    return d
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def structure(p):
    lines=text(p).splitlines()
    if len(lines)<8:return None
    scale=float(lines[1]);assert scale>0
    cell=[[float(v)*scale for v in l.split()[:3]] for l in lines[2:5]]
    elements=lines[5].split();counts=[int(v) for v in lines[6].split()];n=sum(counts)
    i=7;select=lines[i].strip().lower().startswith('s')
    if select:i+=1
    direct=lines[i].strip().lower().startswith('d');i+=1
    coords=[];masks=[]
    for line in lines[i:i+n]:
        q=line.split();v=[float(x) for x in q[:3]]
        coords.append([sum(v[k]*cell[k][j] for k in range(3)) for j in range(3)] if direct else [x*scale for x in v])
        masks.append([x.upper().startswith('T') for x in q[3:6]] if select else [True]*3)
    assert len(coords)==n and all(len(x)==3 for x in masks)
    return {'cell':cell,'elements':elements,'counts':counts,'coords':coords,'mask':masks,'species':sum(([e]*c for e,c in zip(elements,counts)),[])}
def distance(a,b,cell):
    delta=[x-y for x,y in zip(a,b)];det=dot(cell[0],cross(cell[1],cell[2]))
    frac=[dot(delta,cross(cell[1],cell[2]))/det,dot(delta,cross(cell[2],cell[0]))/det,dot(delta,cross(cell[0],cell[1]))/det]
    nearest=[round(x) for x in frac]
    return min(math.sqrt(sum((delta[j]-sum((nearest[k]+t[k])*cell[k][j] for k in range(3)))**2 for j in range(3))) for t in itertools.product((-1,0,1),repeat=3))
def h_neighbours(s):
    if not s:return []
    hs=[i for i,e in enumerate(s['species']) if e=='H'];out=[]
    for i in hs:
        if len(hs)>4 and i not in hs[-2:]:continue
        ns=sorted((distance(s['coords'][i],p,s['cell']),s['species'][j],j+1) for j,p in enumerate(s['coords']) if i!=j)
        out.append({'atom_1based':i+1,'xyz_A':s['coords'][i],'nearest6':ns[:6]})
    return out
def outcar(p,s):
    if not p.is_file():return None
    energies=[];e0=None;force_history=[];inforce=False;fv=[];footer=False;ionic=False;last_ediff=False;n_scf=0;scf_counts=[];last_iter=0
    with p.open(errors='replace') as f:
        for line in f:
            if 'Iteration' in line:
                m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
                if m:
                    it=int(m.group(1));last_iter=int(m.group(2))
                    if it!=n_scf:n_scf=it;last_ediff=False
            if 'aborting loop because EDIFF is reached' in line:last_ediff=True
            if 'free  energy   TOTEN' in line:
                m=re.search(r'TOTEN\s*=\s*([-+0-9.Ee]+)',line)
                if m:energies.append(float(m.group(1)));scf_counts.append(last_iter)
            if 'energy(sigma->0)' in line:
                m=re.search(r'energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)',line)
                if m:e0=float(m.group(1))
            if 'reached required accuracy' in line:ionic=True
            if 'General timing and accounting informations for this job' in line:footer=True
            if 'TOTAL-FORCE (eV/Angst)' in line:inforce=True;fv=[];continue
            if inforce:
                cols=line.split()
                if len(cols)==6:
                    try:fv.append([float(x) for x in cols[3:6]]);continue
                    except ValueError:pass
                if fv:
                    mask=s['mask'] if s and len(s['mask'])==len(fv) else [[True]*3 for v in fv]
                    allf=[math.sqrt(sum(x*x for x in v)) for v in fv]
                    freef=[math.sqrt(sum(x*x for x,m in zip(v,ms) if m)) for v,ms in zip(fv,mask)]
                    order=sorted(range(len(fv)),key=lambda j:freef[j],reverse=True)
                    force_history.append({'all':max(allf),'free':max(freef),'top':[(j+1,s['species'][j] if s else '?',freef[j]) for j in order[:3]]});inforce=False
    return {'F_eV':energies[-1] if energies else None,'E0_eV':e0,'ionic_steps':len(energies),'final_scf_ediff':last_ediff,'last_scf_iterations':scf_counts[-1:] ,'max_scf_iterations':max(scf_counts) if scf_counts else None,'ionic_converged':ionic,'footer':footer,'last_force':force_history[-1:] ,'force_last10':[round(v['free'],6) for v in force_history[-10:]],'F_last10_eV':energies[-10:]}
def model(root):
    d={'path':str(root.relative_to(H)),'metadata':{},'stages':{}}
    for p in sorted(root.iterdir()):
        if p.is_file() and p.suffix in ('.json','.md') and p.stat().st_size<50000:d['metadata'][p.name]=text(p,50000)
    for st in ('01_relax','02_static','01_static'):
        sd=root/st
        if not sd.is_dir():continue
        a=structure(sd/'POSCAR');b=structure(sd/'CONTCAR')
        d['stages'][st]={'incar':incar(sd/'INCAR'),'kpoints':text(sd/'KPOINTS'),'potcar_sha256':sha(sd/'POTCAR'),'input_composition':dict(zip(a['elements'],a['counts'])) if a else None,'cell':a['cell'] if a else None,'n_free_atoms':sum(any(x) for x in a['mask']) if a else None,'mask_sha':hashlib.sha256(json.dumps(a['mask']).encode()).hexdigest() if a else None,'input_sha256':sha(sd/'POSCAR'),'contcar_sha256':sha(sd/'CONTCAR'),'neighbours_in':h_neighbours(a),'neighbours_out':h_neighbours(b),'outcar':outcar(sd/'OUTCAR',a),'wavecar_bytes':(sd/'WAVECAR').stat().st_size if (sd/'WAVECAR').exists() else 0,'chgcar_bytes':(sd/'CHGCAR').stat().st_size if (sd/'CHGCAR').exists() else 0}
    return d
report={'models':{},'project_files':{},'job_logs':{}}
for root in (R,T):
    report['project_files'][root.name]={}
    if not root.is_dir():continue
    for p in sorted(root.iterdir()):
        if p.is_dir() and (re.match(r'^[ABCD](0|-)',p.name) or p.name.startswith('H3_add')):
            try:report['models'][p.name]=model(p)
            except Exception as e:report['models'][p.name]={'parse_error':str(e)}
        if p.is_file() and p.suffix in ('.sh','.slurm','.json','.md','.csv','.tsv') and p.stat().st_size<70000:
            report['project_files'][root.name][p.name]=text(p,70000)
    for base,dirs,files in os.walk(root):
        if Path(base).relative_to(root).parts and len(Path(base).relative_to(root).parts)>2:dirs[:]=[];continue
        for fn in files:
            if any(x in fn for x in ('63535597','63535598')) and not fn.endswith('.xml'):
                p=Path(base)/fn
                with p.open(errors='replace') as f:tail=''.join(deque(f,80))
                report['job_logs'][str(p.relative_to(H))]=tail
p=H/'Mg_LaH3_Hrelay_endpoints_20261008/S0_base_static_retry_20261010'
if p.is_dir():report['models']['Hrelay_S0']=model(p)
q=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%i|%j|%T|%M|%R'],capture_output=True,text=True,timeout=30)
report['queue']=q.stdout
print(json.dumps(report,ensure_ascii=False))
