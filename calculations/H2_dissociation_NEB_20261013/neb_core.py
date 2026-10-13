"""Geometry, VASP input/output, and IDPP helpers for two fixed-cell H2 paths.
This module contains no VASP source or pseudopotential data.
"""
import copy, hashlib, itertools, json, math, re
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
FLOAT=r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?'
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.tmp')
    def enc(x):
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,np.generic): return x.item()
        raise TypeError(type(x).__name__)
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,default=enc));tmp.replace(path)
def parse(raw):
    ls=raw.splitlines();require(len(ls)>8,'Incomplete POSCAR')
    scale=float(ls[1]);require(scale>0,'Only positive single scale supported')
    cell=np.array([list(map(float,l.split()[:3])) for l in ls[2:5]])*scale
    require(np.max(np.abs(cell-np.diag(np.diag(cell))))<1e-7,'This screened path requires the audited orthogonal cells')
    sp=sum(([e]*int(n) for e,n in zip(ls[5].split(),ls[6].split())),[])
    i=7;sd=ls[i].strip().lower().startswith('s');i+=int(sd);direct=ls[i].strip().lower().startswith('d');i+=1
    rows=[l.split() for l in ls[i:i+len(sp)]];require(len(rows)==len(sp),'Atom count mismatch')
    xyz=np.array([list(map(float,r[:3])) for r in rows]);xyz=xyz@cell if direct else xyz*scale
    mask=np.array([[v.upper().startswith('T') for v in r[3:6]] if sd else [True]*3 for r in rows],bool)
    require(mask.shape==xyz.shape and np.isfinite(xyz).all(),'Bad mask/coordinates')
    return {'cell':cell,'sp':sp,'xyz':xyz,'mask':mask}
def clone(s): return {k:v.copy() if isinstance(v,(list,np.ndarray)) else copy.deepcopy(v) for k,v in s.items()}
def poscar(s,title):
    els=list(dict.fromkeys(s['sp']));require(sum(([e]*s['sp'].count(e) for e in els),[])==s['sp'],'Species not grouped')
    ls=[title,'1.0']+[' '.join('%.14f'%v for v in row) for row in s['cell']]+[' '.join(els),' '.join(str(s['sp'].count(e)) for e in els),'Selective dynamics','Cartesian']
    ls+=[' '.join('%.14f'%v for v in row)+' '+' '.join('T' if v else 'F' for v in mask) for row,mask in zip(s['xyz'],s['mask'])]
    return '\n'.join(ls)+'\n'
def mic(delta,cell):
    lengths=np.diag(cell);return delta-np.rint(delta/lengths)*lengths
def distance(a,b,cell):return float(np.linalg.norm(mic(np.asarray(a)-b,cell)))
def summary_geometry(s):
    x=s['xyz'];cell=s['cell'];hs=[i for i,e in enumerate(s['sp']) if e=='H'];met=[i for i,e in enumerate(s['sp']) if e!='H']
    out={'composition':{e:s['sp'].count(e) for e in dict.fromkeys(s['sp'])},'n_atoms':len(x),'n_fixed_atoms':int((~s['mask'].any(axis=1)).sum())}
    if len(hs)==2:out['H_H_A']=distance(x[hs[0]],x[hs[1]],cell)
    out['H_neighbours']=[{'index_1based':i+1,'nearest_metals':sorted([(distance(x[i],x[j],cell),s['sp'][j],j+1) for j in met])[:4]} for i in hs]
    if met:
        out['slab_z_span_A']=float(np.ptp(x[met,2]));out['nominal_vacuum_A']=float(cell[2,2]-np.ptp(x[met,2]))
    return out

def incar(system,phase,nbands,gas=False,neb=None):
    # NBANDS is explicit and identical for each solid's references/endpoints/images.
    d={'SYSTEM':system+' '+phase,'PREC':'Accurate','ENCUT':520,'ISPIN':2,'LDAU':'.FALSE.','ISMEAR':0,'SIGMA':.05,'ALGO':'Normal','ISYM':0,'ISIF':2,'LREAL':'.FALSE.','LASPH':'.TRUE.','ADDGRID':'.TRUE.','EDIFF':'1E-8' if gas else '1E-7','NELM':200,'NELMIN':4,'NCORE':1 if gas else 4,'KPAR':1,'NBANDS':nbands,'ISTART':0,'ICHARG':2,'LWAVE':'.TRUE.','LCHARG':'.FALSE.','LORBIT':11,'IBRION':2 if phase=='relax' else -1,'NSW':(100 if gas else 300) if phase=='relax' else 0,'POTIM':.10,'EDIFFG':-.005 if gas else -.01}
    if gas:d['NUPDOWN']=0
    else:d.update({'LDIPOL':'.TRUE.','IDIPOL':3,'DIPOL':'0.5 0.5 0.5'})
    if neb is not None:
        climb=neb=='ci';d.update({'IBRION':3,'POTIM':0,'IOPT':7,'ICHAIN':0,'IMAGES':5,'SPRING':-5,'LCLIMB':'.TRUE.' if climb else '.FALSE.','LTANGENTOLD':'.FALSE.','NSW':400 if climb else 250,'EDIFFG':-.03 if climb else -.05,'TIMESTEP':.10,'FTIMEMAX':.50,'MAXMOVE':.10})
    return d
def input_text(d,sp):
    d=dict(d);d['MAGMOM']=' '.join('1.0' if e=='Ni' else '0.0' for e in sp)
    return '\n'.join(str(k)+' = '+str(v) for k,v in d.items())+'\n'

def output(folder,mask):
    p=Path(folder)/'OUTCAR';v={'F_eV':None,'E0_eV':None,'last_ediff':False,'normal_end':False,'ionic_converged':False,'ionic_steps':0,'vtst_force_max':None,'lclimb':None,'neb_marker':False}
    forces=[];buf=[];reading=False;it=None
    with p.open(errors='replace') as f:
        for line in f:
            m=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
            if m and it!=m.group(1):it=m.group(1);v['last_ediff']=False
            if 'aborting loop because EDIFF is reached' in line:v['last_ediff']=True
            if 'General timing and accounting informations' in line:v['normal_end']=True
            if 'reached required accuracy' in line:v['ionic_converged']=True
            if 'free  energy   TOTEN' in line:v['F_eV']=float(line.split('=')[1].split()[0]);v['ionic_steps']+=1
            if 'energy(sigma->0)' in line:v['E0_eV']=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
            if 'FORCES: max atom, RMS' in line:
                vals=re.findall(FLOAT,line.split('RMS',1)[1]);v['vtst_force_max']=float(vals[0].replace('D','E')) if vals else None
            if 'NEB:' in line:v['neb_marker']=True
            if 'LCLIMB' in line.upper():
                m=re.search(r'LCLIMB\s*[=:]?\s*\.?([TF])',line.upper())
                if m:v['lclimb']=m.group(1)=='T'
            if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;buf=[];continue
            if reading:
                q=line.split()
                if len(q)==6:
                    try:buf.append(list(map(float,q[3:6])));continue
                    except ValueError:pass
                if buf:forces=buf;reading=False
    if reading and buf:forces=buf
    require(np.asarray(forces).shape==np.asarray(mask).shape,'Missing force block: '+str(folder))
    v['max_free_DFT_force']=float(np.linalg.norm(np.array(forces)*mask,axis=1).max())
    return v

def endpoint_check(s,kind):
    g=summary_geometry(s);hh=g.get('H_H_A')
    if kind=='molecule':require(hh is not None and .60<hh<1.05,'Molecular endpoint lost H2 identity; do not connect blindly')
    if kind=='product':require(hh is not None and 1.50<hh<6.0,'Dissociated endpoint recombined or exceeded selected local path')
    if kind in ('molecule','product'):
        for h in g['H_neighbours']:
            ds=h['nearest_metals'];require(ds and ds[0][0]>(.80 if kind=='product' else .70),'H-metal collision')
            require(ds[0][0]<(2.8 if kind=='product' else 5.5),'H has escaped selected adsorption region')
    return g

def initial_band(initial,final,nimages=5):
    """Conserve indices; choose only the two identical H labels; unwrap periodic images.
    Use a rotating/stretching H-H vector, then bounded IDPP relaxation. No DFT here.
    """
    require(initial['sp']==final['sp'] and np.allclose(initial['cell'],final['cell'],atol=1e-9),'Endpoint composition/cell mismatch')
    require(np.array_equal(initial['mask'],final['mask']),'Endpoint constraint mismatch')
    s=clone(initial);e=clone(final);cell=s['cell'];mask=s['mask'];require(s['sp'][-2:]==['H','H'],'H must be the final two atoms')
    require(np.max(np.abs((e['xyz']-s['xyz'])[~mask]),initial=0)<1e-8,'Fixed atom positions differ')
    s['xyz'][-1]=s['xyz'][-2]+mic(s['xyz'][-1]-s['xyz'][-2],cell)
    e['xyz'][:-2]=s['xyz'][:-2]+mic(e['xyz'][:-2]-s['xyz'][:-2],cell)
    hc0=s['xyz'][-2:].mean(axis=0);candidates=[]
    for swapped in (False,True):
        pair=e['xyz'][-2:].copy()
        if swapped:pair=pair[::-1]
        pair[1]=pair[0]+mic(pair[1]-pair[0],cell);hc=pair.mean(axis=0);pair+=hc0+mic(hc-hc0,cell)-hc
        candidates.append((float(((pair-s['xyz'][-2:])**2).sum()),swapped,pair))
    _,swapped,e['xyz'][-2:]=min(candidates,key=lambda a:a[0]);v0=s['xyz'][-1]-s['xyz'][-2];v1=e['xyz'][-1]-e['xyz'][-2]
    r0=np.linalg.norm(v0);r1=np.linalg.norm(v1);u0=v0/r0;u1=v1/r1
    require(np.dot(u0,u1)>-.95,'Antiparallel H mapping; path needs manual design')
    require(np.linalg.norm(e['xyz'][-2:].mean(axis=0)-hc0)<4.5,'Endpoints too spatially separated for one local dissociation path')
    require(np.max(np.linalg.norm(e['xyz'][:-2]-s['xyz'][:-2],axis=1))<1.5,'Large substrate reconstruction; separate mechanism required')
    n=len(s['sp']);ii,jj=np.triu_indices(n,1)
    def ds(x):return np.linalg.norm(mic(x[ii]-x[jj],cell),axis=1)
    d0=ds(s['xyz']);d1=ds(e['xyz']);images=[s];reports=[]
    for idx in range(1,nimages+1):
        t=idx/(nimages+1);x=(1-t)*s['xyz']+t*e['xyz'];u=(1-t)*u0+t*u1;u/=np.linalg.norm(u);center=(1-t)*hc0+t*e['xyz'][-2:].mean(axis=0);length=(1-t)*r0+t*r1
        x[-2]=center-.5*length*u;x[-1]=center+.5*length*u;seed=x.copy();target=(1-t)*d0+t*d1
        def objective(z):
            a=seed.copy();a[mask]=z;vec=mic(a[ii]-a[jj],cell);r=np.linalg.norm(vec,axis=1);r=np.maximum(r,1e-6);diff=r-target
            value=float(np.sum(diff**2/r**4));fac=2*diff/r**5-4*diff**2/r**6;g=np.zeros_like(a);np.add.at(g,ii,fac[:,None]*vec);np.add.at(g,jj,-fac[:,None]*vec)
            return value,g[mask]
        limits=np.full_like(seed,.35);limits[-2:]=.80
        bounds=list(zip((seed-limits)[mask],(seed+limits)[mask]));res=minimize(objective,x[mask],jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':180,'ftol':1e-11,'gtol':1e-6})
        require(np.isfinite(res.x).all(),'IDPP nonfinite values');x[mask]=res.x
        require(objective(res.x)[0]<=objective(seed[mask])[0]+1e-9,'IDPP increased distance objective')
        cur=clone(s);cur['xyz']=x;g=summary_geometry(cur);require(g['H_H_A']>.55,'Interpolated H-H collision')
        for h in g['H_neighbours']:require(h['nearest_metals'][0][0]>.85,'Interpolated H-metal collision')
        mii,mjj=np.triu_indices(n-2,1);require(np.min(np.linalg.norm(mic(x[mii]-x[mjj],cell),axis=1))>1.45,'Interpolated metal collision')
        require(np.allclose(x[~mask],s['xyz'][~mask],atol=1e-12),'Fixed atoms moved')
        reports.append({'image':idx,'idpp_success':bool(res.success),'idpp_message':str(res.message),'objective':float(res.fun),'H_H_A':g['H_H_A'],'geometry':g})
        images.append(cur)
    images.append(e);return images,{'H_labels_swapped_at_FS':swapped,'image_reports':reports,'interpolation':'pair-rotation seed followed by bounded IDPP; preliminary geometries, not DFT trajectory'}
