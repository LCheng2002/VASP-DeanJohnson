#!/usr/bin/env python3
"""Run on the user's cluster. No network login, credentials, NEB or original edits.
--preview builds six unrelaxed seeds without POTCAR or Slurm.
--submit prepares licensed local POTCARs and submits the fixed dependency graph.
Children are constructed only from verified, converged parents.
"""
import argparse, copy, hashlib, itertools, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())
SHIFTS=list(itertools.product((-1,0,1),repeat=3))
ROOT=Path.home()/PLAN['batch']
SOURCE_ROOT=HERE.parents[1]
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def save(p,obj):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(obj,indent=2,ensure_ascii=False));t.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def norm(v):return math.sqrt(dot(v,v))
def frac(v,c):
 det=dot(c[0],cross(c[1],c[2]));require(abs(det)>1e-8,'Singular cell')
 return [dot(v,cross(c[1],c[2]))/det,dot(v,cross(c[2],c[0]))/det,dot(v,cross(c[0],c[1]))/det]
def cart(v,c):return [sum(v[k]*c[k][j] for k in range(3)) for j in range(3)]
def mic(a,b,c):
 f=frac([x-y for x,y in zip(a,b)],c);r=[round(x) for x in f]
 return min((cart([f[k]-r[k]-t[k] for k in range(3)],c) for t in SHIFTS),key=lambda v:dot(v,v))
def distance(a,b,c):return norm(mic(a,b,c))
def parse(raw):
 l=raw.splitlines();require(len(l)>8,'Incomplete POSCAR');scale=float(l[1]);require(scale>0,'Expected positive scale')
 c=[[float(v)*scale for v in x.split()[:3]] for x in l[2:5]];els=l[5].split();ct=list(map(int,l[6].split()));sp=sum(([e]*n for e,n in zip(els,ct)),[])
 j=7;sd=l[j].strip().lower().startswith('s');j+=int(sd);direct=l[j].strip().lower().startswith('d');j+=1
 xyz=[];mask=[]
 for row in l[j:j+len(sp)]:
  q=row.split();v=list(map(float,q[:3]));xyz.append(cart(v,c) if direct else [x*scale for x in v]);mask.append([v.upper().startswith('T') for v in q[3:6]] if sd else [True]*3)
 require(len(xyz)==len(sp) and all(len(a)==3 for a in mask),'Atom or mask count invalid')
 return {'cell':c,'species':sp,'xyz':xyz,'mask':mask,'ids':list(range(1,len(sp)+1))}
def poscar(s,name):
 els=list(dict.fromkeys(s['species']));require(sum(([e]*s['species'].count(e) for e in els),[])==s['species'],'Ungrouped species')
 lines=[name,'1.0']+[' '.join('%.14f'%v for v in a) for a in s['cell']]+[' '.join(els),' '.join(str(s['species'].count(e)) for e in els),'Selective dynamics','Cartesian']
 return '\n'.join(lines+[' '.join('%.14f'%v for v in a)+' '+' '.join('T' if v else 'F' for v in m) for a,m in zip(s['xyz'],s['mask'])])+'\n'
def subset(s,indices):return {**copy.deepcopy(s),**{k:[copy.deepcopy(s[k][i]) for i in indices] for k in ('species','xyz','mask','ids')}}
def neighbours(s,p,element=None,exclude=None):
 return sorted((distance(p,a,s['cell']),e,j) for a,e,j in zip(s['xyz'],s['species'],s['ids']) if j!=exclude and (element is None or element==e))
def centered(s):
 s=copy.deepcopy(s);fs=[frac(a,s['cell']) for a in s['xyz']];zs=sorted(f[2]%1 for f in fs);g=[zs[i+1]-zs[i] for i in range(len(zs)-1)]+[zs[0]+1-zs[-1]];k=max(range(len(g)),key=g.__getitem__);cut=zs[(k+1)%len(zs)]
 for f in fs:f[2]=(f[2]-cut)%1
 xyz=[cart(f,s['cell']) for f in fs];offset=(s['cell'][2][2]-max(a[2] for a in xyz))/2
 for a in xyz:a[2]+=offset
 s['xyz']=xyz;return s
def source(key):
 info=PLAN['sources'][key];p=SOURCE_ROOT/info['path'];b=p.read_bytes()
 require(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\x00'+b).hexdigest()==info['git_blob_sha1'],'Source changed: '+info['path'])
 return centered(parse(b.decode()))
def initial_models():
 l,c,a=source('L'),source('C'),source('A')
 require(l['species']==['Mg']*12+['La']*16+['H']*72,'La source composition')
 require(c['species']==['Mg']*30+['Ni']*6+['H']*36,'C source composition')
 out={'L_parent':(l,'L',{}),'L_MgH2_reference':(subset(l,[i for i,j in enumerate(l['ids']) if j<=12 or 29<=j<=52]),'L',{'finite_cleavage_control':True}),'Mg_parent':(subset(a,list(range(24))),'Mg',{'finite_cleavage_control':True})}
 c2={'cell':[c['cell'][0],[2*v for v in c['cell'][1]],c['cell'][2]],'species':[],'xyz':[],'ids':[],'mask':[]}
 for e in ('Mg','Ni','H'):
  for rep in range(2):
   for i,sp in enumerate(c['species']):
    if sp==e:
     c2['species'].append(e);c2['xyz'].append([c['xyz'][i][j]+rep*c['cell'][1][j] for j in range(3)]);c2['ids'].append(c['ids'][i]+72*rep);c2['mask'].append(c['mask'][i])
 pool=PLAN['low_H_proxy']['candidate_pool']
 for suffix,key in [('A','arrangement_A_pool_indices'),('B','arrangement_B_pool_indices')]:
  s=copy.deepcopy(c2);sites=[]
  for hidx,pidx in enumerate(PLAN['low_H_proxy'][key]):
   ni,*v=pool[pidx];p=[s['xyz'][s['ids'].index(ni)][j]+v[j] for j in range(3)]
   require(neighbours(s,p,'Mg')[0][0]>1.60 and neighbours(s,p,'Ni')[0][0]>1.45 and neighbours(s,p,'H')[0][0]>1.75,'Low-H initial collision')
   s['species'].append('H');s['xyz'].append(p);s['ids'].append(145+hidx);s['mask'].append([True]*3);sites.append({'anchor':ni,'vector':v})
  require(len(s['ids'])==147 and s['species'].count('H')==75,'Low-H count')
  out['X025_'+suffix]=(s,'X',{'x_initial':.25,'literature_x':.3,'proxy_not_exact':True,'site_choices':sites,'background_H_ids':[145,146,147],'transferred_H_ids':[]})
 h={'cell':[[20.,0.,0.],[0.,20.,0.],[0.,0.,20.]],'species':['H','H'],'xyz':[[10.,10.,9.625],[10.,10.,10.375]],'mask':[[True]*3]*2,'ids':[1,2]}
 out['H2_reference']=(h,'H2',{})
 return out

def incar(s,name,stage):
 gas=name=='H2_reference';vals={'SYSTEM':name+' '+stage,'PREC':'Accurate','ENCUT':520,'ISPIN':2,'MAGMOM':' '.join(str(1.0 if e=='Ni' else .1 if e=='La' else 0.0) for e in s['species']),'LDAU':'.FALSE.','ISMEAR':0,'SIGMA':.05,'ALGO':'Normal','ISYM':0,'ISIF':2,'LREAL':'.FALSE.','LASPH':'.TRUE.','ADDGRID':'.TRUE.','EDIFF':'1E-8' if gas else '1E-7','NELM':200,'NCORE':1 if gas else 8,'ISTART':0,'ICHARG':2,'LWAVE':'.FALSE.','LCHARG':'.FALSE.','IBRION':2 if stage=='01_relax' else -1,'NSW':(100 if gas else 400) if stage=='01_relax' else 0,'POTIM':.1,'EDIFFG':-.005 if gas else -.015}
 # Neutral two-surface endpoint models. Dipole and slab-size sensitivity remain explicit limitations.
 return '\n'.join(k+' = '+str(v) for k,v in vals.items())+'\n'
def local_potentials():
 base=Path.home();p=base/'Mg2Ni_Mg_H_thermo_validation_20261010/A-H-Ni_Mg_Mg2Ni/01_relax/POTCAR'
 require(sha(p)=='72b607bd5434dcb71996bb13fd478af3592115762200d3c53e0f47c68b0a84ae','Canonical Mg/Ni/H POTCAR changed')
 def split(path):
  blocks=[];buf=[]
  for line in path.read_bytes().splitlines(keepends=True):
   buf.append(line)
   if line.strip()==b'End of Dataset':blocks.append(b''.join(buf));buf=[]
  result={}
  for b in blocks:
   titles=re.findall(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',b);require(len(titles)==1,'PBE POTCAR title required')
   symbol=titles[0].decode().split('_')[0];require(symbol not in result,'Duplicate dataset');result[symbol]=b
  return result
 pots=split(p);require(set(pots)=={'Mg','Ni','H'},'Canonical element roster')
 lafile=base/'Mg_LaH3_Hrelay_endpoints_20261008/S0_base_static_retry_20261010/01_static/POTCAR'
 lap=split(lafile);require('La' in lap,'Licensed La dataset absent');pots['La']=lap['La']
 return pots

def write_model(name,s,group,metadata,pots=None):
 d=ROOT/name;d.mkdir(exist_ok=False)
 require(len(set(s['ids']))==len(s['ids']),'Duplicate atom map')
 kp='Common group mesh\n0\nGamma\n'+' '.join(map(str,PLAN['kpoints'][group]))+'\n0 0 0\n'
 for st in ('01_relax','02_static'):
  p=d/st;p.mkdir();(p/'INCAR').write_text(incar(s,name,st));(p/'KPOINTS').write_text(kp)
  if pots is not None:(p/'POTCAR').write_bytes(b''.join(pots[e] for e in dict.fromkeys(s['species'])))
 (d/'01_relax/POSCAR').write_text(poscar(s,name+' initial configuration; not a relaxed result'))
 save(d/'model.json',{'name':name,'group':group,'ids':s['ids'],'composition':{e:s['species'].count(e) for e in dict.fromkeys(s['species'])},'metadata':metadata,'initial_poscar_sha256':sha(d/'01_relax/POSCAR')})

def jobfile(mode,wall='24:00:00',array=None,gas=False):
 ranks=4 if gas else 64
 s='#!/bin/bash\n#SBATCH --partition=64c512g\n#SBATCH --nodes=1\n#SBATCH --ntasks='+str(ranks)+'\n#SBATCH --time='+wall+'\n'
 s+= '#SBATCH --mem=16G\n' if gas else '#SBATCH --exclusive\n'
 if array:s+='#SBATCH --array='+array+'\n'
 s+='#SBATCH --output=logs/%x-%A_%a.out\n#SBATCH --error=logs/%x-%A_%a.err\nset -euo pipefail\ncd "$SLURM_SUBMIT_DIR"\nmodule purge\nmodule load oneapi vasp/6.3.2-intel-2021.4.0\nexport OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\nulimit -s unlimited\nulimit -l unlimited\npython3 runner.py '+mode+'\n'
 return s

def prepare(preview=False):
 require(not ROOT.exists(),'Target exists; inspect it, never overwrite an earlier attempt')
 models=initial_models();pots=None if preview else local_potentials();ROOT.mkdir(parents=True);(ROOT/'logs').mkdir()
 for n,(s,g,meta) in models.items():write_model(n,s,g,meta,pots)
 shutil.copy2(HERE/'plan.json',ROOT/'plan.json');shutil.copy2(Path(__file__),ROOT/'runner.py')
 # Child construction reuses these six archived seeds; no repository access is needed on compute nodes.
 save(ROOT/'seed_structures.json',{n:s for n,(s,_,_) in models.items()})
 for fn,mode,wall,arr,gas in [('parents.slurm','--run-parent "$SLURM_ARRAY_TASK_ID"','24:00:00','0-4%2',False),('H2.slurm','--run-H2','00:30:00',None,True),('La.slurm','--run-La','24:00:00',None,False),('Mg.slurm','--run-Mg','12:00:00',None,False),('Hx.slurm','--run-Hx','24:00:00',None,False)]:
  p=ROOT/fn;p.write_text(jobfile(mode,wall,arr,gas));subprocess.run(['bash','-n',str(p)],check=True)
 save(ROOT/'manifest.json',{'status':'PREVIEW_ONLY_NO_POTCAR' if preview else 'PREPARED_NOT_SUBMITTED','initial_models':list(models),'maximum_models':15,'preview':preview,'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in ROOT.rglob('*') if p.is_file()}})
 print('PREPARED|'+str(ROOT)+'|6 seeds;9 dependent models not yet constructed',flush=True)

def read_output(d):
 s=parse((d/'POSCAR').read_text());F=E0=None;footer=ionic=ediff=False;it=None;on=False;buf=[];forces=[];n=0
 with (d/'OUTCAR').open(errors='replace') as f:
  for line in f:
   v=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
   if v and v.group(1)!=it:it=v.group(1);ediff=False
   if 'aborting loop because EDIFF is reached' in line:ediff=True
   if 'reached required accuracy' in line:ionic=True
   if 'General timing and accounting informations' in line:footer=True
   if 'free  energy   TOTEN' in line:F=float(line.split('=')[1].split()[0]);n+=1
   if 'energy(sigma->0)' in line:E0=float(line.split('energy(sigma->0)')[1].split('=')[1].split()[0])
   if 'TOTAL-FORCE (eV/Angst)' in line:on=True;buf=[];continue
   if on:
    cols=line.split()
    if len(cols)==6:
     try:buf.append(list(map(float,cols[3:6])));continue
     except ValueError:pass
    if buf:forces=buf;on=False
 require(len(forces)==len(s['ids']),'Missing final forces')
 fm=max(math.sqrt(sum(v*v for v,m in zip(a,mask) if m)) for a,mask in zip(forces,s['mask']))
 return {'F_eV':F,'E0_eV':E0,'last_scf_ediff':ediff,'normal_end':footer,'ionic_converged':ionic,'max_free_force_eV_A':fm,'ionic_steps':n}
def finished(name):
 d=ROOT/name;require((d/'VERIFIED_FINISHED').exists(),'Unverified parent '+name)
 meta=json.loads((d/'model.json').read_text());s=parse((d/'02_static/CONTCAR').read_text());s['ids']=meta['ids']
 return s,meta

def run_model(name):
 d=ROOT/name
 if (d/'VERIFIED_FINISHED').exists():return
 require(int(os.environ.get('SLURM_NTASKS','0'))==(4 if name=='H2_reference' else 64),'Only run in the intended Slurm allocation')
 (d/'.execution_started').mkdir(exist_ok=False);expected=parse((d/'01_relax/POSCAR').read_text());meta=json.loads((d/'model.json').read_text())
 for st in ('01_relax','02_static'):
  p=d/st
  if st=='02_static':shutil.copy2(d/'01_relax/CONTCAR',p/'POSCAR')
  with (p/'vasp.stdout').open('w') as o,(p/'vasp.stderr').open('w') as e:
   proc=subprocess.run(['mpirun','-np',os.environ['SLURM_NTASKS'],'vasp_std'],cwd=p,stdout=o,stderr=e)
  require(proc.returncode==0,'VASP failed: '+name+'/'+st)
  v=read_output(p);s=parse((p/'CONTCAR').read_text());s['ids']=meta['ids']
  require(s['species']==expected['species'] and s['mask']==expected['mask'],'Atom roster/constraints changed')
  require(max(abs(a-b) for row,ref in zip(s['cell'],expected['cell']) for a,b in zip(row,ref))<1e-7,'Cell changed')
  lim=.005 if name=='H2_reference' else (.015 if st=='01_relax' else .02)
  v['verified']=bool(v['normal_end'] and v['last_scf_ediff'] and all(v[k] is not None and math.isfinite(v[k]) for k in ('F_eV','E0_eV','max_free_force_eV_A')) and v['max_free_force_eV_A']<lim and (st=='02_static' or v['ionic_converged']))
  hs=[i for i,e in enumerate(s['species']) if e=='H'];v['min_H_H_A']=min((distance(s['xyz'][i],s['xyz'][j],s['cell']) for i,j in itertools.combinations(hs,2)),default=None)
  watch=meta['metadata'].get('background_H_ids',[])+meta['metadata'].get('transferred_H_ids',[])
  v['watched_H']=[{'id':j,'nearest_metals':[(dd,e,k) for dd,e,k in neighbours(s,s['xyz'][s['ids'].index(j)],exclude=j) if e!='H'][:5]} for j in watch if j in s['ids']]
  save(p/'validation.json',v);print('STAGE|'+name+'|'+st+'|'+json.dumps(v),flush=True);require(v['verified'],'Numerical check failed: '+name+'/'+st)
 (d/'VERIFIED_FINISHED').touch()

def copied_pots(name):
 p=ROOT/name/'01_relax/POTCAR';raw=p.read_bytes();blocks=[];buf=[]
 for line in raw.splitlines(keepends=True):
  buf.append(line)
  if line.strip()==b'End of Dataset':blocks.append(b''.join(buf));buf=[]
 out={}
 for b in blocks:
  m=re.search(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',b);require(m is not None,'Unrecognized potential');out[m.group(1).decode().split('_')[0]]=b
 return out

def make_vacancy(parent,name,deleted):
 s,meta=finished(parent);j=s['ids'].index(deleted);require(s['species'][j]=='H','Remove H only');v={'deleted_id':deleted,'xyz':s['xyz'][j]}
 if name=='L_vLa':
  _,_,anchor=neighbours(s,s['xyz'][j],'La')[0];v.update({'anchor_id':anchor,'vector':mic(s['xyz'][j],s['xyz'][s['ids'].index(anchor)],s['cell'])})
 new=subset(s,[i for i in range(len(s['ids'])) if i!=j]);write_model(name,new,meta['group'],{'parent':parent,'deleted_site':v},copied_pots(parent));run_model(name)

def La_pipeline():
 for parent,name,j in [('L_parent','L_vMg_near',39),('L_parent','L_vMg_far',31),('L_parent','L_vLa',55),('L_MgH2_reference','L_MgH2_vMg',39)]:make_vacancy(parent,name,j)
 s,meta=finished('L_vLa');site=meta['metadata']['deleted_site'];anchor=site['anchor_id'];a=s['xyz'][s['ids'].index(anchor)];target=[a[k]+site['vector'][k] for k in range(3)]
 cand=[]
 for donor in range(29,53):
  if donor not in s['ids']:continue
  j=s['ids'].index(donor);hh=[distance(target,x,s['cell']) for x,e,oid in zip(s['xyz'],s['species'],s['ids']) if e=='H' and oid!=donor]
  met=[distance(target,x,s['cell']) for x,e in zip(s['xyz'],s['species']) if e!='H']
  if min(hh)>1.10 and min(met)>1.25 and neighbours(s,s['xyz'][j],'Mg')[0][0]<2.4:cand.append((distance(target,s['xyz'][j],s['cell']),donor))
 require(cand,'No safe MgH2 donor for La-side vacancy');travel,donor=min(cand);s['xyz'][s['ids'].index(donor)]=target
 write_model('L_transfer',s,'L',{'parent':'L_vLa','transferred_H_ids':[donor],'acceptor_deleted_H':55,'displacement_A':travel,'no_H_added':True},copied_pots('L_vLa'));run_model('L_transfer')

def Mg_pipeline():
 s,meta=finished('Mg_parent');anchor=max(range(len(s['ids'])),key=lambda i:s['xyz'][i][2])
 for name,sep,height in [('Mg_H2',.75,2.30),('Mg_2H',2.0,1.65)]:
  n=copy.deepcopy(s);a=n['xyz'][anchor];points=[[a[0],a[1]+sign*sep/2,a[2]+height] for sign in (1,-1)]
  for p in points:require(neighbours(n,p,'Mg')[0][0]>1.35,'Mg endpoint collision')
  n['species']+=['H','H'];n['xyz']+=points;n['ids']+=[25,26];n['mask']+=[[True]*3,[True]*3]
  pots=copied_pots('Mg_parent');pots['H']=copied_pots('H2_reference')['H']
  write_model(name,n,'Mg',{'parent':'Mg_parent','anchor_id':s['ids'][anchor],'initial_H_H_A':sep,'not_a_barrier':True},pots);run_model(name)

def lowH_transfer(parent,name):
 s,meta=finished(parent);m=copy.deepcopy(meta['metadata']);moved=m.get('transferred_H_ids',[]);seed=json.loads((ROOT/'seed_structures.json').read_text())['X025_A'];candidates=[]
 for anchor,dx,dy,dz in PLAN['low_H_proxy']['candidate_pool']:
  ai=s['ids'].index(anchor);require(s['species'][ai]=='Ni','Bad acceptor anchor')
  if distance(s['xyz'][ai],seed['xyz'][seed['ids'].index(anchor)],s['cell'])>1.2:continue
  target=[s['xyz'][ai][k]+v for k,v in enumerate((dx,dy,dz))]
  for donor in list(range(37,73))+list(range(109,145)):
   if donor not in s['ids'] or donor in moved:continue
   j=s['ids'].index(donor)
   if neighbours(s,s['xyz'][j],'Ni')[0][0]<2.2 or neighbours(s,s['xyz'][j],'Mg')[0][0]>2.4:continue
   hh=[distance(target,x,s['cell']) for x,e,oid in zip(s['xyz'],s['species'],s['ids']) if e=='H' and oid!=donor]
   mg=neighbours(s,target,'Mg');ni=neighbours(s,target,'Ni')
   if min(hh)>1.40 and mg[0][0]>1.45 and mg[2][0]<2.70 and ni[0][0]>1.30 and ni[0][0]<2.10:candidates.append((distance(target,s['xyz'][j],s['cell']),donor,target,anchor))
 require(candidates,'No safe low-H receiving site; reassess instead of forcing H transfer')
 travel,donor,target,anchor=min(candidates,key=lambda a:a[0]);s['xyz'][s['ids'].index(donor)]=target
 m.update({'parent':parent,'transferred_H_ids':moved+[donor],'last_transfer':{'donor_id':donor,'anchor_Ni':anchor,'distance_A':travel},'nominal_receiver_x':(4+len(moved))/12,'total_H_conserved':75})
 require(s['species'].count('H')==75,'H conservation failed');write_model(name,s,'X',m,copied_pots(parent));run_model(name)

def Hx_pipeline():
 values={}
 for n in ('X025_A','X025_B'):
  s,meta=finished(n)
  for j in (145,146,147):require(neighbours(s,s['xyz'][s['ids'].index(j)],'Ni')[0][0]<2.35,'Low-H proxy did not retain its intended Ni-neighbor H; no forced continuation')
  values[n]=json.loads((ROOT/n/'02_static/validation.json').read_text())
 winner=min(values,key=lambda n:values[n]['F_eV']);require(winner==min(values,key=lambda n:values[n]['E0_eV']),'F/E0 low-H parent rankings disagree')
 save(ROOT/'low_H_selection.json',{'winner':winner,'energies':values,'limited_configurational_screen':True})
 lowH_transfer(winner,'X_transfer_one');lowH_transfer('X_transfer_one','X_transfer_two')

def submit():
 require(sys.platform.startswith('linux') and shutil.which('sbatch'),'Submit must run on a Slurm cluster, not the MacBook')
 if not (ROOT/'manifest.json').exists():prepare(False)
 manifest=json.loads((ROOT/'manifest.json').read_text());require(not manifest['preview'],'Preview has no licensed potentials; use a fresh production directory')
 if (ROOT/'submission.json').exists():print('EXISTING_SUBMISSION|'+(ROOT/'submission.json').read_text());return
 for rel,h in manifest['input_hashes'].items():require(sha(ROOT/rel)==h,'Pinned input changed: '+rel)
 (ROOT/'.submission_attempted').mkdir(exist_ok=False);jobs={}
 def send(key,opts):
  proc=subprocess.run(['sbatch','--parsable','--job-name=mechext_'+key]+opts,cwd=ROOT,capture_output=True,text=True,timeout=45)
  require(proc.returncode==0,'sbatch failed '+key+': '+proc.stderr[:300]);jid=proc.stdout.strip();require(re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?',jid),'Ambiguous sbatch response')
  jobs[key]=jid;save(ROOT/'submission.json',{'jobs':jobs,'all_submitted':False});print('SUBMITTED|'+key+'|'+jid,flush=True);return jid.split(';')[0]
 p=send('parents',['parents.slurm']);h=send('H2',['H2.slurm'])
 # Bounded second wave: no more than three solid jobs after all five parents succeed.
 for key,filename in [('La','La.slurm'),('Mg','Mg.slurm'),('Hx','Hx.slurm')]:
  dep='afterok:'+p+(':'+h if key=='Mg' else '')
  send(key,['--dependency='+dep,'--kill-on-invalid-dep=yes',filename])
 save(ROOT/'submission.json',{'jobs':jobs,'all_submitted':True,'maximum_models':15,'submitted_at':subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()})

def status():
 p=ROOT/'submission.json'
 if not p.exists():print('NOT_SUBMITTED');return
 d=json.loads(p.read_text());print(json.dumps(d,indent=2));ids=','.join(j.split(';')[0] for j in d['jobs'].values())
 subprocess.run(['sacct','-X','-n','-P','-j',ids,'-o','JobID,State,ExitCode,Elapsed,Start,End'],check=False)
 for name in PLAN['parents']+['H2_reference']+PLAN['La_children']+PLAN['Mg_children']+PLAN['Hx_children']:
  q=ROOT/name/'02_static/validation.json';print('MODEL|'+name+'|'+(q.read_text() if q.exists() else 'NO_FINAL_VALIDATION'))

def main():
 global ROOT,SOURCE_ROOT
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
 for flag in ('preview','prepare','submit','status','run-H2','run-La','run-Mg','run-Hx'):g.add_argument('--'+flag,action='store_true')
 g.add_argument('--run-parent',type=int);p.add_argument('--root',type=Path);p.add_argument('--source-root',type=Path);a=p.parse_args()
 if a.root:ROOT=a.root.expanduser().resolve()
 if a.source_root:SOURCE_ROOT=a.source_root.expanduser().resolve()
 if a.preview:prepare(True)
 elif a.prepare:prepare(False)
 elif a.submit:submit()
 elif a.status:status()
 elif a.run_parent is not None:
  require(0<=a.run_parent<5,'Parent index out of range');run_model(PLAN['parents'][a.run_parent])
 elif a.run_H2:run_model('H2_reference')
 elif a.run_La:La_pipeline()
 elif a.run_Mg:Mg_pipeline()
 elif a.run_Hx:Hx_pipeline()
if __name__=='__main__':main()
