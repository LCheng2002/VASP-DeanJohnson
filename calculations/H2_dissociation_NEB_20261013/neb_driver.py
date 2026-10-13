#!/usr/bin/env python3
"""Two requested H2 dissociation paths only. No existing calculation is changed.
Fresh common references/endpoints -> five-image NEB -> CI-NEB. No implicit retries.
"""
import argparse, importlib.util, json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import neb_core as c
HERE=Path(__file__).resolve().parent
PLAN=json.loads((HERE/'plan.json').read_text())
ROOT=Path.home()/PLAN['batch']
SYSTEMS=['Mg','Mg2Ni_region']
MODULE='vasp/6.3.2-vtst-intel-2021.4.0'

def source(key):
    spec=PLAN['sources'][key];p=Path.home()/spec['path'];c.require(c.sha(p)==spec['sha256'],'Source geometry changed: '+key)
    return c.parse(p.read_text())
def center_and_fix(s):
    s=c.clone(s);x=s['xyz'];s['xyz'][:,2]+=(s['cell'][2,2]-x[:,2].max()-x[:,2].min())/2
    pool=[i for i,e in enumerate(s['sp']) if e=='Mg'];fixed=sorted(pool,key=lambda i:s['xyz'][i,2])[:8]
    s['mask'][:]=True;s['mask'][fixed]=False
    c.require(len(s['sp']) in (24,42) and len(fixed)==8,'Unexpected substrate roster')
    c.require(c.summary_geometry(s)['nominal_vacuum_A']>14.,'Insufficient nominal vacuum')
    return s

def seeds():
    mg=center_and_fix(source('Mg_clean'));ni=center_and_fix(source('Ni_clean'))
    c.require(mg['sp']==['Mg']*24 and ni['sp']==['Mg']*36+['Ni']*6,'Wrong chemical system')
    c.require(np.allclose(mg['cell'],ni['cell'],atol=1e-8),'This comparison expects matched in-plane cells')
    h={'cell':np.eye(3)*20,'sp':['H','H'],'xyz':np.array([[10.,10.,9.625],[10.,10.,10.375]]),'mask':np.ones((2,3),bool)}
    return {'Mg':mg,'Mg2Ni_region':ni,'H2':h}

def licensed_potentials():
    p=Path.home()/PLAN['potcar']['path'];c.require(c.sha(p)==PLAN['potcar']['sha256'],'Canonical Mg/Ni/H potential changed')
    blocks={};buf=[]
    for line in p.read_bytes().splitlines(keepends=True):
        buf.append(line)
        if line.strip()==b'End of Dataset':
            raw=b''.join(buf);buf=[];m=re.search(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',raw);c.require(m is not None,'Unrecognized potential')
            key=m.group(1).decode();c.require(key in ('Mg','Ni','H') and key not in blocks,'Wrong/duplicate potential variant');blocks[key]=raw
    c.require(set(blocks)=={'Mg','Ni','H'},'POTCAR element roster changed');return blocks

def backend():
    p=shutil.which('vasp_std');c.require(p is not None,'VTST VASP not on PATH')
    p=Path(p).resolve();b=p.read_bytes();flags={w.decode():w in b for w in (b'LCLIMB',b'ICHAIN',b'IOPT')}
    c.require(all(flags.values()),'Loaded executable does not expose required VTST CI-NEB options; no submission')
    return {'module':MODULE,'path':str(p),'sha256':c.sha(p),'feature_markers':flags}

def nbands(system):return 8 if system=='H2' else 48 if system=='Mg' else 96

def write_calc(path,s,system,phase,pots,neb=None):
    path.mkdir(parents=True,exist_ok=False);gas=system=='H2'
    opts=c.incar(system,phase,nbands(system),gas,neb)
    (path/'INCAR').write_text(c.input_text(opts,s['sp']))
    (path/'KPOINTS').write_text('Common comparison mesh\n0\nGamma\n'+('1 1 1' if gas else '3 6 1')+'\n0 0 0\n')
    (path/'POSCAR').write_text(c.poscar(s,system+' '+phase))
    if pots is not None:(path/'POTCAR').write_bytes(b''.join(pots[e] for e in dict.fromkeys(s['sp'])))
    c.save(path/'input_metadata.json',{'system':system,'phase':phase,'geometry':c.summary_geometry(s),'incar':opts,'poscar_sha256':c.sha(path/'POSCAR'),'potcar_sha256':c.sha(path/'POTCAR') if (path/'POTCAR').exists() else None})

def split_local_pots(system):
    return {e:(ROOT/'licensed'/e/'POTCAR').read_bytes() for e in ('Mg','Ni','H')}

def slurm(mode,ntasks,wall,array=None,gas=False):
    x='#!/bin/bash\n#SBATCH --partition=64c512g\n#SBATCH --nodes=1\n#SBATCH --ntasks='+str(ntasks)+'\n#SBATCH --time='+wall+'\n'
    x+='#SBATCH --mem=16G\n' if gas else '#SBATCH --exclusive\n'
    if array:x+='#SBATCH --array='+array+'\n'
    x+='#SBATCH --output=logs/%x-%A_%a.out\n#SBATCH --error=logs/%x-%A_%a.err\nset -euo pipefail\ncd "$SLURM_SUBMIT_DIR"\nmodule purge\nmodule load oneapi '+MODULE+'\nexport OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\nulimit -s unlimited\nulimit -l unlimited\npython3 neb_driver.py '+mode+'\n'
    return x

def prepare(preview=False):
    c.require(not ROOT.exists(),'Target already exists; do not overwrite any earlier submission')
    ss=seeds();caps=None if preview else backend();pots=None if preview else licensed_potentials()
    ni_IS=source('Ni_molecule');ni_FS=source('Ni_product');c.endpoint_check(ni_IS,'molecule');c.endpoint_check(ni_FS,'product')
    ROOT.mkdir(parents=True);(ROOT/'logs').mkdir();(ROOT/'sources').mkdir()
    for key in ('Mg_clean','Ni_clean','Ni_molecule','Ni_product'):
        (ROOT/'sources'/(key+'.vasp')).write_text(c.poscar(source(key),key+' audited geometry'))
    for fn in ('neb_driver.py','neb_core.py','plan.json'):shutil.copy2(HERE/fn,ROOT/fn)
    if pots is not None:
        for e,b in pots.items():p=ROOT/'licensed'/e;p.mkdir(parents=True);(p/'POTCAR').write_bytes(b)
    for system in SYSTEMS+['H2']:write_calc(ROOT/system/'clean'/'relax',ss[system],system,'relax',pots)
    jobs=[('references.slurm','--reference "$SLURM_ARRAY_TASK_ID"',64,'12:00:00','0-1',False),('H2.slurm','--reference-H2',4,'00:30:00',None,True),('endpoints.slurm','--endpoints "$SLURM_ARRAY_TASK_ID"',64,'24:00:00','0-1',False),('NEB_Mg.slurm','--neb 0',60,'24:00:00',None,False),('NEB_Mg2Ni.slurm','--neb 1',60,'24:00:00',None,False)]
    for fn,mode,n,w,arr,gas in jobs:
        p=ROOT/fn;p.write_text(slurm(mode,n,w,arr,gas));subprocess.run(['bash','-n',str(p)],check=True)
    c.save(ROOT/'manifest.json',{'preview':preview,'backend':caps,'system_geometries':{k:c.summary_geometry(s) for k,s in ss.items()},'sources':PLAN['sources'],'explicit_scope':'2 paths;2 clean slabs+H2;4 endpoints;5 moving images per path;NEB then CI-NEB','input_hashes':{str(p.relative_to(ROOT)):c.sha(p) for p in ROOT.rglob('*') if p.is_file()}})
    print('PREPARED|2 references+H2;4 endpoints deferred;2 CI-NEB paths gated',flush=True)

def check_runtime(ranks):
    c.require(int(os.environ.get('SLURM_NTASKS','0'))==ranks,'Wrong Slurm allocation')
    old=json.loads((ROOT/'manifest.json').read_text())['backend'];new=backend();c.require(new['sha256']==old['sha256'],'Executable changed from preflight')

def vasp(folder,ranks):
    with (folder/'vasp.stdout').open('w') as o,(folder/'vasp.stderr').open('w') as e:
        p=subprocess.run(['mpirun','-np',str(ranks),'vasp_std'],cwd=folder,stdout=o,stderr=e)
    c.require(p.returncode==0,'VASP execution failed: '+str(folder))

def geometry_constant(s,ref):
    c.require(s['sp']==ref['sp'] and np.allclose(s['cell'],ref['cell'],atol=1e-8) and np.array_equal(s['mask'],ref['mask']),'Composition, cell or mask changed')
    c.require(np.max(np.abs(c.mic(s['xyz']-ref['xyz'],s['cell'])[~ref['mask']]),initial=0)<1e-7,'Fixed substrate atoms moved')

def run_relax_static(base,system,kind):
    if (base/'VERIFIED_FINISHED').exists():return
    (base/'.execution_started').mkdir(exist_ok=False);ranks=4 if system=='H2' else 64;check_runtime(ranks)
    ref=c.parse((base/'relax/POSCAR').read_text());pots=split_local_pots(system)
    for phase in ('relax','static'):
        folder=base/phase
        if phase=='static':
            last=c.parse((base/'relax/CONTCAR').read_text());write_calc(folder,last,system,'static',pots)
        vasp(folder,ranks);result=c.output(folder,ref['mask']);end=c.parse((folder/'CONTCAR').read_text());geometry_constant(end,ref)
        limit=.005 if system=='H2' else .01
        result['verified']=bool(result['normal_end'] and result['last_ediff'] and result['F_eV'] is not None and result['E0_eV'] is not None and np.isfinite(result['F_eV']) and result['max_free_DFT_force']<limit and (phase=='static' or result['ionic_converged']))
        result['geometry']=c.summary_geometry(end);c.save(folder/'validation.json',result)
        print('STAGE|'+system+'|'+str(base.relative_to(ROOT))+'|'+phase+'|'+json.dumps({k:v for k,v in result.items() if k!='geometry'}),flush=True)
        c.require(result['verified'],'Endpoint/reference numerical verification failed: '+str(folder))
    if kind in ('molecule','product'):c.endpoint_check(end,kind)
    (base/'VERIFIED_FINISHED').touch()

def verified_geometry(base):
    c.require((base/'VERIFIED_FINISHED').exists(),'Required parent has not passed: '+str(base))
    return c.parse((base/'static/CONTCAR').read_text())

def endpoint_seed(system,kind,parent):
    s=c.clone(parent)
    if system=='Mg2Ni_region':
        key='Ni_molecule' if kind=='molecule' else 'Ni_product';old=c.parse((ROOT/'sources'/(key+'.vasp')).read_text());anchor=40
        c.require(old['sp'][anchor]=='Ni','Ni anchor identity changed')
        pair=old['xyz'][-2:].copy();pair[1]=pair[0]+c.mic(pair[1]-pair[0],old['cell']);center=pair.mean(axis=0)
        shift=parent['xyz'][anchor]+c.mic(center-old['xyz'][anchor],parent['cell'])-center;pair+=shift
    else:
        x=parent['xyz'];top=np.where(x[:,2]>x[:,2].max()-.15)[0];target=np.array([parent['cell'][0,0]/2,parent['cell'][1,1]/2,x[:,2].max()])
        anchor=int(min(top,key=lambda i:np.linalg.norm(c.mic(x[i]-target,parent['cell'])[:2])))
        sep,height=(.75,2.30) if kind=='molecule' else (2.0,1.65)
        center=x[anchor]+[0,0,height];pair=np.array([center+[0,-sep/2,0],center+[0,sep/2,0]])
    s['sp']+=['H','H'];s['xyz']=np.vstack([s['xyz'],pair]);s['mask']=np.vstack([s['mask'],np.ones((2,3),bool)])
    c.endpoint_check(s,kind);return s

def endpoints(index):
    system=SYSTEMS[index];check_runtime(64);base=ROOT/system;parent=verified_geometry(base/'clean');pots=split_local_pots(system)
    (base/'.endpoints_started').mkdir(exist_ok=False)
    for kind in ('molecule','product'):
        s=endpoint_seed(system,kind,parent);write_calc(base/kind/'relax',s,system,'relax',pots);run_relax_static(base/kind,system,kind)
    (base/'ENDPOINTS_VERIFIED').touch()

def stage_band(folder,images,system,mode):
    pots=split_local_pots(system);write_calc(folder,images[0],system,'neb',pots,mode)
    for i,s in enumerate(images):
        q=folder/('%02d'%i);q.mkdir();(q/'POSCAR').write_text(c.poscar(s,system+' '+mode+' image '+str(i)))
    for index,kind in ((0,'molecule'),(6,'product')):
        for fn in ('OUTCAR','OSZICAR'):shutil.copy2(ROOT/system/kind/'static'/fn,folder/('%02d'%index)/fn)

def band_results(folder,images,mode):
    values=[];limit=.03 if mode=='ci' else .05
    for i in range(1,6):
        q=folder/('%02d'%i);v=c.output(q,images[i]['mask']);s=c.parse((q/'CONTCAR').read_text());geometry_constant(s,images[i])
        v['geometry']=c.summary_geometry(s)
        v['verified']=bool(v['normal_end'] and v['last_ediff'] and v['neb_marker'] and v['lclimb']==(mode=='ci') and v['vtst_force_max'] is not None and v['vtst_force_max']<limit and v['F_eV'] is not None and np.isfinite(v['F_eV']))
        c.save(q/'validation.json',v);values.append(v)
    report={'mode':mode,'force_tolerance_eV_A':limit,'images':values,'verified':all(v['verified'] for v in values)}
    c.save(folder/'validation.json',report);c.require(report['verified'],'NEB/CI-NEB did not meet all image checks; results retained, no barrier claim')
    return values

def neb(index):
    system=SYSTEMS[index];base=ROOT/system;check_runtime(60);c.require((base/'ENDPOINTS_VERIFIED').exists(),'Endpoints unverified')
    (base/'.neb_started').mkdir(exist_ok=False);a=verified_geometry(base/'molecule');b=verified_geometry(base/'product');c.endpoint_check(a,'molecule');c.endpoint_check(b,'product')
    images,interp=c.initial_band(a,b,5);c.save(base/'path_preflight.json',interp)
    coarse=base/'neb_coarse';stage_band(coarse,images,system,'coarse');vasp(coarse,60);band_results(coarse,images,'coarse')
    refined=[images[0]]+[c.parse((coarse/('%02d'%i)/'CONTCAR').read_text()) for i in range(1,6)]+[images[-1]]
    fine=base/'ci_neb';stage_band(fine,refined,system,'ci');vasp(fine,60);values=band_results(fine,refined,'ci')
    ends=[json.loads((base/k/'static/validation.json').read_text()) for k in ('molecule','product')]
    F=[ends[0]['F_eV']]+[v['F_eV'] for v in values]+[ends[1]['F_eV']];E0=[ends[0]['E0_eV']]+[v['E0_eV'] for v in values]+[ends[1]['E0_eV']]
    top=int(np.argmax(F));c.require(0<top<6,'No resolved interior activation maximum; do not equate with a zero barrier')
    c.require(values[top-1]['max_free_DFT_force']<.04,'Highest image not stationary in raw atomic force')
    ref=json.loads((base/'clean/static/validation.json').read_text());gasfile=ROOT/'H2/clean/static/validation.json';gas=json.loads(gasfile.read_text()) if gasfile.is_file() else None
    r={'status':'CI_NEB_NUMERICALLY_VERIFIED_TS_VIBRATIONS_NOT_DONE','system':system,'image_F_eV':F,'image_E0_eV':E0,'highest_image':top,'Ea_from_adsorbed_F_eV':F[top]-F[0],'Ea_reverse_F_eV':F[top]-F[-1],'reaction_F_eV':F[-1]-F[0],'Ea_from_adsorbed_E0_eV':E0[top]-E0[0],'Ea_from_gas_F_eV':F[top]-ref['F_eV']-gas['F_eV'] if gas and gas.get('verified') else None,'Ea_from_gas_E0_eV':E0[top]-ref['E0_eV']-gas['E0_eV'] if gas and gas.get('verified') else None,'limitations':PLAN['limitations']}
    c.save(base/'barrier_summary.json',r);(base/'NEB_VERIFIED').touch();print('BARRIER_REPORT|'+json.dumps(r),flush=True)

def submit():
    if not (ROOT/'manifest.json').exists():prepare(False)
    m=json.loads((ROOT/'manifest.json').read_text());c.require(not m['preview'],'Cannot submit preview without licensed potentials')
    ledger=ROOT/'submission.json'
    if ledger.exists():
        d=json.loads(ledger.read_text());print('EXISTING_SUBMISSION|'+json.dumps(d),flush=True);c.require(d.get('all_submitted'),'Partial prior submission requires manual review');return
    for rel,h in m['input_hashes'].items():c.require(c.sha(ROOT/rel)==h,'Pinned input changed: '+rel)
    (ROOT/'.submission_attempted').mkdir(exist_ok=False);jobs={}
    def send(key,opts):
        c.save(ledger,{'jobs':jobs,'attempting':key,'all_submitted':False})
        p=subprocess.run(['sbatch','--parsable','--job-name=h2neb_'+key]+opts,cwd=ROOT,capture_output=True,text=True,timeout=45)
        c.require(p.returncode==0,'sbatch rejected '+key+': '+p.stderr[:400]);jid=p.stdout.strip();c.require(re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?',jid),'Ambiguous Slurm response')
        jobs[key]=jid;c.save(ledger,{'jobs':jobs,'all_submitted':False});print('SUBMITTED|'+key+'|'+jid,flush=True);return jid.split(';')[0]
    refs=send('refs',['references.slurm']);gas=send('H2',['H2.slurm'])
    ep=send('endpoints',['--dependency=afterok:'+refs,'--kill-on-invalid-dep=yes','endpoints.slurm'])
    send('Mg',['--dependency=afterok:'+ep+'_0:'+gas,'--kill-on-invalid-dep=yes','NEB_Mg.slurm'])
    send('Mg2Ni',['--dependency=afterok:'+ep+'_1:'+gas,'--kill-on-invalid-dep=yes','NEB_Mg2Ni.slurm'])
    c.save(ledger,{'jobs':jobs,'all_submitted':True,'submitted_at':subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip(),'path_count':2,'moving_images_each':5,'neb_mpi_ranks':60,'no_original_jobs_modified':True})

def status():
    print('CHECKED_AT|'+subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip())
    p=ROOT/'submission.json'
    if not p.is_file():print('NO_SUBMISSION');return
    d=json.loads(p.read_text());print('SUBMISSION|'+json.dumps(d));ids=','.join(x.split(';')[0] for x in d['jobs'].values())
    for cmd in (['squeue','-r','-h','-j',ids,'-o','%i|%j|%T|%M|%R'],['sacct','-X','-n','-P','-j',ids,'-o','JobID,JobName%25,State,ExitCode,Elapsed,Start,End']):
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=30);print(p.stdout)
    for system in SYSTEMS+['H2']:
        for stage in ('clean','molecule','product'):
            p=ROOT/system/stage/'static/validation.json'
            if p.is_file():
                v=json.loads(p.read_text());print('ENDPOINT|'+system+'|'+stage+'|'+json.dumps({k:v.get(k) for k in ('F_eV','max_free_DFT_force','verified')}))
        p=ROOT/system/'barrier_summary.json'
        if p.is_file():print('BARRIER|'+p.read_text().replace('\n',''))
        else:print('BARRIER|'+system+'|NOT_YET_VERIFIED')

def main():
    global ROOT
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    for v in ('prepare','preview','submit','status','reference-H2'):g.add_argument('--'+v,action='store_true')
    for v in ('reference','endpoints','neb'):g.add_argument('--'+v,type=int,choices=[0,1])
    p.add_argument('--preview-root',type=Path);a=p.parse_args()
    if a.preview_root:c.require(a.preview,'Custom root only allowed for preview');ROOT=a.preview_root.resolve()
    if a.prepare:prepare()
    elif a.preview:prepare(True)
    elif a.submit:submit()
    elif a.status:status()
    elif a.reference_H2:run_relax_static(ROOT/'H2/clean','H2','clean')
    elif a.reference is not None:
        s=SYSTEMS[a.reference];run_relax_static(ROOT/s/'clean',s,'clean')
    elif a.endpoints is not None:endpoints(a.endpoints)
    elif a.neb is not None:neb(a.neb)
if __name__=='__main__':main()
