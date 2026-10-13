#!/usr/bin/env python3
"""Resume ONLY the two cancelled La/Mg branches of MgNi_LaH3_extensions_20261012.

No re-optimization of completed parents, no Hx submission/restart, no new physics
or changes to VASP inputs, and no reuse of invalid parent-array afterok dependency.
Records unambiguous sbatch acknowledgments. If re-run, never blindly resubmits.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
ROOT = Path.home()/'MgNi_LaH3_extensions_20261012'
FOLDER = ROOT/'recovery_La_Mg_only_20261013'
LEDGER = FOLDER/'submission.json'
NAMES = {
    'La': ['L_vMg_near','L_vMg_far','L_vLa','L_MgH2_vMg','L_transfer'],
    'Mg': ['Mg_H2','Mg_2H']
}
REQUIRED = {
    'La': ['L_parent','L_MgH2_reference'],
    'Mg': ['Mg_parent','H2_reference']
}
ORIGINAL_IDS = {'parents':'63585911','H2':'63585912','La':'63585913','Mg':'63585914','Hx':'63585915'}

def require(ok,reason):
    if not ok: raise RuntimeError('PRECHECK_FAIL|'+reason)
def now():
    return datetime.datetime.now().astimezone().isoformat()
def ledger_save(d):
    temp=LEDGER.with_suffix('.json.tmp')
    temp.write_text(json.dumps(d,indent=2,ensure_ascii=False))
    temp.replace(LEDGER)
def sacct(ids):
    p=subprocess.run(['sacct','-X','-n','-P','-j',','.join(ids),'-o','JobID,JobName%25,State,ExitCode,Elapsed,Start,End'],capture_output=True,text=True,timeout=30)
    require(p.returncode==0,'sacct failed: '+p.stderr[:250])
    return p.stdout
def records(raw):
    return {r.split('|')[0]:r.split('|') for r in raw.splitlines() if len(r.split('|'))>=4}
def validated_model(name):
    p=ROOT/name
    require((p/'VERIFIED_FINISHED').is_file(),name+' missing VERIFIED_FINISHED')
    meta=json.loads((p/'model.json').read_text())
    stages={}
    for phase,force_lim in [('01_relax',.015),('02_static',.02)]:
        d=p/phase
        v=json.loads((d/'validation.json').read_text())
        require(bool(v.get('verified')) and bool(v.get('last_scf_ediff')) and bool(v.get('normal_end')),name+'/'+phase+' electronic convergence/normal end gate')
        if phase=='01_relax': require(bool(v.get('ionic_converged')),name+' ionic relaxation gate')
        require(v['max_free_force_eV_A']<force_lim,name+'/'+phase+' force gate')
        require((d/'CONTCAR').is_file() and (d/'CONTCAR').stat().st_size>100,name+'/'+phase+' missing CONTCAR')
        stages[phase]={'F_eV':v.get('F_eV'),'force_eV_A':v['max_free_force_eV_A']}
    return {'composition':meta['composition'],'stages':stages}
def check():
    require(ROOT.is_dir(),'Original production directory missing')
    original=json.loads((ROOT/'submission.json').read_text())
    require(original.get('all_submitted') and all(original['jobs'].get(k)==v for k,v in ORIGINAL_IDS.items()),'Original ledger differs from expected job IDs')
    manifest=json.loads((ROOT/'manifest.json').read_text())
    require(not manifest.get('preview'),'Original production was preview only')
    driver=ROOT/'runner.py'
    require(driver.is_file(),'Missing original production runner')
    pinned=manifest['input_hashes'].get('runner.py')
    require(pinned and hashlib.sha256(driver.read_bytes()).hexdigest()==pinned,'Original production driver changed')
    for fn in ('La.slurm','Mg.slurm'):
        require((ROOT/fn).is_file(),fn+' missing')
        raw=(ROOT/fn).read_text()
        require('runner.py --run-'+fn.split('.')[0] in raw,fn+' runs unexpected branch')
        require('Hx' not in raw and '--dependency' not in raw,fn+' unexpected Hx/dependency contents')
    snapshot=sacct(['63585911','63585912','63585913','63585914','63585915'])
    state=records(snapshot)
    print('ORIGINAL_SACCT_CHECK|'+snapshot.replace('\n',' || '),flush=True)
    for ident in ('63585911_0','63585911_1','63585911_2','63585912'):
        v=state.get(ident)
        require(v and v[2]=='COMPLETED' and v[3]=='0:0','Required parent Slurm not successful: '+ident+' '+str(v))
    for ident in ('63585913','63585914','63585915'):
        v=state.get(ident)
        require(v and v[2]=='CANCELLED','Original downstream not cancelled: '+ident+' '+str(v))
    verification={}
    for name in sorted(set(sum(list(REQUIRED.values()),[]))):
        verification[name]=validated_model(name)
    return {'checked_at':now(),'slurm':snapshot,'parents_verified':verification}
def submit():
    if FOLDER.exists():
        print('RECOVERY_EXISTS|'+str(FOLDER))
        if LEDGER.is_file():print('RECOVERY_SUBMISSION|'+LEDGER.read_text())
        else:print('RECOVERY_LOCK_WITHOUT_LEDGER|Inspect previous uncertain submission; no duplicate sbatch')
        return
    state=check()
    for branch in ('La','Mg'):
        for name in NAMES[branch]:
            require(not (ROOT/name).exists(),'Child already exists: '+name+'; do not blindly duplicate')
    FOLDER.mkdir(mode=0o700)
    d={'scope':'La and Mg only; NO Hx','created_at':now(),'original':ORIGINAL_IDS,'preflight':state,'jobs':{},'completed_submission':False}
    ledger_save(d)
    print('VALIDATED_PARENTS|'+json.dumps(state['parents_verified'],ensure_ascii=False),flush=True)
    for branch in ('La','Mg'):
        fn=branch+'.slurm'
        cmd=['sbatch','--parsable','--job-name=mechext_'+branch+'_recovery','--chdir='+str(ROOT),fn]
        p=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=45)
        require(p.returncode==0,'sbatch '+branch+': '+p.stderr[:350])
        ident=p.stdout.strip().split(';')[0]
        require(re.fullmatch(r'[0-9]+',ident),'Invalid sbatch response for '+branch+': '+p.stdout[:100])
        d['jobs'][branch]=ident
        ledger_save(d)
        print('SUBMITTED|'+branch+'|'+ident,flush=True)
    d['completed_submission']=True
    d['submitted_at']=now()
    ledger_save(d)
    print('SUBMISSION_COMPLETE|'+json.dumps({'jobs':d['jobs'],'Hx_submitted':False,'root':str(ROOT)},ensure_ascii=False),flush=True)
def status():
    d=check()
    print('ORIGINAL_PARENT_AND_CANCELLED_STATUS|'+d['slurm'])
    if LEDGER.is_file():
        j=json.loads(LEDGER.read_text())
        print('RECOVERY_SUBMISSION|'+json.dumps({'jobs':j['jobs'],'complete':j['completed_submission'],'scope':j['scope']},ensure_ascii=False))
        if j['jobs']:print('RECOVERY_SACCT|'+sacct(list(j['jobs'].values())))
    else:print('NO_RECOVERY_LEDGER')
    for branch in ('La','Mg'):
        for name in NAMES[branch]:
            p=ROOT/name
            st='VERIFIED_FINISHED' if (p/'VERIFIED_FINISHED').is_file() else 'EXISTS_UNVERIFIED' if p.exists() else 'NOT_YET_CONSTRUCTED'
            print('CHILD|'+name+'|'+st)
    print('HX_BRANCH|NOT_SUBMITTED_BY_RECOVERY')
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    gr=parser.add_mutually_exclusive_group(required=True)
    gr.add_argument('--submit',action='store_true')
    gr.add_argument('--status',action='store_true')
    opts=parser.parse_args()
    if opts.submit:submit()
    else:status()
