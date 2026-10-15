#!/usr/bin/env python3
"""READ ONLY: preflight the timed-out La branch. Print metadata only; no POTCAR data."""
import datetime,hashlib,importlib.util,json,math,re,subprocess
from pathlib import Path
HOME=Path.home()
ROOT=HOME/'MgNi_LaH3_extensions_20261012'
REQ=['L_vMg_near','L_vMg_far','L_vLa','L_MgH2_vMg','L_transfer']
def req(flag,msg):
 if not flag:raise RuntimeError("AUDIT_FAILED|"+msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_runner():
 path=ROOT/'runner.py';req(path.exists(),'missing old runner.py')
 p=json.loads((ROOT/'manifest.json').read_text())
 req(sha(path)==p['input_hashes']['runner.py'],'old driver hash mismatch')
 spec=importlib.util.spec_from_file_location('la_original_runner',path)
 r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
 return r
def slab_summary(path,runner,ids):
 raw=path.read_text(errors='replace')
 s=runner.parse(raw);s['ids']=ids
 req(len(s['ids'])==len(s['species']),'IDs mismatch for '+str(path))
 return {'elements':{e:s['species'].count(e) for e in sorted(set(s['species']))},
         'atoms':len(s['species']),'sha256':sha(path),
         'cell_vectors':s['cell'],'z_min_max':[min(v[2] for v in s['xyz']),max(v[2] for v in s['xyz'])],
         'max_coordinate_abs':max(abs(v) for xyz in s['xyz'] for v in xyz),
         'bottom_fixed_count':sum(not any(m) for m in s['mask'])}
def out_status(p):
 d={'size':p.stat().st_size if p.exists() else None}; ecount=0
 if not p.is_file():return d
 import collections
 tail=collections.deque(maxlen=260);iterations=[];final_force=None;ions=0
 with p.open(errors='replace') as f:
  force=False;buf=[]
  for line in f:
   match=re.search(r'Iteration\s+(\d+)\(\s*(\d+)\)',line)
   if match:ions=max(ions,int(match.group(1)))
   if 'free  energy   TOTEN' in line:
    ecount+=1
    try:d['last_F_eV']=float(line.split('=')[1].split()[0])
    except Exception:pass
   if 'reached required accuracy' in line:d['ionic_converged_marker']=True
   if 'General timing and accounting informations' in line:d['normal_end']=True
   if 'TOTAL-FORCE (eV/Angst)' in line:force=True;buf=[];continue
   if force:
    cols=line.split()
    if len(cols)==6:
     try:buf.append([float(x) for x in cols[3:6]]);continue
     except ValueError:pass
    if buf:final_force=buf;force=False
   tail.append(line)
 d['max_ionic_iter_seen']=ions;d['ionic_energy_rows']=ecount
 if final_force:
  d['force_vector_rows']=len(final_force)
  d['max_total_force_A']=max(math.sqrt(sum(v*v for v in x)) for x in final_force)
 d['last_scf_ediff_seen_in_tail']='aborting loop because EDIFF is reached' in ''.join(tail)
 return d
def main():
 runner=read_runner()
 ledger=ROOT/'recovery_La_Mg_only_20261013/submission.json'
 req(ledger.is_file(),'missing recovery job ledger')
 rec=json.loads(ledger.read_text())
 req(rec.get('jobs',{}).get('La')=='63636883','Unexpected La job id')
 acc=subprocess.run(['sacct','-X','-n','-P','-j','63636883','-o','JobID,JobName%38,State,ExitCode,Elapsed,Timelimit,Start,End'],capture_output=True,text=True,timeout=25)
 req(acc.returncode==0,'sacct unavailable')
 status={x.split('|')[0]:x.split('|') for x in acc.stdout.splitlines() if len(x.split('|'))>=4}
 req('63636883' in status,'No old job accounting')
 req(status['63636883'][2]=='TIMEOUT','Old La branch did not timeout; refuse automatic continuation')
 queue=subprocess.run(['squeue','-h','-u',str(ROOT.owner()) if False else 'deanfilan','-o','%i|%j|%T'],capture_output=True,text=True,timeout=25)
 req(queue.returncode==0,'queue unavailable')
 active=[x for x in queue.stdout.splitlines() if 'mechext_La' in x or 'LaH3_La_resume' in x]
 req(not active,'Another active La continuation exists')
 parents={}
 for n in ('L_parent','L_MgH2_reference'):
  d=ROOT/n; req((d/'VERIFIED_FINISHED').exists(),n+' missing VERIFIED_FINISHED')
  meta=json.loads((d/'model.json').read_text())
  parents[n]={'metadata_ids_count':len(meta['ids']),'composition':meta.get('composition')}
  parents[n]['h39_present']=39 in meta['ids'] and runner.parse((d/'02_static/CONTCAR').read_text())['species'][meta['ids'].index(39)]=='H'
  req(parents[n]['h39_present'] or n!='L_MgH2_reference','L_MgH2_reference H39 not mapped')
 models={}
 for n in REQ:
  d=ROOT/n;entry={'exists':d.exists(),'finished':(d/'VERIFIED_FINISHED').is_file()}
  if n in ('L_vMg_near','L_vMg_far'):
   req(entry['finished'],'Expected finished model not verified '+n)
  elif n in ('L_MgH2_vMg','L_transfer'):req(not d.exists(),'Unstarted model unexpectedly exists: '+n)
  if d.exists():
   meta=json.loads((d/'model.json').read_text());entry['model_json_sha256']=sha(d/'model.json')
   entry['metadata']=meta.get('metadata',{})
   entry['ids_count']=len(meta['ids'])
   for stage in ('01_relax','02_static'):
    sd=d/stage; info={'exists':sd.exists()}
    if sd.is_dir():
     for filename in ('INCAR','KPOINTS','POTCAR','POSCAR','CONTCAR','WAVECAR','CHGCAR','OUTCAR','OSZICAR','vasp.stdout','vasp.stderr'):
      p=sd/filename
      if p.is_file():info[filename]={'size':p.stat().st_size}
     if (sd/'INCAR').is_file():info['INCAR_text']=(sd/'INCAR').read_text()[:3800]
     if (sd/'CONTCAR').is_file():info['geometry']=slab_summary(sd/'CONTCAR',runner,meta['ids'])
     if (sd/'OUTCAR').is_file():info['outcar']=out_status(sd/'OUTCAR')
    entry[stage]=info
  models[n]=entry
 report={'checked_at':datetime.datetime.now().astimezone().isoformat(),'old_slurm':acc.stdout.strip(),
         'active_related_queue':active,'parents':parents,'models':models,'old_driver_sha256':sha(ROOT/'runner.py')}
 print(json.dumps(report,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
