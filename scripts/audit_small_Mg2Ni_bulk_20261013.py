#!/usr/bin/env python3
"""Read-only audit of a known converged Mg2Ni reference, H2 and cluster state.
No licensed POTCAR bytes and no auth secrets are read or exported.
"""
from pathlib import Path
import hashlib,json,datetime,re,subprocess,os
home=Path.home()
ROOT=home/'Mg2Ni_H_transfer_bulk_screen_20261009_supplement_v2/supplement_20261009/02_Mg2Ni'
H2=home/'MgNi_LaH3_extensions_20261012/H2_reference'
POTS=home/'Mg2Ni_Mg_H_thermo_validation_20261010/A-H-Ni_Mg_Mg2Ni/01_relax/POTCAR'
report={'checked_at':datetime.datetime.now().astimezone().isoformat(),'bulk_root':str(ROOT.relative_to(home)),'bulk':{},'h2':{},'potential':{},'small_batch_exists':(home/'Mg2Ni_lowH_bulk_20261013').exists()}
for name,folder,dest in [('Mg2Ni_01_relax',ROOT/'01_relax',report['bulk']),('Mg2Ni_02_static',ROOT/'02_static',report['bulk']),('H2_02_static',H2/'02_static',report['h2'])]:
 val={}
 for f in ('POSCAR','CONTCAR','INCAR','KPOINTS','validation.json'):
  p=folder/f
  if p.is_file() and p.stat().st_size<40000:
   raw=p.read_text(errors='replace')
   val[f]={'sha256':hashlib.sha256(raw.encode()).hexdigest(),'text':raw}
   if f=='POSCAR':val['composition']={'symbols':raw.splitlines()[5].split(),'counts':raw.splitlines()[6].split()}
 p=folder/'OUTCAR'
 if p.is_file():
  last=None;scf=False;end=False;ion=False;iters=0;force=None
  with p.open(errors='replace') as q:
   reading=False; forces=[]
   for line in q:
    if 'free  energy   TOTEN' in line:last=float(line.split('=')[1].split()[0]);iters+=1
    if 'aborting loop because EDIFF is reached' in line:scf=True
    if 'General timing and accounting informations' in line:end=True
    if 'reached required accuracy' in line:ion=True
    if 'TOTAL-FORCE (eV/Angst)' in line:reading=True;forces=[];continue
    if reading:
     row=line.split()
     if len(row)==6:
      try:forces.append([float(v) for v in row[3:]]);continue
      except ValueError:pass
     if forces:force=max(sum(x*x for x in v)**.5 for v in forces);reading=False
  val['OUTCAR_check']={'F_eV':last,'SCF_ediff_mark':scf,'normal_end':end,'ionic_converged':ion,'ionic_steps':iters,'max_force_all_eV_A':force}
 dest[name]=val
if POTS.is_file():
 with POTS.open(errors='replace') as p:titles=[line.strip() for line in p if 'TITEL' in line]
 report['potential']={'exists':True,'sha256':hashlib.sha256(POTS.read_bytes()).hexdigest(),'titles':titles}
else:report['potential']={'exists':False}
q=subprocess.run(['squeue','-u',os.environ['USER'],'-h','-o','%i|%j|%T|%M'],text=True,capture_output=True,timeout=25)
report['queue']={'rc':q.returncode,'tail':q.stdout[-2000:]}
print('AUDIT|'+json.dumps(report,ensure_ascii=False))
