#!/usr/bin/env python3
"""Copy authorized local PAW datasets into a dedicated potential cache.

Read only POTCAR files under the exact historical Ca job work directory and
known reference directories. Do not read POSCAR/CONTCAR/OUTCAR, do not export
potential contents, and never download a potential from the internet.
"""
import hashlib, json, os, re, subprocess, sys, time
from pathlib import Path
HOME=Path.home()
LABELS={'Mg':'Mg','Ca':'Ca_sv','O':'O','H':'H'}
CACHE=HOME/'VASP/potential_sets/Ca_corrosion_fresh_20261014'

def blocks(p):
    if not p.is_file() or p.stat().st_size>6000000:return {}
    chunks=[];buf=[]
    for line in p.read_bytes().splitlines(keepends=True):
        buf.append(line)
        if line.strip()==b'End of Dataset':chunks.append(b''.join(buf));buf=[]
    d={}
    for b in chunks:
        t=re.findall(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)',b)
        if len(t)!=1:continue
        label=t[0].decode()
        if label not in LABELS.values():continue
        en=re.findall(rb'ENMAX\s*=\s*([\d.]+)',b)
        if not en or 1.3*float(en[0])>520.00001:continue
        d[label]=(b,float(en[0]))
    return d

def main():
    if len(sys.argv)!=2:raise SystemExit('Usage: prepare_potentials.py path/to/fresh.py')
    fresh=Path(sys.argv[1]).resolve()
    if fresh.name!='fresh.py' or not fresh.is_file():raise RuntimeError('Expected reviewed fresh.py')
    roots=[]
    # Known exact project hint and the explicitly referenced Ca job only.
    q=HOME/'next_stage_63131486_20261002'
    if q.is_dir():roots.append(q)
    r=subprocess.run(['sacct','-X','-n','-P','-j','63131486','-o','WorkDir%300'],capture_output=True,text=True,timeout=25)
    if r.returncode==0:
        for line in r.stdout.splitlines():
            s=line.split('|',1)[0].strip()
            if not s:continue
            p=Path(s)
            try:relative=p.relative_to(HOME)
            except ValueError:continue
            if not relative.parts or not p.is_dir():continue
            roots.append(p)
            # At most one parent, never the home directory itself.
            if len(relative.parts)>=2:roots.append(p.parent)
    direct=[HOME/'VASP/potcar_profiles/Mg/POTCAR',
            HOME/'Mg2Ni_Mg_H_thermo_validation_20261010/A-H-Ni_Mg_Mg2Ni/01_relax/POTCAR']
    found={};provenance={}
    def read(p):
        for label,(b,en) in blocks(p).items():
            if label not in found:
                found[label]=b
                provenance[label]={'sha256':hashlib.sha256(b).hexdigest(),'ENMAX_eV':en,'source_path_local_only':str(p)}
    for p in direct:read(p)
    visited=0;files=0;start=time.monotonic()
    for root in dict.fromkeys(roots):
        for base,dirs,names in os.walk(root,followlinks=False):
            visited+=1;p=Path(base);depth=len(p.relative_to(root).parts)
            dirs[:]=[] if depth>=5 else [x for x in dirs if not x.startswith('.') and x not in ('__pycache__','node_modules')]
            if 'POTCAR' in names:read(p/'POTCAR');files+=1
            if set(found)>=set(LABELS.values()) or visited>=600 or files>=100 or time.monotonic()-start>20:break
        if set(found)>=set(LABELS.values()) or visited>=600 or files>=100 or time.monotonic()-start>20:break
    print('POTENTIAL_METADATA_ONLY|'+json.dumps({'labels_found':list(found),'files_checked':files,'directories_checked':visited}),flush=True)
    missing=set(LABELS.values())-set(found)
    if missing:raise RuntimeError('Missing licensed datasets '+','.join(sorted(missing))+'. No model or Slurm job was created; no old geometry was read.')
    CACHE.mkdir(parents=True,exist_ok=True)
    for label,b in found.items():
        p=CACHE/label;p.mkdir(exist_ok=True);f=p/'POTCAR'
        if f.exists() and f.read_bytes()!=b:raise RuntimeError('Existing cache differs; do not overwrite '+label)
        if not f.exists():f.write_bytes(b);f.chmod(0o600)
    (CACHE/'provenance.json').write_text(json.dumps(provenance,indent=2));(CACHE/'provenance.json').chmod(0o600)
    env=dict(os.environ);env['CA_PBE_POTENTIAL_ROOT']=str(CACHE)
    print('LICENSED_POTENTIALS_READY|Mg,Ca_sv,O,H|no potential exported',flush=True)
    p=subprocess.run([sys.executable,str(fresh),'--submit'],env=env)
    raise SystemExit(p.returncode)
if __name__=='__main__':main()
