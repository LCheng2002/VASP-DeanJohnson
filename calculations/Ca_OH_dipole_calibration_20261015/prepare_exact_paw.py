#!/usr/bin/env python3
"""Locate exact original PAW datasets in a bounded set of authorized Ca paths.
Only inspect POTCAR datasets locally. Never export their bytes or change old files.
If any exact hash is missing, do not call the calibration submission driver.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

HOME = Path.home()
HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE/'plan.json').read_text())
EXPECTED = PLAN['exact_PAW']
CACHE = HOME/'VASP/potential_sets/Ca_OH_exact_20261015'


def main():
    found = {}; origins = {}; seen = set(); candidates = []; visited = 0
    def read(path):
        path = Path(path)
        if path in seen or not path.is_file() or path.stat().st_size > 6000000:
            return
        seen.add(path); buffer = []
        for line in path.read_bytes().splitlines(keepends=True):
            buffer.append(line)
            if line.strip() != b'End of Dataset': continue
            block = b''.join(buffer); buffer = []
            checksum = hashlib.sha256(block).hexdigest()
            labels = re.findall(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)', block)
            if len(labels) != 1: continue
            label = labels[0].decode()
            if label == 'O':
                candidates.append({'label':label,'sha256':checksum,'relative_path':str(path.relative_to(HOME))})
            for element, spec in EXPECTED.items():
                if checksum == spec['sha256'] and label == spec['label']:
                    if element not in found:
                        found[element] = block
                        origins[element] = str(path.relative_to(HOME))
    for label in ('Mg','Ca_sv','O','H'):
        read(CACHE/label/'POTCAR')
        read(HOME/'VASP/potential_sets/Ca_corrosion_fresh_20261014'/label/'POTCAR')
    roots = []
    for job in ('63131486','63306398'):
        p = subprocess.run(['sacct','-X','-n','-P','-j',job,'-o','WorkDir%300'],
                           capture_output=True,text=True,timeout=25)
        if p.returncode:
            print('OLD_JOB_LOCATION_QUERY_FAILED|'+job+'|'+str(p.returncode),flush=True)
            continue
        for row in p.stdout.splitlines():
            text = row.split('|')[0].strip()
            if not text or text in ('Unknown','None'): continue
            q = Path(text)
            try: relative = q.relative_to(HOME)
            except ValueError: continue
            if not relative.parts or not q.is_dir() or q.is_symlink(): continue
            roots.append(q)
            if len(relative.parts) >= 2 and not q.parent.is_symlink(): roots.append(q.parent)
    for rel in PLAN['optional_historical_potcar_roots']:
        p = HOME/rel
        if p.is_dir() and not p.is_symlink(): roots.append(p)
    roots = list(dict.fromkeys(roots))
    print('EXACT_PAW_SEARCH_ROOTS|'+json.dumps([str(p.relative_to(HOME)) for p in roots]),flush=True)
    start = time.monotonic()
    for root in roots:
        if len(found) == 4: break
        for base, dirs, files in os.walk(root,followlinks=False):
            visited += 1; p = Path(base); depth = len(p.relative_to(root).parts)
            dirs[:] = [] if depth >= 6 else [n for n in dirs if not n.startswith('.') and n not in ('node_modules','__pycache__')]
            if 'POTCAR' in files: read(p/'POTCAR')
            if len(found)==4 or visited>=1500 or len(seen)>=100 or time.monotonic()-start>30: break
        if visited>=1500 or len(seen)>=100 or time.monotonic()-start>30: break
    print('EXACT_PAW_SEARCH_RESULT|'+json.dumps({'matched':sorted(found),'files':len(seen),'directories':visited,'oxygen_candidates':candidates[:25]}),flush=True)
    if set(found) != set(EXPECTED):
        raise RuntimeError('Original PAW not located: '+','.join(sorted(set(EXPECTED)-set(found)))+'. No sbatch executed.')
    # Do not overwrite the independent cache if it was created by another run.
    CACHE.mkdir(parents=True,exist_ok=True)
    for element, block in found.items():
        path = CACHE/EXPECTED[element]['label']; path.mkdir(exist_ok=True); f=path/'POTCAR'
        if f.exists() and hashlib.sha256(f.read_bytes()).hexdigest()!=EXPECTED[element]['sha256']:
            raise RuntimeError('Existing exact-calibration cache differs; stop without overwriting')
        if not f.exists():
            f.write_bytes(block); f.chmod(0o600)
    provenance={'data_hashes':{e:EXPECTED[e]['sha256'] for e in found},'source_paths_local_only':origins}
    p=CACHE/'provenance.json';p.write_text(json.dumps(provenance,indent=2));p.chmod(0o600)
    print('EXACT_ORIGINAL_PAW_READY|Mg,Ca_sv,O,H|no original files changed',flush=True)
    if PLAN['potcar_cache_relpath']!='VASP/potential_sets/Ca_OH_exact_20261015':
        raise RuntimeError('Calibration plan does not point at the separate exact cache')
    result=subprocess.run([sys.executable,str(HERE/'calibrate.py'),'--submit'])
    raise SystemExit(result.returncode)

if __name__=='__main__':
    main()
