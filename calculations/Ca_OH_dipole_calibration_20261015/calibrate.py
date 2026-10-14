#!/usr/bin/env python3
"""Four single-point dipole calibrations on the user's verified OH/H2O slabs.
No geometry optimization, new adsorbates, automatic retry, or CaFresh restart.
POTCAR datasets stay inside the licensed cluster account. Standard library only.
"""
import argparse
import datetime
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE / 'plan.json').read_text())
ROOT = Path.home() / PLAN['batch']
SPECS = {m['name']: m for m in PLAN['models']}
NAMES = list(SPECS)
FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
SHIFTS = tuple(itertools.product((-1, 0, 1), repeat=3))


def require(ok, message):
    if not ok:
        raise RuntimeError('CHECK_FAILED|' + message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temp.replace(path)


def now():
    return datetime.datetime.now().astimezone().isoformat()


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def cart(f, cell):
    return [sum(f[k] * cell[k][j] for k in range(3)) for j in range(3)]


def fractional(x, cell):
    a, b, c = cell
    det = dot(a, cross(b, c))
    require(abs(det) > 1e-9, 'Singular cell')
    return [dot(x, cross(b,c))/det, dot(x, cross(c,a))/det, dot(x, cross(a,b))/det]


def mic_length(x, y, cell):
    f = fractional([a-b for a,b in zip(x,y)], cell)
    nearest = [round(v) for v in f]
    return min(math.sqrt(dot(v, v)) for shift in SHIFTS
               for v in [cart([f[k]-nearest[k]-shift[k] for k in range(3)], cell)])


def parse_structure(raw):
    lines = raw.splitlines()
    require(len(lines) > 9, 'Incomplete POSCAR/CONTCAR')
    scale = float(lines[1]); require(scale > 0, 'Expected positive single lattice scale')
    cell = [[float(v)*scale for v in line.split()[:3]] for line in lines[2:5]]
    elements = lines[5].split(); counts = list(map(int, lines[6].split()))
    require(len(elements) == len(counts), 'Invalid species roster')
    species = sum(([e]*n for e,n in zip(elements,counts)), [])
    k = 7; selective = lines[k].strip().lower().startswith('s'); k += int(selective)
    direct = lines[k].strip().lower().startswith('d'); k += 1
    xyz, masks = [], []
    for line in lines[k:k+len(species)]:
        cols = line.split(); vals = list(map(float, cols[:3]))
        require(len(vals) == 3, 'Incomplete atomic coordinate')
        xyz.append(cart(vals, cell) if direct else [v*scale for v in vals])
        mask = [v.upper().startswith('T') for v in cols[3:6]] if selective else [True]*3
        require(len(mask) == 3, 'Missing selective dynamics mask')
        masks.append(mask)
    require(len(xyz) == len(species), 'Atomic count mismatch')
    require(all(math.isfinite(v) for p in xyz for v in p), 'Nonfinite geometry')
    geometry_only = '\n'.join(lines[:k+len(species)]) + '\n'
    return {'cell': cell, 'species': species, 'xyz': xyz, 'masks': masks,
            'geometry_text': geometry_only, 'composition': dict(zip(elements, counts))}


def dipole_spec(s):
    c = s['cell']
    require(max(abs(c[0][2]),abs(c[1][2]),abs(c[2][0]),abs(c[2][1])) < 1e-8,
            'Calibration assumes z normal to the basal plane')
    length = c[2][2]; require(length > 25, 'Unexpected vacuum-cell length')
    f = [fractional(x, c) for x in s['xyz']]
    zs = sorted(v[2] % 1 for v in f)
    gaps = [zs[i+1]-zs[i] for i in range(len(zs)-1)] + [zs[0]+1-zs[-1]]
    i = max(range(len(gaps)), key=gaps.__getitem__)
    origin = zs[(i+1) % len(zs)]
    masses = [PLAN['atomic_masses_for_DIPOL'][e] for e in s['species']]
    total = sum(masses)
    centre = [sum(m*v[k] for m,v in zip(masses,f))/total for k in (0,1)]
    centre += [origin + sum(m*((v[2]-origin) % 1) for m,v in zip(masses,f))/total]
    centre = [v % 1 for v in centre]
    cut = (centre[2] + .5) % 1
    clearance = min(abs(((v[2]-cut+.5) % 1)-.5)*length for v in f)
    require(clearance >= PLAN['minimum_atomic_clearance_to_dipole_cut_A'],
            'Dipole discontinuity candidate is too close to atoms')
    return {'DIPOL_direct': centre, 'cut_z_fractional': cut,
            'minimum_atom_to_cut_A': clearance, 'largest_atom_free_gap_A': max(gaps)*length,
            'definition': 'Same mass-centre rule for all four states; direct lattice coordinates'}


def parse_incar(raw):
    result = {}
    for line in raw.splitlines():
        for part in re.split(r';', re.split(r'[!#]', line)[0]):
            if '=' not in part: continue
            key, val = (x.strip() for x in part.split('=',1)); key = key.upper()
            require(key not in result or result[key] == val, 'Conflicting duplicate INCAR '+key)
            result[key] = val
    return result


def verify_sources(source):
    require(digest(source/'AUDIT.json') == PLAN['source_AUDIT_sha256'], 'AUDIT checksum differs')
    audit = json.loads((source/'AUDIT.json').read_text())
    verified = {}
    common_cell = None
    for name in NAMES:
        spec = SPECS[name]; d = source/'models'/name
        files = {'CONTCAR': d/'results/01_relax/CONTCAR',
                 'INCAR': d/'calculation_inputs/INCAR.static',
                 'KPOINTS': d/'calculation_inputs/KPOINTS'}
        for kind, path in files.items():
            require(path.is_file() and digest(path) == spec['source_'+kind+'_sha256'],
                    'Source hash mismatch '+name+'/'+kind)
        a = audit['models'][name]
        for stage in ('01_relax','02_static'):
            v = a[stage]
            require(v['normal_end'] and v['electronic_ediff_last'], 'Unverified archived SCF '+name)
            if stage == '01_relax': require(v['ionic_convergence_marker'], 'Unconverged archived parent')
        require(abs(a['02_static']['F_eV']-spec['old_F_eV']) < 1e-8, 'Wrong reference energy')
        s = parse_structure(files['CONTCAR'].read_text())
        require(s['composition'] == spec['composition'] and len(s['species']) == spec['expected_atoms'],
                'Wrong source composition '+name)
        require(sum(not any(m) for m in s['masks']) == 32, 'Expected 32 fixed metal atoms')
        require(all(e in ('Mg','Ca') for e,m in zip(s['species'],s['masks']) if not any(m)),
                'Adsorbate incorrectly frozen')
        if common_cell is None: common_cell = s['cell']
        require(max(abs(a-b) for r,t in zip(common_cell,s['cell']) for a,b in zip(r,t)) < 1e-10,
                'The four old cells differ')
        old = parse_incar(files['INCAR'].read_text())
        numeric = {'ENCUT':520.,'EDIFF':1e-7,'ISMEAR':1.,'SIGMA':.2,'NELM':240.,
                   'NCORE':8.,'KPAR':1.,'ISYM':0.,'NSW':0.,'IBRION':-1.,'ISIF':2.}
        for key, val in numeric.items(): require(float(old[key]) == val, 'Old setting mismatch '+key)
        require(old['GGA']=='PE' and old['LREAL'].upper()=='.FALSE.', 'Different PBE/LREAL baseline')
        require(int(old.get('ISPIN','1'))==1, 'Different spin baseline')
        require(old.get('LDAU','.FALSE.').upper() in ('.FALSE.','F'), 'Unexpected Hubbard correction')
        new = dict(old)
        dipol = dipole_spec(s)
        new.update({'SYSTEM': 'CaP0_'+name, 'ISPIN':'1', 'LDAU':'.FALSE.',
                    'LDIPOL':'.TRUE.', 'IDIPOL':'3',
                    'DIPOL':' '.join('%.12f'%v for v in dipol['DIPOL_direct']),
                    'NBANDS':str(spec['NBANDS']), 'ISTART':'0', 'ICHARG':'2',
                    'LWAVE':'.FALSE.', 'LCHARG':'.TRUE.', 'LVTOT':'.TRUE.'})
        verified[name] = {'structure':s, 'incar':'\n'.join(k+' = '+v for k,v in new.items())+'\n',
                          'kpoints':files['KPOINTS'].read_text(), 'dipole':dipol,
                          'original_ldipol':old.get('LDIPOL','.FALSE.'),
                          'original_dipol':old.get('DIPOL')}
    return verified


def matching_potentials():
    """Use exact block hashes, not just element names; read at most named Ca roots."""
    home = Path.home(); found = {}; paths_seen = set()
    def read(path):
        if path in paths_seen or not path.is_file() or path.stat().st_size > 6000000: return
        paths_seen.add(path)
        buffer = []
        for line in path.read_bytes().splitlines(keepends=True):
            buffer.append(line)
            if line.strip() != b'End of Dataset': continue
            data = b''.join(buffer); buffer=[]; h=hashlib.sha256(data).hexdigest()
            for e, setting in PLAN['exact_PAW'].items():
                if h == setting['sha256']: found[e] = data
    cache = home/PLAN['potcar_cache_relpath']
    for setting in PLAN['exact_PAW'].values(): read(cache/setting['label']/'POTCAR')
    start = time.monotonic(); visited = 0
    for relative in PLAN['optional_historical_potcar_roots']:
        if len(found)==4: break
        root = home/relative
        if not root.is_dir() or root.is_symlink(): continue
        for base, dirs, filenames in os.walk(root, followlinks=False):
            visited+=1; d=Path(base); depth=len(d.relative_to(root).parts)
            dirs[:] = [] if depth>=6 else [n for n in dirs if not n.startswith('.')]
            if 'POTCAR' in filenames: read(d/'POTCAR')
            if len(found)==4 or visited>=400 or len(paths_seen)>=60 or time.monotonic()-start>20: break
        if visited>=400 or len(paths_seen)>=60 or time.monotonic()-start>20: break
    print('EXACT_POTENTIAL_METADATA|'+json.dumps({'matched_elements':sorted(found),'files_read':len(paths_seen)}),flush=True)
    require(set(found)==set(PLAN['exact_PAW']), 'Original PAW datasets unavailable; no job submitted')
    return found


def prepare(source, preview=False):
    require(not ROOT.exists(), 'Target directory exists; never overwrite a previous attempt')
    models = verify_sources(source)
    blocks = {} if preview else matching_potentials()
    if not preview:
        for name,m in models.items():
            els = list(dict.fromkeys(m['structure']['species']))
            p = b''.join(blocks[e] for e in els)
            require(hashlib.sha256(p).hexdigest()==SPECS[name]['source_POTCAR_sha256'],
                    'Combined original POTCAR mismatch '+name)
    ROOT.mkdir(parents=True); (ROOT/'logs').mkdir()
    for name,m in models.items():
        d=ROOT/name; d.mkdir(); s=m['structure']
        (d/'POSCAR').write_text(s['geometry_text'])
        (d/'INCAR').write_text(m['incar']); (d/'KPOINTS').write_text(m['kpoints'])
        if not preview:
            data=b''.join(blocks[e] for e in dict.fromkeys(s['species']))
            (d/'POTCAR').write_bytes(data); (d/'POTCAR').chmod(0o600)
        model={**SPECS[name], 'dipole':m['dipole'], 'original_ldipol':m['original_ldipol'],
               'original_dipol':m['original_dipol'], 'fixed_atoms':32,
               'new_geometry_equals_source':True, 'mode':'static calibration, not a relaxed structure'}
        save(d/'model.json',model)
        print('MODEL_PREPARED|'+name+'|DIPOL='+m['incar'].split('DIPOL = ')[-1].split('\n')[0]+
              '|atoms='+str(len(s['species']))+'|cut_clearance_A='+str(m['dipole']['minimum_atom_to_cut_A']),flush=True)
    shutil.copy2(HERE/'plan.json',ROOT/'plan.json'); shutil.copy2(Path(__file__),ROOT/'calibrate.py')
    job='''#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load oneapi vasp/6.3.2-intel-2021.4.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited
ulimit -l unlimited
python3 calibrate.py --run-index "$SLURM_ARRAY_TASK_ID"
'''
    (ROOT/'calibration.slurm').write_text(job)
    subprocess.run(['bash','-n',str(ROOT/'calibration.slurm')],check=True)
    hashes={str(p.relative_to(ROOT)):digest(p) for p in ROOT.rglob('*') if p.is_file()}
    save(ROOT/'manifest.json',{'created_at':now(),'preview':preview,'models':NAMES,
                              'input_hashes':hashes,'max_models':4,'no_CaFresh_restart':True})


def read_output(path, structure):
    text=Path(path).read_text(errors='replace')
    last = lambda pat: (re.findall(pat,text)[-1] if re.findall(pat,text) else None)
    f=last(r'free\s+energy\s+TOTEN\s*=\s*('+FLOAT+')')
    e0=last(r'energy\(sigma->0\)\s*=\s*('+FLOAT+')')
    iterations=list(re.finditer(r'Iteration\s+(\d+)\(\s*(\d+)\)',text))
    first_last_cycle=next((m.start() for m in iterations if m.group(1)==iterations[-1].group(1)),0) if iterations else 0
    scf='aborting loop because EDIFF is reached' in text[first_last_cycle:]
    part=text.rsplit('TOTAL-FORCE (eV/Angst)',1)
    coords=[]; forces=[]
    if len(part)==2:
        for line in part[-1].splitlines():
            cols=line.split()
            if len(cols)==6:
                try:
                    vals=list(map(float,cols)); coords.append(vals[:3]); forces.append(vals[3:]); continue
                except ValueError: pass
            if forces: break
    require(len(forces)==len(structure['species']), 'Final force table missing or incomplete')
    maxfree=max(math.sqrt(sum(v*v for v,free in zip(force,mask) if free))
                for force,mask in zip(forces,structure['masks']))
    displacement=max(mic_length(a,b,structure['cell']) for a,b in zip(coords,structure['xyz']))
    require(displacement<2e-4,'Static geometry differs from input')
    parameters={}
    for key in ('ISPIN','NBANDS','NELECT','IDIPOL','ENCUT','ISMEAR','SIGMA'):
        val=last(r'\b'+key+r'\s*=\s*('+FLOAT+')')
        parameters[key]=float(val.replace('D','E')) if val else None
    parameters['LDIPOL']=last(r'\bLDIPOL\s*=\s*([TF])')
    return {'F_eV':float(f) if f else None,'E0_eV':float(e0) if e0 else None,
            'last_scf_ediff':scf,'normal_end':'General timing and accounting informations' in text,
            'max_free_force_eV_A':maxfree,'final_force_atoms':len(forces),
            'max_geometry_difference_A':displacement,'effective_parameters':parameters}


def planar_grid(path):
    """Return first-grid planar averages only; do not export raw charge or potential files."""
    with Path(path).open(errors='replace') as f:
        header=[next(f) for _ in range(7)]
        counts=list(map(int,header[6].split())); n=sum(counts)
        line=next(f)
        if line.strip().lower().startswith('s'): line=next(f)
        require(line.strip().lower().startswith(('d','c','k')), 'Invalid grid header')
        for _ in range(n): next(f)
        for line in f:
            cols=line.split()
            if not cols: continue
            require(len(cols)==3 and all(v.isdigit() for v in cols), 'Grid dimension line missing')
            nx,ny,nz=map(int,cols); break
        else: raise RuntimeError('No grid found')
        require(nx>0 and ny>0 and nz>0 and nx*ny*nz<200000000,'Unexpected grid dimensions')
        plane=nx*ny; total=plane*nz; means=[0.0]*nz; k=0
        for line in f:
            for val in line.split():
                if k==total: break
                means[k//plane]+=float(val.replace('D','E'))/plane; k+=1
            if k==total: break
        require(k==total,'Grid incomplete')
    return (nx,ny,nz),means


def run(index):
    require(0<=index<4,'Invalid model index')
    require(os.environ.get('SLURM_NTASKS')=='64','VASP must run inside prescribed 64-rank allocation')
    name=NAMES[index]; d=ROOT/name
    if (d/'CALIBRATION_COMPLETE').exists(): print('ALREADY_CALIBRATED|'+name); return
    manifest=json.loads((ROOT/'manifest.json').read_text())
    require(not manifest['preview'],'Cannot run a preview')
    for file in ('POSCAR','INCAR','KPOINTS','POTCAR','model.json'):
        require(digest(d/file)==manifest['input_hashes'][name+'/'+file],'Pinned input changed '+name+'/'+file)
    require(digest(ROOT/'calibrate.py')==manifest['input_hashes']['calibrate.py'],'Production driver changed')
    (d/'.run_started').mkdir(exist_ok=False)
    with (d/'vasp.stdout').open('w') as out,(d/'vasp.stderr').open('w') as err:
        proc=subprocess.run(['mpirun','-np','64','vasp_std'],cwd=d,stdout=out,stderr=err)
    require(proc.returncode==0,'VASP returned nonzero for '+name)
    s=parse_structure((d/'POSCAR').read_text()); v=read_output(d/'OUTCAR',s)
    expected={'ISPIN':1.,'NBANDS':float(SPECS[name]['NBANDS']),'NELECT':SPECS[name]['NELECT'],
              'IDIPOL':3.,'ENCUT':520.,'ISMEAR':1.,'SIGMA':.2}
    p=v['effective_parameters']
    require(all(p.get(k) is not None and abs(p[k]-x)<1e-6 for k,x in expected.items()) and p['LDIPOL']=='T',
            'Effective VASP settings do not match calibrated plan')
    require(v['last_scf_ediff'] and v['normal_end'] and all(v[k] is not None and math.isfinite(v[k]) for k in ('F_eV','E0_eV')),
            'Electronic single-point calculation not converged')
    meta=json.loads((d/'model.json').read_text())
    gridinfo={}; profiles={}
    for fn,key in (('CHGCAR','charge'),('LOCPOT','potential')):
        shape,avg=planar_grid(d/fn); profiles[key]=avg
        gridinfo[key]={'shape':shape,'mean':sum(avg)/len(avg)}
        (d/(key+'_planar.tsv')).write_text('z_fractional\tplanar_average\n'+'\n'.join('%.10f\t%.12g'%(i/len(avg),x) for i,x in enumerate(avg))+'\n')
    rho=profiles['charge']; n=len(rho); cut=meta['dipole']['cut_z_fractional']; j=round(cut*n)%n
    ratio=max(abs(rho[(j+k)%n]) for k in (-1,0,1))/max(max(abs(x) for x in rho),1e-20)
    require(abs(gridinfo['charge']['mean']-SPECS[name]['NELECT']) < .01,'CHGCAR electron-count normalization mismatch')
    v.update({'model':name,'electronic_calibration_verified':True,
              'geometry_force_gate_passed':v['max_free_force_eV_A']<.02,
              'needs_consistent_relaxation':v['max_free_force_eV_A']>=.02,
              'old_F_eV':SPECS[name]['old_F_eV'],'old_E0_eV':SPECS[name]['old_E0_eV'],
              'delta_F_eV':v['F_eV']-SPECS[name]['old_F_eV'],
              'delta_E0_eV':v['E0_eV']-SPECS[name]['old_E0_eV'],
              'dipole':meta['dipole'],'grid_metadata':gridinfo,'density_ratio_near_cut_to_peak':ratio,
              'vacuum_cut_density_small':ratio<1e-4,
              'additional_dipole_center_or_vacuum_sensitivity_may_be_required':True,
              'no_relaxation_or_followup_auto_submitted':True})
    save(d/'validation.json',v); (d/'CALIBRATION_COMPLETE').touch()
    print('CALIBRATION_RESULT|'+json.dumps(v),flush=True)
    print('NEXT_GATE|'+name+'|'+('SHORT_RELAXATION_NEEDED' if v['needs_consistent_relaxation'] else 'FORCE_GATE_PASSED'),flush=True)


def submit(source):
    require(shutil.which('sbatch') is not None,'Slurm not available')
    if (ROOT/'submission.json').exists():
        print('EXISTING_SUBMISSION|'+(ROOT/'submission.json').read_text()); return
    if not ROOT.exists(): prepare(source,False)
    m=json.loads((ROOT/'manifest.json').read_text())
    require(not m['preview'] and m['models']==NAMES and m['max_models']==4,'Invalid production manifest')
    for rel,h in m['input_hashes'].items(): require(digest(ROOT/rel)==h,'Changed production file '+rel)
    (ROOT/'.submission_attempted').mkdir(exist_ok=False)
    p=subprocess.run(['sbatch','--parsable','--job-name=CaOH_P0_calibration','calibration.slurm'],cwd=ROOT,
                     capture_output=True,text=True,timeout=45)
    require(p.returncode==0,'sbatch rejected: '+p.stderr[:400])
    jid=p.stdout.strip().split(';')[0]; require(re.fullmatch(r'\d+',jid),'Uncertain submission response; do not retry blindly')
    rec={'job_id':jid,'array':'0-3','models':NAMES,'submitted_at':now(),'four_statics_only':True,
         'cancelled_CaFresh_restarted':False,'automatic_followups':False}
    save(ROOT/'submission.json',rec); print('SUBMITTED|CaOH_P0|'+jid+'|array=0-3',flush=True)


def status():
    print('CHECKED_AT|'+now())
    if not (ROOT/'submission.json').exists(): print('NOT_SUBMITTED'); return
    rec=json.loads((ROOT/'submission.json').read_text()); print('SUBMISSION|'+json.dumps(rec))
    for cmd in (['squeue','-r','-h','-j',rec['job_id'],'-o','%i|%j|%T|%M|%R'],
                ['sacct','-X','-n','-P','-j',rec['job_id'],'-o','JobID,State,ExitCode,Elapsed,Start,End']):
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
        print('SLURM_CHECK_RC|'+str(p.returncode)); print(p.stdout.strip())
        if p.returncode: print('SLURM_READ_ERROR|'+p.stderr[:250])
    vals={}
    for name in NAMES:
        p=ROOT/name/'validation.json'
        if p.exists() and (ROOT/name/'CALIBRATION_COMPLETE').exists(): vals[name]=json.loads(p.read_text())
        print('MODEL|'+name+'|'+json.dumps(vals.get(name,'NOT_YET_CALIBRATED')))
    if len(vals)==4:
        dd=lambda key: vals['M1_OH_H2O'][key]-vals['M1_OH'][key]-vals['M0_OH_H2O'][key]+vals['M0_OH'][key]
        summary={'new_double_difference_F_eV':dd('F_eV'),'new_double_difference_E0_eV':dd('E0_eV'),
                 'old_double_difference_F_eV':-.21814632,'old_double_difference_E0_eV':-.21739583,
                 'all_force_gates_passed':all(v['geometry_force_gate_passed'] for v in vals.values()),
                 'all_cut_densities_small':all(v['vacuum_cut_density_small'] for v in vals.values()),
                 'not_an_absolute_adsorption_energy_or_corrosion_rate':True}
        print('FOUR_POINT_COMPARISON|'+json.dumps(summary))


def main():
    global ROOT
    ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--preview',action='store_true'); g.add_argument('--submit',action='store_true')
    g.add_argument('--run-index',type=int); g.add_argument('--status',action='store_true')
    ap.add_argument('--source-root',type=Path); ap.add_argument('--preview-root',type=Path)
    args=ap.parse_args()
    if args.preview_root:
        require(args.preview,'Custom output root is preview-only'); ROOT=args.preview_root.resolve()
    source=args.source_root or HERE.parents[1]/PLAN['source_repo_directory']
    if args.preview: prepare(source,True)
    elif args.submit: submit(source)
    elif args.status: status()
    else: run(args.run_index)


if __name__=='__main__':
    main()
