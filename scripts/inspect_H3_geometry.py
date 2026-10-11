#!/usr/bin/env python3
"""Read one explicitly identified relaxed atom structure; no jobs submitted."""
import json, hashlib
from pathlib import Path
p=Path.home()/'TaskA_Mg17La2_H_enrichment_20261010/H3_addLa/01_relax/CONTCAR'
b=p.read_bytes()
if len(b)>20000:raise RuntimeError('Unexpected structure size')
assert hashlib.sha256(b).hexdigest()=='7fdf86b6153683dda3ec22231ddcb039db5a48fb6bda0caf5a63308ce66e3b24','Parent changed'
s=b.decode('utf-8').splitlines()
assert s[5].split()==['Mg','La','H'] and list(map(int,s[6].split()))==[68,8,3]
i=7
if s[i].lower().startswith('s'):i+=1
i+=1
print(json.dumps({'parent_sha256':hashlib.sha256(b).hexdigest(),'structure':'\n'.join(s[:i+79])+'\n'}))
