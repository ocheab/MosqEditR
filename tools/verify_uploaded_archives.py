#!/usr/bin/env python3
"""Read-only source archiving preflight. Does NOT execute the analysis or publish files."""
import hashlib,json,sys
from pathlib import Path
from zipfile import ZipFile

def digest(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for part in iter(lambda:f.read(8*1024*1024),b''):h.update(part)
 return h.hexdigest()

if len(sys.argv)!=4:
 sys.exit('Usage: python tools/verify_uploaded_archives.py PACKAGE.zip PHASE2.zip REVDOCS.zip')
manifest=json.loads((Path(__file__).resolve().parents[1]/'analysis/provenance/source_archive_sha256.json').read_text())
for label,path in zip(manifest,sys.argv[1:]):
 p=Path(path)
 if not p.is_file():sys.exit('Missing: '+str(p))
 got=digest(p)
 expected=manifest[label]['sha256']
 print(f'{label}: {"PASS" if got==expected else "FAIL"} sha256={got}')
 if got!=expected:sys.exit('Archive mismatch '+str(p))
 with ZipFile(p) as z: bad=z.testzip()
 if bad:sys.exit('Bad ZIP member: '+bad)
print('PASS: byte identity of all three archived input ZIPs verified; analysis reexecution NOT tested.')
