#!/usr/bin/env python3
"""Simple subset checks; absence of variants is NOT proof of reference conservation."""
import gzip, sys, csv, pathlib, collections
if len(sys.argv)!=2:
    raise SystemExit('Usage: python 03_check_subset.py output/PILOT_2R.vcf.gz')
p=pathlib.Path(sys.argv[1]); samples=[]; variants=0; chromosome=collections.Counter(); missing_gt=0
with gzip.open(p,'rt') as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            h=line.rstrip('\r\n').split('\t'); samples=h[9:]; continue
        if not line.strip(): continue
        col=line.rstrip('\r\n').split('\t'); variants+=1; chromosome[col[0]]+=1
        if len(col)!=len(samples)+9: raise SystemExit('STOP: VCF row has unexpected number of columns')
print('File:',p)
print('Samples:',len(samples))
print('Variant rows:',variants)
print('Reference contigs:',dict(chromosome))
if not samples: raise SystemExit('STOP: No individual genotype columns present in this subset!')
metadata=next((p for p in (pathlib.Path('downloads/samples.meta.txt'), pathlib.Path('../downloads/samples.meta.txt')) if p.exists()), None)
if metadata is not None:
    with metadata.open() as fh:
        meta=csv.DictReader(fh,delimiter='\t')
        if meta.fieldnames and 'ox_code' in meta.fieldnames:
            known={r['ox_code'] for r in meta if r.get('ox_code')}
            print('Samples matched in metadata:',sum(x in known for x in samples), '/',len(samples))
        else: print('Metadata format: inspect manually (may not be tab-delimited).')
else: print('Metadata file not in kit. Original download is stored in parent project downloads/ folder.')
