#!/usr/bin/env python3
"""Generate 0-based half-open BED intervals from frozen MosqEdit-R site manifest.
No networking; stdlib only. Refuses malformed coordinates/identifiers.
"""
from pathlib import Path
import csv, collections, hashlib
ROOT = Path(__file__).resolve().parent
manifest = ROOT / '13_population_site_manifest_FROZEN.csv'
out = ROOT / 'regions'
out.mkdir(exist_ok=True)
chroms = ('2L','2R','3L','3R','X')
rows = list(csv.DictReader(manifest.open(newline='', encoding='utf-8-sig')))
if len(rows) != 530 or len({r['population_site_id'] for r in rows}) != len(rows):
    raise SystemExit('STOP: Frozen panel not 530 unique sites. Do not proceed.')
coords = collections.defaultdict(list)
for r in rows:
    chrom = r['genomic_seqid']; s=int(r['genomic_start']); e=int(r['genomic_end'])
    if s < 1 or e < s or (e-s+1)!=23:
        raise SystemExit(f'STOP: Unexpected site length or coordinates: {r["population_site_id"]}')
    if chrom in chroms:
        coords[chrom].append((s-1,e))  # BED 0 based, half-open
    elif chrom != 'UNKN':
        raise SystemExit(f'STOP: Unknown chromosome identifier: {chrom}')
with (out/'UNKN_not_automatically_queried.tsv').open('w',newline='') as fh:
    wr=csv.writer(fh,delimiter='\t'); wr.writerow(['population_site_id','gene_id','seqid','start_1based','end_1based'])
    for r in rows:
        if r['genomic_seqid']=='UNKN': wr.writerow([r['population_site_id'],r['gene_id'],r['genomic_seqid'],r['genomic_start'],r['genomic_end']])
for c in chroms:
    intervals=sorted(coords[c]); merged=[]
    for s,e in intervals:
        if merged and s <= merged[-1][1]: merged[-1] = (merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    with (out/f'{c}_merged_0based.bed').open('w',newline='') as fh:
        for s,e in merged: fh.write(f'{c}\t{s}\t{e}\n')
    with (out/f'{c}_individual_0based.bed').open('w',newline='') as fh:
        for s,e in intervals: fh.write(f'{c}\t{s}\t{e}\n')
    print(f'{c}: {len(intervals)} frozen sites; {len(merged)} merged query intervals')
print('UNKN:',sum(r['genomic_seqid']=='UNKN' for r in rows),'sites excluded pending coordinate mapping')
print('Manifest SHA256:',hashlib.sha256(manifest.read_bytes()).hexdigest())
