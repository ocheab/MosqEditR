#!/usr/bin/env python3
"""MosqEdit-R frozen-guide 23-bp FASTA concordance QC (not a conservation test).

Inputs relative to current extraction-kit directory:
  13_population_site_manifest_FROZEN.csv
  output/AR1_PASS_SNP_REF_5arms.tsv
  ../downloads/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz

Uses only Python standard library. One FASTA chromosome in memory at a time.
Compares every base of 515 chromosome-arm guides after strand transformation,
and the 2,006 PASS VCF REF alleles against the AgamP4 FASTA.
The 15 UNKN sites are explicitly left unresolved.
"""
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED_SHA = '04687e2ea261caecee18065dcd2523313aca345ac865eb0483261976f703d85c'
ARMS = ('2L', '2R', '3L', '3R', 'X')
LENS = {'2L': 49364325, '2R': 61545105, '3L': 41963435,
        '3R': 53200684, 'X': 24393108}
EXPECTED_VCF = {'2L': 727, '2R': 395, '3L': 239, '3R': 576, 'X': 69}
BASES = set('ACGT')


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load_inputs(manifest_path, vcf_path):
    digest = sha256(manifest_path)
    if digest != EXPECTED_SHA:
        raise ValueError('STOP: manifest SHA256 differs from frozen version: ' + digest)
    with manifest_path.open(newline='', encoding='utf-8-sig') as f:
        sites = list(csv.DictReader(f))
    if len(sites) != 530 or len({r['population_site_id'] for r in sites}) != 530:
        raise ValueError('STOP: manifest count or site IDs are invalid')
    by_arm = defaultdict(list)
    for r in sites:
        arm = r['genomic_seqid']
        if arm not in (*ARMS, 'UNKN'):
            raise ValueError('STOP: invalid contig: ' + arm)
        start, end = int(r['genomic_start']), int(r['genomic_end'])
        proto, pam = r['protospacer_20nt'].upper(), r['pam'].upper()
        strand = r['guide_genomic_strand']
        if end - start + 1 != 23 or start < 1 or strand not in ('+', '-'):
            raise ValueError('STOP: invalid guide coordinates / strand for ' + r['population_site_id'])
        if len(proto) != 20 or len(pam) != 3 or not set(proto + pam).issubset(BASES):
            raise ValueError('STOP: invalid guide or PAM sequence for ' + r['population_site_id'])
        r['_start'] = start
        r['_end'] = end
        r['_guide'] = proto + pam
        r['_expected_plus'] = r['_guide'] if strand == '+' else revcomp(r['_guide'])
        by_arm[arm].append(r)
    if any(len(by_arm[a]) != n for a, n in {'2L':185,'2R':100,'3L':60,'3R':140,'X':30,'UNKN':15}.items()):
        raise ValueError('STOP: frozen chromosome site counts differ')
    with vcf_path.open(newline='', encoding='utf-8') as f:
        variants = list(csv.DictReader(f, delimiter='\t'))
    if len(variants) != 2006:
        raise ValueError('STOP: expected 2,006 VCF REF entries')
    var_by_arm = defaultdict(list)
    for r in variants:
        arm = r['chromosome']
        if arm not in ARMS or not r['ref'].upper() or set(r['ref'].upper()) - BASES:
            raise ValueError('STOP: unexpected contig or VCF REF allele')
        r['_pos'] = int(r['pos'])
        var_by_arm[arm].append(r)
    if any(len(var_by_arm[a]) != n for a,n in EXPECTED_VCF.items()):
        raise ValueError('STOP: unexpected per-arm VCF count')
    return digest, sites, by_arm, var_by_arm


def fasta_iter(path):
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'rt', encoding='ascii') as f:
        name = None
        chunks = []
        for raw in f:
            if raw.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = raw[1:].split()[0]
                chunks = []
            elif name is None:
                raise ValueError('FASTA not correctly formatted; missing header')
            elif name in ARMS:
                chunks.append(raw.strip().upper())
        if name is not None:
            yield name, ''.join(chunks)


def write_csv(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def run(root=Path('.'), fasta_override=None):
    root = Path(root).resolve()
    manifest = root / '13_population_site_manifest_FROZEN.csv'
    vcf = root / 'output/AR1_PASS_SNP_REF_5arms.tsv'
    fasta = Path(fasta_override).resolve() if fasta_override else root / '../downloads/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz'
    if not manifest.is_file() or not vcf.is_file() or not fasta.is_file():
        raise FileNotFoundError(f'Missing required files: manifest={manifest.is_file()}, VCF REF={vcf.is_file()}, AgamP4 FASTA={fasta.is_file()}')
    digest, sites, by_arm, var_by_arm = load_inputs(manifest, vcf)
    ref_rows, variants_out, mismatches = [], [], []
    observed = set()
    for arm, seq in fasta_iter(fasta):
        if arm not in ARMS:
            continue
        if arm in observed:
            raise ValueError('STOP: duplicate chromosome in FASTA ' + arm)
        observed.add(arm)
        if len(seq) != LENS[arm]:
            raise ValueError(f'STOP: FASTA length of {arm}={len(seq)}, expected {LENS[arm]}')
        for r in by_arm[arm]:
            observed_plus = seq[r['_start']-1:r['_end']]
            expected_plus = r['_expected_plus']
            diff = [i for i in range(23) if observed_plus[i] != expected_plus[i]]
            observed_oriented = observed_plus if r['guide_genomic_strand'] == '+' else revcomp(observed_plus)
            ref_rows.append({'population_site_id':r['population_site_id'], 'gene_id':r['gene_id'],
                             'chromosome':arm, 'start_1based':r['_start'], 'end_1based':r['_end'],
                             'strand':r['guide_genomic_strand'], 'expected_guide_oriented_23bp':r['_guide'],
                             'observed_fasta_guide_oriented_23bp':observed_oriented,
                             'matches_23bp':int(not diff), 'n_matching_bases':23-len(diff),
                             'n_mismatching_bases':len(diff), 'mismatching_genomic_positions':';'.join(str(r['_start']+i) for i in diff),
                             'qc_status': 'ALL_23BP_MATCH' if not diff else 'FASTA_MISMATCH_STOP'})
            for i in diff:
                mismatches.append({'type':'guide_23bp', 'site_id':r['population_site_id'], 'chromosome':arm,
                                   'genomic_pos':r['_start']+i, 'expected':expected_plus[i], 'fasta':observed_plus[i]})
        for v in var_by_arm[arm]:
            pos = v['_pos']
            ref = v['ref'].upper()
            if pos < 1 or pos + len(ref) - 1 > len(seq):
                raise ValueError(f'VCF REF out of bounds: {arm}:{pos}')
            actual = seq[pos-1:pos-1+len(ref)]
            ok = ref == actual
            variants_out.append({'chromosome':arm, 'pos':pos, 'vcf_REF':ref,
                                 'fasta_REF':actual, 'ref_matches':int(ok)})
            if not ok:
                mismatches.append({'type':'vcf_ref', 'site_id':'', 'chromosome':arm,
                                   'genomic_pos':pos, 'expected':ref, 'fasta':actual})
        print(f'{arm}: {len(by_arm[arm])} guide sites + {len(var_by_arm[arm])} VCF REF rows checked; FASTA length {len(seq):,}', flush=True)
        del seq
    if observed != set(ARMS):
        raise ValueError('STOP: FASTA missing required chromosome(s): ' + ','.join(set(ARMS) - observed))
    for r in by_arm['UNKN']:
        ref_rows.append({'population_site_id':r['population_site_id'], 'gene_id':r['gene_id'],
                         'chromosome':'UNKN', 'start_1based':r['_start'], 'end_1based':r['_end'],
                         'strand':r['guide_genomic_strand'], 'expected_guide_oriented_23bp':r['_guide'],
                         'observed_fasta_guide_oriented_23bp':'', 'matches_23bp':'', 'n_matching_bases':'',
                         'n_mismatching_bases':'', 'mismatching_genomic_positions':'',
                         'qc_status':'UNKN_NOT_EVALUABLE'})
    ref_rows.sort(key=lambda r: next(i for i,s in enumerate(sites) if s['population_site_id']==r['population_site_id']))
    out = root/'output'; out.mkdir(exist_ok=True)
    write_csv(out/'AR1_frozen_site_FULL23bp_FASTA_QC_530.csv',
              ['population_site_id','gene_id','chromosome','start_1based','end_1based','strand',
               'expected_guide_oriented_23bp','observed_fasta_guide_oriented_23bp','matches_23bp',
               'n_matching_bases','n_mismatching_bases','mismatching_genomic_positions','qc_status'],ref_rows)
    write_csv(out/'AR1_PASS_2006_VCF_vs_AgamP4_FASTA_REF.csv',
              ['chromosome','pos','vcf_REF','fasta_REF','ref_matches'],variants_out)
    write_csv(out/'AR1_reference_mismatches.csv',
              ['type','site_id','chromosome','genomic_pos','expected','fasta'],mismatches)
    n_full = sum(r['qc_status']=='ALL_23BP_MATCH' for r in ref_rows)
    n_pass = sum(r['ref_matches'] for r in variants_out)
    summary={'manifest_sha256':digest,'reference_name':'Ensembl Metazoa AgamP4 toplevel genome FASTA',
             'reference_fasta_path':str(fasta),'reference_fasta_sha256':sha256(fasta),
             'sites_frozen':530,'sites_evaluable':515,'sites_UNKN_unresolved':15,
             'full_23bp_guide_matches':n_full,'full_23bp_guide_mismatches':515-n_full,
             'guide_bases_compared':11845,'guide_bases_matching':sum(int(r['n_matching_bases']) for r in ref_rows if r['chromosome'] != 'UNKN'),
             'VCF_REF_rows_compared':len(variants_out),'VCF_REF_matches':n_pass,
             'VCF_REF_mismatches':len(variants_out)-n_pass,
             'qc_pass': n_full==515 and n_pass==2006}
    (out/'AR1_full23bp_FASTA_QC_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(f'FULL-23bp reference guide concordance: {n_full}/515; mismatching sites: {515-n_full}')
    print(f'PASS VCF REF vs FASTA concordance: {n_pass}/2006; mismatches: {2006-n_pass}')
    print('Outputs: output/AR1_frozen_site_FULL23bp_FASTA_QC_530.csv; output/AR1_PASS_2006_VCF_vs_AgamP4_FASTA_REF.csv; output/AR1_reference_mismatches.csv; output/AR1_full23bp_FASTA_QC_summary.json')
    print('CAUTION: FASTA concordance does not measure population guide conservation, indels or individual genotype callability.')
    if not summary['qc_pass']:
        raise ValueError('STOP: reference discordances present. Review mismatch table before population variation analysis.')
    return summary


if __name__=='__main__':
    try:
        run(Path(sys.argv[1]) if len(sys.argv)>1 else Path('.'),
            Path(sys.argv[2]) if len(sys.argv)>2 else None)
    except (ValueError, FileNotFoundError, OSError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        sys.exit(2)
