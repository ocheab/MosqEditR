#!/usr/bin/env python3
"""MosqEdit-R population-stratified *PASS-SNP* variation, frozen 530-guide panel.

Standard Python 3 only. Does NOT infer haplotypes, whole-23bp sample callability,
indels, observed CRISPR editing efficacy, or invariance from absent SNP records.
Reads the five previously extracted, genotype-bearing VCFs and the existing
Ag1000G Phase 2 sample metadata and frozen guide/accessibility/reference QCs.

Outputs:
  output/AR1_SNP_variant_by_population.csv
  output/AR1_SNP_site_by_population.csv
  output/AR1_SNP_site_summary_530.csv
  output/AR1_SNP_QC_summary.json

Primary: all completely observed GT calls, no newly invented GQ/DP thresholds.
Sensitivity: GT with GQ>=20 and DP>=10 (illustrative, adjustable CLI).
Unphased calls support SNP-wise allele frequencies, not intact haplotype frequencies.
"""
import argparse, csv, gzip, hashlib, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

ARMS = ('2L','2R','3L','3R','X')
EXPECTED_SHA = '04687e2ea261caecee18065dcd2523313aca345ac865eb0483261976f703d85c'
EXPECTED_COUNTS = {'2L':727,'2R':395,'3L':239,'3R':576,'X':69}
SITES_PER_ARM = {'2L':185,'2R':100,'3L':60,'3R':140,'X':30,'UNKN':15}
BASES = set('ACGT')


def stop(msg):
    raise ValueError('STOP: '+str(msg))


def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def read_csv(p):
    if not p.is_file(): stop('missing required input '+str(p))
    with p.open(newline='',encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    path.parent.mkdir(exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def load(root, metadata):
    manifest=root/'13_population_site_manifest_FROZEN.csv'
    digest=sha256(manifest)
    if digest != EXPECTED_SHA: stop('FROZEN MANIFEST SHA mismatch '+digest)
    sites=read_csv(manifest)
    if len(sites)!=530 or len({s['population_site_id'] for s in sites})!=530:
        stop('frozen guide count or duplicate IDs')
    for arm,n in SITES_PER_ARM.items():
        if sum(s['genomic_seqid']==arm for s in sites)!=n: stop('unexpected site count '+arm)
    for s in sites:
        a,b=int(s['genomic_start']),int(s['genomic_end'])
        if a<1 or b-a+1!=23 or s['guide_genomic_strand'] not in ('+','-'):
            stop('bad guide position/strand '+s['population_site_id'])
        s['_s']=a;s['_e']=b;s['_nvar']=0
        s['_cat']=Counter()
    acc=read_csv(root/'output/AR1_frozen_site_accessibility_530.csv')
    ref=read_csv(root/'output/AR1_frozen_site_FULL23bp_FASTA_QC_530.csv')
    accmap={x['population_site_id']:x for x in acc}
    refmap={x['population_site_id']:x for x in ref}
    if len(accmap)!=530 or len(refmap)!=530: stop('expected QC tables of 530 rows')
    for s in sites:
        sid=s['population_site_id']
        if sid not in accmap or sid not in refmap: stop('missing QC site '+sid)
        s['_accessible']=accmap[sid]['fully_accessible_23bp'] if s['genomic_seqid']!='UNKN' else ''
        if s['genomic_seqid']!='UNKN' and refmap[sid]['qc_status']!='ALL_23BP_MATCH':
            stop('guide reference QC not passing '+sid)
    if not metadata.is_file(): stop('sample metadata missing: '+str(metadata))
    with metadata.open(newline='',encoding='utf-8-sig') as f:
        header=f.readline().rstrip('\r\n').split('\t')
        if not {'ox_code','population','country'}.issubset(header):
            stop('metadata missing ox_code/population/country headers; got '+str(header[:8]))
        f.seek(0)
        records=list(csv.DictReader(f,delimiter='\t'))
    byid={r['ox_code']:r for r in records}
    if len(byid)!=len(records): stop('duplicate sample IDs in metadata')
    return sites,byid,digest


def genomic_region(site,pos):
    # Guide orientation is [20-nt protospacer][3-nt NGG PAM].
    # A genomic minus-strand guide has PAM at leftmost (low-coordinate) 3 bp.
    gpos=pos-site['_s']+1 if site['guide_genomic_strand']=='+' else site['_e']-pos+1
    if not 1<=gpos<=23: stop('SNP position outside guide '+site['population_site_id'])
    if gpos>=22: return 'PAM_ESSENTIAL_GG'
    if gpos==21: return 'PAM_N_NONESSENTIAL'
    if gpos>=11: return 'PROTOSPACER_PAM_PROXIMAL_10'
    return 'PROTOSPACER_DISTAL_10'


def parse_gt(field, gt_idx, gq_idx, dp_idx, n_alts, min_gq, min_dp):
    tokens=field.split(':')
    gt=tokens[gt_idx] if gt_idx < len(tokens) else '.'
    phase='|' if '|' in gt else '/' if '/' in gt else ''
    parts=gt.split(phase) if phase else [gt]
    if len(parts)>2: stop('Unexpected ploidy > 2 in GT='+gt)
    if not parts or any(not x.isdigit() for x in parts):
        return 0,0,False,False,len(parts),phase,False,False
    alleles=[int(x) for x in parts]
    if any(x>n_alts for x in alleles): stop('GT allele index exceeds ALT count: '+gt)
    gq=tokens[gq_idx] if gq_idx is not None and gq_idx<len(tokens) else '.'
    dp=tokens[dp_idx] if dp_idx is not None and dp_idx<len(tokens) else '.'
    has_gq=gq.isdecimal();has_dp=dp.isdecimal()
    quality_ok=(has_gq and has_dp and int(gq)>=min_gq and int(dp)>=min_dp)
    return sum(a>0 for a in alleles),len(alleles),True,quality_ok,len(alleles),phase,has_gq,has_dp


def main(root,metadata,min_gq,min_dp):
    sites,byid,digest=load(root,metadata)
    byarm={arm:[s for s in sites if s['genomic_seqid']==arm] for arm in ARMS}
    variant_rows=[]; sample_ids=None; pop_idx=None; pop_country=None
    nvariants=Counter(); genotype_qc=Counter(); sample_n=0
    # For each site store bytearray of per-sample flags over *observed PASS SNP positions*.
    # Bit 1 any nonREF; 2 missing GT; 4 any essential PAM GG nonREF; 8 PAM-proximal 10-nt nonREF.
    # Bit 16 at least one GT below illustrative QC threshold (including missing GQ/DP).
    state={}
    for arm in ARMS:
        path=root/f'output/AR1_frozen_{arm}_genotypes.vcf.gz'
        if not path.is_file():stop('missing genotype VCF '+str(path))
        print('Reading',arm,flush=True)
        count=0
        with gzip.open(path,'rt',encoding='ascii') as f:
            for line in f:
                if line.startswith('##'):continue
                if line.startswith('#CHROM'):
                    columns=line.rstrip('\r\n').split('\t')
                    current_ids=columns[9:]
                    if len(current_ids)!=1142 or len(set(current_ids))!=1142:
                        stop('unexpected VCF sample columns '+arm)
                    if sample_ids is None:
                        sample_ids=current_ids;sample_n=len(sample_ids)
                        if any(s not in byid for s in sample_ids):stop('unmatched sample ID')
                        populations=sorted({byid[s]['population'] for s in sample_ids})
                        if any(not p for p in populations): stop('missing population label')
                        pop_idx={p:[i for i,s in enumerate(sample_ids) if byid[s]['population']==p] for p in populations}
                        pop_country={p:','.join(sorted({byid[s]['country'] for s in sample_ids if byid[s]['population']==p})) for p in populations}
                        for s in sites:
                            if s['genomic_seqid'] in ARMS:
                                state[s['population_site_id']]=bytearray(sample_n)
                    elif sample_ids!=current_ids:
                        stop('VCF sample identifiers or order mismatch '+arm)
                    continue
                if not line.strip():continue
                c=line.rstrip('\r\n').split('\t');count+=1
                if len(c)!=9+sample_n:stop('wrong VCF column count '+arm)
                chrom,pos,ref,alt=c[0],int(c[1]),c[3].upper(),c[4].upper()
                if chrom!=arm:stop('unexpected chromosome '+chrom)
                alts=alt.split(',')
                if len(ref)!=1 or ref not in BASES or any(len(x)!=1 or x not in BASES for x in alts):
                    stop('Not a SNP record: '+f'{arm}:{pos} {ref}>{alt}; investigate separately')
                formatkeys=c[8].split(':')
                if 'GT' not in formatkeys:stop('missing GT format '+arm)
                gt_i=formatkeys.index('GT');gq_i=formatkeys.index('GQ') if 'GQ' in formatkeys else None
                dp_i=formatkeys.index('DP') if 'DP' in formatkeys else None
                # Decode samples once, sharing values across overlapping target sites.
                values=[parse_gt(x,gt_i,gq_i,dp_i,len(alts),min_gq,min_dp) for x in c[9:]]
                # Per-pop SNP frequency: counts every callable allele; excludes partial GTs.
                for p,ids in pop_idx.items():
                    vv=[values[i] for i in ids]
                    an=sum(v[1] for v in vv if v[2]);ac=sum(v[0] for v in vv if v[2]);ncalled=sum(v[2] for v in vv)
                    n_nonref=sum(v[0]>0 for v in vv if v[2])
                    anq=sum(v[1] for v in vv if v[2] and v[3]);acq=sum(v[0] for v in vv if v[2] and v[3]);nq=sum(v[2] and v[3] for v in vv)
                    variant_rows.append({'chromosome':arm,'pos_1based':pos,'REF':ref,'ALT':alt,'population':p,
                       'country':pop_country[p],'n_samples':len(ids),'n_GT_called':ncalled,
                       'n_GT_missing_or_partial':len(ids)-ncalled,'n_alt_carriers':n_nonref,
                       'n_alleles_called':an,'n_nonref_alleles':ac,
                       'nonref_AF':round(ac/an,8) if an else '',
                       'n_GT_called_strict_sensitivity':nq,'n_alleles_strict_sensitivity':anq,
                       'n_nonref_alleles_strict_sensitivity':acq,
                       'nonref_AF_strict_sensitivity':round(acq/anq,8) if anq else '',
                       'n_observed_haploid_GT':sum(v[2] and v[4]==1 for v in vv),
                       'n_observed_phased_GT':sum(v[2] and v[5]=='|' for v in vv)})
                for v in values:
                    if v[2]:
                        genotype_qc['called']+=1
                        genotype_qc['haploid' if v[4]==1 else 'diploid']+=1
                        genotype_qc['phased' if v[5]=='|' else 'unphased_or_haploid']+=1
                        if not v[6]: genotype_qc['called_missing_GQ']+=1
                        if not v[7]: genotype_qc['called_missing_DP']+=1
                        if v[3]:genotype_qc['strict_pass']+=1
                    else: genotype_qc['missing_or_partial']+=1
                hits=[s for s in byarm[arm] if s['_s']<=pos<=s['_e']]
                if not hits:stop('VCF locus outside all frozen sites '+f'{arm}:{pos}')
                for s in hits:
                    s['_nvar']+=1
                    region=genomic_region(s,pos)
                    s['_cat'][region]+=1
                    flags=state[s['population_site_id']]
                    for i,v in enumerate(values):
                        alt_n,an,called,strict,ploidy,phased,hasgq,hasdp=v
                        if not called:flags[i]|=2
                        elif alt_n:
                            flags[i]|=1
                            if region=='PAM_ESSENTIAL_GG':flags[i]|=4
                            if region=='PROTOSPACER_PAM_PROXIMAL_10':flags[i]|=8
                        if not strict:flags[i]|=16
        if count!=EXPECTED_COUNTS[arm]:stop(f'{arm} expected {EXPECTED_COUNTS[arm]} variants, got {count}')
        nvariants[arm]=count
        print(f'  {count} PASS SNP records, {sample_n} matched samples',flush=True)
    site_rows=[]; site_pop_rows=[]
    for s in sites:
        sid=s['population_site_id'];arm=s['genomic_seqid'];nvar=s['_nvar']
        allacc=s['_accessible']
        site_rows.append({'population_site_id':sid,'gene_id':s['gene_id'], 'chromosome':arm,
           'start_1based':s['_s'],'end_1based':s['_e'],'strand':s['guide_genomic_strand'],
           'fully_accessible_23bp':allacc, 'n_observed_PASS_SNP_records':nvar if arm!='UNKN' else '',
           'n_essential_PAM_GG_SNP_positions':s['_cat']['PAM_ESSENTIAL_GG'] if arm!='UNKN' else '',
           'n_PAM_N_SNP_positions':s['_cat']['PAM_N_NONESSENTIAL'] if arm!='UNKN' else '',
           'n_protospacer_PAM_proximal10_SNP_positions':s['_cat']['PROTOSPACER_PAM_PROXIMAL_10'] if arm!='UNKN' else '',
           'n_protospacer_distal10_SNP_positions':s['_cat']['PROTOSPACER_DISTAL_10'] if arm!='UNKN' else '',
           'interpretation_status':'UNKN_NOT_EVALUABLE' if arm=='UNKN' else
             'NO_REPORTED_PASS_SNP_NOT_PROOF_OF_INVARIANCE' if nvar==0 else
             'INCOMPLETE_SITE_ACCESSIBILITY' if allacc!='1' else 'PASS_SNP_LOCI_OBSERVED_NOT_FULL_SITE_CALLABILITY'})
        for p,ids in pop_idx.items():
            flags=state.get(sid)
            result={'population_site_id':sid,'gene_id':s['gene_id'],'chromosome':arm,'population':p,
                'country':pop_country[p],'n_samples':len(ids),'fully_accessible_23bp':allacc,
                'n_observed_PASS_SNP_records':nvar if arm!='UNKN' else ''}
            if arm=='UNKN' or not nvar:
                result.update({f:'' for f in ('n_complete_GT_at_reported_SNP_loci','n_missing_GT_at_reported_SNP_loci',
                    'n_any_nonref_carriers','n_no_nonref_at_reported_SNP_loci',
                    'fraction_any_nonref_among_complete_GT','n_essential_PAM_GG_nonref_carriers',
                    'n_PAM_proximal10_nonref_carriers','n_complete_GT_strict_sensitivity',
                    'n_any_nonref_strict_sensitivity','fraction_any_nonref_strict_sensitivity')})
                result['interpretation_status']='UNKN_NOT_EVALUABLE' if arm=='UNKN' else 'NO_PASS_SNP_RECORDS_NO_GT_EVIDENCE'
            else:
                complete=[i for i in ids if not(flags[i]&2)]
                complete_strict=[i for i in ids if not(flags[i]&(2|16))]
                nalt=sum(bool(flags[i]&1) for i in complete)
                naltq=sum(bool(flags[i]&1) for i in complete_strict)
                result.update({'n_complete_GT_at_reported_SNP_loci':len(complete),
                    'n_missing_GT_at_reported_SNP_loci':len(ids)-len(complete),
                    'n_any_nonref_carriers':nalt,
                    'n_no_nonref_at_reported_SNP_loci':len(complete)-nalt,
                    'fraction_any_nonref_among_complete_GT':round(nalt/len(complete),8) if complete else '',
                    'n_essential_PAM_GG_nonref_carriers':sum(bool(flags[i]&4) for i in complete),
                    'n_PAM_proximal10_nonref_carriers':sum(bool(flags[i]&8) for i in complete),
                    'n_complete_GT_strict_sensitivity':len(complete_strict),
                    'n_any_nonref_strict_sensitivity':naltq,
                    'fraction_any_nonref_strict_sensitivity':round(naltq/len(complete_strict),8) if complete_strict else '',
                    'interpretation_status':'PASS_SNP_LOCI_ONLY' if allacc=='1' else 'PARTIALLY_INACCESSIBLE_PASS_SNP_LOCI_ONLY'})
            site_pop_rows.append(result)
    n_snpsites=sum(s['_nvar']>0 for s in sites if s['genomic_seqid']!='UNKN')
    if n_snpsites != 479:
        stop(f'expected 479 guides overlapping observed PASS SNPs from prior reference QC, found {n_snpsites}')
    out=root/'output'
    fields_site=list(site_rows[0]);fields_pop=list(site_pop_rows[0]);fields_var=list(variant_rows[0])
    write_csv(out/'AR1_SNP_site_summary_530.csv',site_rows,fields_site)
    write_csv(out/'AR1_SNP_site_by_population.csv',site_pop_rows,fields_pop)
    write_csv(out/'AR1_SNP_variant_by_population.csv',variant_rows,fields_var)
    summary={'source':'Ag1000G Phase 2 AR1, PASS SNP genotype-bearing VCFs',
      'frozen_manifest_SHA256':digest,'n_frozen_sites':530,'n_arm_sites':515,
      'n_UNKN_unresolved':15,'n_fully_accessible_23bp':sum(s['_accessible']=='1' for s in sites),
      'n_sites_with_observed_PASS_SNP':sum(s['_nvar']>0 for s in sites),
      'n_sites_without_observed_PASS_SNP':sum(s['_nvar']==0 for s in sites if s['genomic_seqid']!='UNKN'),
      'n_SNP_records_by_arm':dict(nvariants),'n_unique_mosquito_samples':sample_n,
      'n_populations':len(pop_idx),'sample_counts_by_population':{p:len(ids) for p,ids in pop_idx.items()},
      'genotype_QC_over_SNP_sample_cells':dict(genotype_qc),
      'strict_sensitivity_thresholds':{'GQ_at_least':min_gq,'DP_at_least':min_dp},
      'interpretation_limits':['SNPs only, no indel validation','Absence of a PASS SNP does not prove reference genotype at every base',
         'Site-level genotypes cover observed PASS SNP loci only, not per-sample 23-base callability',
         'Unphased GT does not establish complete guide haplotypes','X-ploidy and sex annotation require separate review',
         'No inference of experimental CRISPR editing efficacy','15 UNKN sites unresolved']}
    (out/'AR1_SNP_QC_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(f"PASS SNPs: {sum(nvariants.values())}; frozen guides with >=1 observed PASS SNP: {summary['n_sites_with_observed_PASS_SNP']}/515")
    print(f"Populations: {len(pop_idx)}; sample IDs linked: {sample_n}/{sample_n}; fully accessible: {summary['n_fully_accessible_23bp']}/515")
    print('Genotype-call QC:',dict(genotype_qc))
    print('Outputs in output/: AR1_SNP_site_summary_530.csv, AR1_SNP_site_by_population.csv, AR1_SNP_variant_by_population.csv, AR1_SNP_QC_summary.json')
    print('CAUTION: PASS SNP variation at reported loci != complete target conservation; GTs may be unphased.')


if __name__=='__main__':
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    arg.add_argument('--metadata',type=Path,default=None)
    arg.add_argument('--min-gq',type=int,default=20,help='Sensitivity-only illustrative threshold (primary uses all complete GT)')
    arg.add_argument('--min-dp',type=int,default=10,help='Sensitivity-only illustrative threshold (primary uses all complete GT)')
    a=arg.parse_args()
    root=a.root.resolve()
    metadata=a.metadata if a.metadata else next((p for p in (root/'../downloads/samples.meta.txt',root/'downloads/samples.meta.txt') if p.is_file()),root/'../downloads/samples.meta.txt')
    try: main(root,metadata,a.min_gq,a.min_dp)
    except (ValueError,KeyError,OSError) as e:
        print(e,file=sys.stderr);sys.exit(2)
