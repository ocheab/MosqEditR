#!/usr/bin/env python3
"""Executed deterministic R2.Q3 sensitivity, frozen Step 04 outputs.
Do not call these analyses independent RNA-seq validation or PU-model refits.
"""
from pathlib import Path
import zipfile,io,json,hashlib,argparse
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)
parser=argparse.ArgumentParser()
parser.add_argument('--processed-zip',default=None,help='Optional original data_processed.zip; otherwise uses packaged inputs/ CSVs')
args=parser.parse_args()
source_path=Path(args.processed_zip).resolve() if args.processed_zip else None
Z=zipfile.ZipFile(source_path) if source_path else None
def load(n):
    if Z is not None:
        return pd.read_csv(io.BytesIO(Z.read('data_processed/'+n+'.csv')),low_memory=False)
    file=ROOT/'inputs'/f'{n}.csv'
    if not file.exists():raise FileNotFoundError(f'Missing input {file}; use --processed-zip /path/to/data_processed.zip')
    return pd.read_csv(file,low_memory=False)

d=load('04_preliminary_genomewide_ranking')
t=load('09_novel_candidate_ranking')
assert len(d)==13845 and d.gene_id.is_unique
assert all(x in d.columns for x in ['fris','functional_evidence_score','score_primary','evidence_domain_coverage','rank_primary','preliminary_rank'])
E=pd.to_numeric(d['fris'],errors='coerce').to_numpy(dtype=float)
F=pd.to_numeric(d['functional_evidence_score'],errors='coerce').to_numpy(dtype=float)
genes=d['gene_id'].astype(str).to_numpy()
assert np.isfinite(E).sum()==9228 and np.isfinite(F).sum()==8537
np.savez_compressed(OUT/'sensitivity_base_values.npz',expr=E,function=F) # internal numeric reproducibility only

def rank_scores(vals):
    # Original Step04 uses decreasing score, then evidence completeness, gene_id; keep tie handling exact.
    ok=np.isfinite(vals)
    completeness=pd.to_numeric(d['evidence_feature_completeness'],errors='coerce').to_numpy()
    order=np.lexsort((genes,-np.nan_to_num(completeness,nan=-np.inf),-np.where(ok,vals,-np.inf)))
    ranks=np.full(len(vals),np.nan)
    ranks[order]=np.arange(1,len(vals)+1)
    ranks[~ok]=np.nan
    return ranks

def score(w,E=E,F=F):
    mE=np.isfinite(E);mF=np.isfinite(F)
    used=w*mE + (1-w)*mF
    raw=(w*np.where(mE,E,0)+(1-w)*np.where(mF,F,0))/np.where(used>0,used,1)
    return np.where(used>0,raw*(.75+.25*used),np.nan)

ref=score(.60)
mask=np.isfinite(ref)&np.isfinite(d['score_primary'])
assert mask.sum()==10925
absdiff=np.max(np.abs(ref[mask]-d.loc[mask,'score_primary'].to_numpy()))
assert absdiff<1e-9,('original Step04 score not reproduced',absdiff)
ref_rank=rank_scores(ref)
original_rank=d['rank_primary'].to_numpy(dtype=float)
ref_corr=spearmanr(ref_rank[np.isfinite(original_rank)], original_rank[np.isfinite(original_rank)]).statistic
assert ref_corr>0.999,ref_corr

def top_ids(rank,k,valid=None):
    pos=np.isfinite(rank)
    if valid is not None:pos=pos & valid
    idx=np.where(pos)[0]
    return set(genes[idx[np.argsort(rank[idx],kind='stable')[:k]]])

base20=top_ids(original_rank,20);base50=top_ids(original_rank,50);base100=top_ids(original_rank,100);base200=top_ids(original_rank,200)
score_rows=[];gene_rows=[]
for w in [0.0,.15,.30,.45,.60,.75,.90,1.0]:
    s=score(w);r=rank_scores(s)
    common=np.isfinite(original_rank)&np.isfinite(r)
    row={'expression_weight':w, 'functional_weight':round(1-w,2),
         'n_ranked_genes':int(np.isfinite(r).sum()),'n_common_ranked':int(common.sum()),
         'spearman_common_with_original':round(float(spearmanr(original_rank[common],r[common]).statistic),6)}
    for k,base in [(20,base20),(50,base50),(100,base100),(200,base200)]:
        now=top_ids(r,k)
        row[f'top{k}_overlap_n']=len(base&now)
        row[f'top{k}_overlap_frac']=len(base&now)/k
    score_rows.append(row)
    if w in [0,.30,.60,.90]:
        for i in range(len(genes)):
            gene_rows.append({'gene_id':genes[i], 'expression_weight':w,'score':s[i], 'rank':r[i], 'original_rank':original_rank[i], 'expression_observed':bool(np.isfinite(E[i])), 'functional_observed':bool(np.isfinite(F[i]))})

sr=pd.DataFrame(score_rows)
sr.to_csv(OUT/'R2Q3_step04_domain_weight_sensitivity.csv',index=False)
pd.DataFrame(gene_rows).to_csv(OUT/'R2Q3_step04_gene_rank_scenarios.csv',index=False)

# Candidate-specific comparisons within fixed original targeted panel, NOT rerun Step05 PU or step09.
gene_idx={g:i for i,g in enumerate(genes)}
t=t.sort_values('novel_candidate_rank',na_position='last')
original_target=set(t.loc[t['is_external_benchmark']==False,'gene_id'].astype(str))
assert len(original_target)==199
tr=[]
for _,row in t.iterrows():
    g=row.gene_id;i=gene_idx[g]
    out={'gene_id':g,'is_external_benchmark':bool(row['is_external_benchmark']), 'original_integrated_rank':row['final_integrated_rank'], 'novel_candidate_rank':row['novel_candidate_rank'], 'pu_label':row['pu_label'], 'moz_expression_available':np.isfinite(E[i]),'n_probes_mapped':None,'fris':E[i], 'original_preliminary_rank':d.iloc[i]['preliminary_rank']}
    for w in [.0,.30,.60,.90]:out[f'step04_rank_expr_weight_{int(w*100)}']=rank_scores(score(w))[i]
    tr.append(out)
tr=pd.DataFrame(tr);tr.to_csv(OUT/'R2Q3_fixed_target_panel_preliminary_rank_sensitivity.csv',index=False)

# Missingness strata
p=load('02_mozatlas_gene_expression')[['gene_id','n_probes_mapped']]
d2=d[['gene_id','is_prespecified_benchmark','preliminary_rank','fris','functional_evidence_score']].merge(p,on='gene_id',how='left',validate='one_to_one')
groups={'Whole annotated gene universe':d2,'Ranked preliminary top 100':d2.loc[d2.preliminary_rank.le(100)],'Original top 100 non-benchmark integrated':d2.loc[d2.gene_id.isin(t[t['novel_candidate_rank'].le(100)]['gene_id'])], 'Original targeted 199 non-benchmarks':d2.loc[d2.gene_id.isin(original_target)]}
cover=[]
for group,sub in groups.items():
    cover.append({'cohort':group,'n_genes':len(sub), 'n_mozatlas_expression':int(sub.fris.notna().sum()),'pct_mozatlas_expression':round(100*sub.fris.notna().mean(),3),'n_flybase_functional':int(sub.functional_evidence_score.notna().sum()),'n_without_expression':int(sub.fris.isna().sum()),'n_single_probe':int(sub.n_probes_mapped.eq(1).sum()), 'n_multiple_probes':int(sub.n_probes_mapped.gt(1).sum()),'n_probes_unknown':int(sub.n_probes_mapped.isna().sum())})
pd.DataFrame(cover).to_csv(OUT/'R2Q3_expression_coverage_by_cohort.csv',index=False)

# Dropout is an *illustrative availability stress test*; not re-fit or RNAseq substitution.
rng=np.random.default_rng(20261031)
obs_idx=np.where(np.isfinite(E))[0]
mask_rows=[];replicates=150
for loss in [0.10,0.25,0.50]:
    for rep in range(replicates):
        idx=rng.choice(obs_idx,int(round(loss*len(obs_idx))),replace=False)
        e=E.copy();e[idx]=np.nan
        ss=score(.6,E=e); rr=rank_scores(ss)
        common=np.isfinite(rr)&np.isfinite(original_rank)
        row={'fraction_of_available_expression_masked':loss,'replicate':rep+1, 'n_expression_removed':len(idx),'n_remaining_ranked':int(np.isfinite(rr).sum()),'rank_spearman_common':float(spearmanr(original_rank[common],rr[common]).statistic)}
        for k,base in [(20,base20),(50,base50),(100,base100),(200,base200)]:row[f'top{k}_overlap_n']=len(top_ids(rr,k)&base)
        mask_rows.append(row)
mask_data=pd.DataFrame(mask_rows)
mask_data.to_csv(OUT/'R2Q3_expression_masking_replicates.csv',index=False)
summary=[]
for fraction,sub in mask_data.groupby('fraction_of_available_expression_masked'):
    dct={'fraction_of_available_expression_masked':fraction,'n_replicates':len(sub)}
    for key in ['rank_spearman_common','top20_overlap_n','top50_overlap_n','top100_overlap_n','top200_overlap_n']:
        dct[key+'_median']=float(sub[key].median());dct[key+'_p05']=float(sub[key].quantile(.05));dct[key+'_p95']=float(sub[key].quantile(.95))
    summary.append(dct)
pd.DataFrame(summary).to_csv(OUT/'R2Q3_expression_masking_summary.csv',index=False)

# source evidence selection, not analyzed yet
sources=[{'accession':'GSE21689','year':2011,'assay':'Affymetrix microarray','tissue': '15 adult tissue/sex groups including ovaries/testis','scope':'Original expression features','status':'analyzed_original'},{'accession':'GSE239600','year':2023,'assay':'bulk RNA-seq','tissue':'WT ovary n=4, WT fat body n=2, WT whole body n=3; blood fed','scope':'External WT ovary expression concordance only; not whole-tissue replacement','status':'verified_public_processed_counts_not_downloaded'},{'accession':'GSE243215','year':2024,'assay':'bulk RNA-seq (Salmon transcript quantification)','tissue':'salivary glands infected and non-infected','scope':'Insufficient ovarian coverage for primary R2Q3 test','status':'verified_not_used'},{'accession':'GSE284722','year':2025,'assay':'bulk RNA-seq','tissue':'female aging and temperature conditions','scope':'Biologically nonmatched context; possible extension','status':'verified_not_used'}]
pd.DataFrame(sources).to_csv(OUT/'R2Q3_external_expression_source_inventory.csv',index=False)

meta={'script':'01_audit_sensitivity.py','original_source':str(source_path) if source_path else 'inputs/ original extracted CSV files','source_zip_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest() if source_path else None,'input_csv_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'inputs').glob('*.csv')},'n_universe':len(d),'n_expression_mapped':int(np.isfinite(E).sum()),'n_functional_evidence':int(np.isfinite(F).sum()),'n_both_domains':int((np.isfinite(E)&np.isfinite(F)).sum()),'n_no_domain':int((~np.isfinite(E)&~np.isfinite(F)).sum()),'original_score_max_abs_diff':float(absdiff),'original_rank_spearman_check':float(ref_corr),'n_targeted_nonbenchmark':len(original_target),'expression_weight_primary':0.60,'functional_weight_primary':0.40,'completeness_floor':0.75,'masking_reps_per_fraction':replicates,'external_RNAseq_results': 'PENDING: download publicly available GSE239600 processed gene-expression table'}
(OUT/'R2Q3_audit_summary.json').write_text(json.dumps(meta,indent=2))

print('SOURCE AUDIT',json.dumps(meta,indent=2))
print('\nCOVERAGE\n',pd.DataFrame(cover).to_string(index=False))
print('\nWEIGHT\n',sr.to_string(index=False))
print('\nEXPRESSION DROPOUT\n',pd.DataFrame(summary).to_string(index=False))
