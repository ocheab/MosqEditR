#!/usr/bin/env python3
"""MosqEdit-R R1.3/R2.Q2: fixed-label, fold-matched descriptive PU-label benchmarks.
No model refitting, no experimental Anopheles outcome labels, no independence claim.
Uses data_processed.zip as original analysis provenance. Python executes numeric QC;
R companion independently reproduces the main tables and generates publication figure.
"""
import io, zipfile, csv, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, auc
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent
ARCHIVE=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/data_processed.zip')
OUT=HERE/'tables'; OUT.mkdir(exist_ok=True)
z=zipfile.ZipFile(ARCHIVE)
def read(n,cols=None):return pd.read_csv(io.BytesIO(z.read('data_processed/'+n)),usecols=cols,low_memory=False)

pred=read('08_pu_cv_predictions.csv')
model='M3_expression_orthology_annotation'; algo='xgboost'
sel=pred[(pred.model_set==model)&(pred.algorithm==algo)].copy()
assert len(sel)==41517 and sel.gene_id.nunique()==13839
m3=read('07_pu_M3_expression_orthology_annotation.csv')
features=['expr_ovary','ovary_specificity_delta','female_male_delta','tau_tissue_specificity','dmel_ortholog_count','has_dmel_ortholog']
assert m3.gene_id.nunique()==len(m3)==13839
assert all(c in m3.columns for c in features)
assert not m3.gene_id.duplicated().any()
assert m3.pu_label.isin([0,1]).all()
assert (m3.pu_label.sum()==718)

# All transformations are explicitly outcome-unaware. No FlyBase-derived predictors.
# Missing values are neutral (middle percentile, 0.5); no label-specific imputation.
for c in features:
 col=pd.to_numeric(m3[c],errors='coerce')
 if c=='has_dmel_ortholog':
  # univariate indicator (0 or 1); not a label-independent validation feature
  m3['rank_'+c]=col.fillna(0.5)
 else:
  m3['rank_'+c]=col.rank(method='average',pct=True).fillna(0.5)

m3['simple_expression_mean4']=m3[['rank_expr_ovary','rank_ovary_specificity_delta','rank_female_male_delta','rank_tau_tissue_specificity']].mean(axis=1)
m3['simple_expression_orthology_mean5']=(m3['simple_expression_mean4']*4+m3['rank_dmel_ortholog_count'])/5
m3['simple_expression_orthology_mean5']=m3['simple_expression_orthology_mean5'].astype(float)
scorers={
  'M3_XGBoost_PU_out_of_fold':'pu_score',
  'Single_ovary_expression':'rank_expr_ovary',
  'Single_ovary_specificity':'rank_ovary_specificity_delta',
  'Single_female_male_delta':'rank_female_male_delta',
  'Single_Dmel_ortholog_count':'rank_dmel_ortholog_count',
  'Simple_4_expression_mean':'simple_expression_mean4',
  'Simple_4_expression_plus_orthology_mean5':'simple_expression_orthology_mean5',
  'Ortholog_presence_diagnostic':'rank_has_dmel_ortholog',
}

P=sel.merge(m3[['gene_id','pu_label']+list(set(scorers.values())-{'pu_score'})],on='gene_id',validate='many_to_one',suffixes=('','_source'))
assert (P.pu_label==P.pu_label_source).all()
P.drop(columns='pu_label_source',inplace=True)
assert len(P)==41517 and P[['repeat_id','fold','gene_id']].duplicated().sum()==0

def metrics(y,score):
 y=np.asarray(y,dtype=int);s=np.asarray(score,dtype=float)
 assert np.isfinite(s).all() and set(np.unique(y))=={0,1}
 n=len(y); pos=int(y.sum());prev=pos/n
 idx=np.argsort(-s,kind='stable'); ys=y[idx]
 tp=np.cumsum(ys)
 r=np.r_[0,tp/pos]; precision=np.r_[1,tp/np.arange(1,n+1)]
 auprc=float(np.sum(np.diff(r)*(precision[:-1]+precision[1:])/2))
 # Group tied scores at their unique thresholds. The original line-by-line
 # PR-AUC is spuriously sensitive to the input row order under heavy ties
 # (notably a binary ortholog-present ranking). Use tie-aware metrics in
 # scientific comparisons; preserve original legacy metric for M3 audit.
 prec_t,rec_t,_=precision_recall_curve(y,s)
 tieaware=float(auc(rec_t,prec_t))
 res=dict(auprc_trapezoidal=auprc,auprc_tieaware_trapezoidal=tieaware,average_precision_tieaware=float(average_precision_score(y,s)),auroc=float(roc_auc_score(y,s)),n=n,n_positive=pos,prevalence=prev)
 for k in [10,20,50,100,200]:
  nk=min(k,n); cutoff=s[idx[nk-1]]
  higher=s>cutoff; at_cutoff=s==cutoff
  hits=float(y[higher].sum()+((nk-int(higher.sum()))/int(at_cutoff.sum()))*y[at_cutoff].sum())
  p_at_k=hits/nk
  res[f'positive_at_{k}']=hits
  res[f'precision_at_{k}']=p_at_k
  res[f'recall_at_{k}']=hits/pos
  res[f'enrichment_at_{k}']=p_at_k/prev
 return res

records=[]
from numpy.random import default_rng
rng=default_rng(20261031)
# Permuted out-of-fold score reference distribution, 300 draws/fold.
# Explicit Monte-Carlo random benchmark rather than treating prevalence as
# precisely equal to the ORIGINAL trapezoidal PR-AUC formula.
for (repeat,fold),df in P.groupby(['repeat_id','fold'],sort=True):
 y=df.pu_label.to_numpy()
 for name,column in scorers.items():
  s=df[column].to_numpy(dtype=float)
  records.append(dict(repeat_id=repeat,fold=fold,method=name,**metrics(y,s)))
 model_s=df.pu_score.to_numpy()
 for j in range(300):
  randomized=rng.permutation(model_s)
  records.append(dict(repeat_id=repeat,fold=fold,method='Permuted_scores_null',null_draw=j+1,**metrics(y,randomized)))
fold_df=pd.DataFrame(records)
fold_df.to_csv(OUT/'R13_all_fold_baseline_results_with_random_permutations.csv',index=False)
actual=fold_df[fold_df.method!='Permuted_scores_null'].copy()
null=fold_df[fold_df.method=='Permuted_scores_null'].copy()
num_fields=['auprc_trapezoidal','auprc_tieaware_trapezoidal','average_precision_tieaware','auroc','prevalence']+[f'{metric}_at_{k}' for metric in ['positive','precision','recall','enrichment'] for k in [10,20,50,100,200]]
res=actual.groupby('method',sort=False)[num_fields].agg(['mean','std'])
res.columns=[a+'_'+b for a,b in res.columns]
res=res.reset_index()
res_null=pd.DataFrame([dict(method='Permuted_scores_null',**{field+'_mean':null[field].mean() for field in num_fields},**{field+'_sd':null[field].std() for field in num_fields})])
result=pd.concat([res,res_null],ignore_index=True)
result.to_csv(OUT/'R13_fold_matched_comparison_summary.csv',index=False)

# Same original M3 held-out folds, no alternate fold splitting.
main=actual[actual.method=='M3_XGBoost_PU_out_of_fold'].sort_values(['repeat_id','fold'])
paired=[]
for meth in scorers:
 if meth=='M3_XGBoost_PU_out_of_fold':continue
 row=actual[actual.method==meth].sort_values(['repeat_id','fold'])
 assert ((main[['repeat_id','fold']].to_numpy()==row[['repeat_id','fold']].to_numpy()).all())
 for metric in ['auprc_tieaware_trapezoidal','average_precision_tieaware','auroc','precision_at_50','enrichment_at_50','precision_at_100','enrichment_at_100']:
  diff=main[metric].to_numpy()-row[metric].to_numpy()
  paired.append(dict(comparator=meth,metric=metric,n_fold_assessments=len(diff),mean_m3_minus_comparator=float(diff.mean()),sd_over_15_reused_folds=float(diff.std(ddof=1)),n_differences_positive=int((diff>0).sum()),caution='Repeated folds within same gene pool: not independent; descriptive only'))
pd.DataFrame(paired).to_csv(OUT/'R13_descriptive_paired_differences_15_folds.csv',index=False)

# Top-k OOF gene mean; model OOF averaging across 3 repeat-specific held-out predictions.
pergene=sel.groupby('gene_id',as_index=False).agg(mean_M3_crossfit_PU_score=('pu_score','mean'),n_oof=('pu_score','size'),pu_label=('pu_label','first'))
assert (pergene.n_oof==3).all()
G=pergene.merge(m3[['gene_id','has_dmel_ortholog']+list(set(scorers.values())-{'pu_score'})],on='gene_id',validate='one_to_one')
G['has_dmel_ortholog']=pd.to_numeric(G['has_dmel_ortholog'],errors='coerce')
G.to_csv(OUT/'R13_all_13839_genes_scores_and_comparators.csv',index=False)
for method,scorecol in [('M3_XGBoost_mean_OOF','mean_M3_crossfit_PU_score'),('Single_ovary_expression','rank_expr_ovary'),('Simple_4_expression_mean','simple_expression_mean4'),('Simple_4_expression_plus_orthology_mean5','simple_expression_orthology_mean5')]:
 d=G.sort_values([scorecol,'gene_id'],ascending=[False,True],kind='stable').reset_index(drop=True)
 for k in [20,50,100,200,500,1000]:
  top=d.iloc[:k];hits=int(top.pu_label.sum());prev=float(G.pu_label.mean());
  unlabeled=k-hits
  # "recall on unlabeled" is unidentifiable without experimentally confirmed outcomes.
  records_tmp=dict(method=method,k=k,n_eligible=len(G),n_observed_positive=int(G.pu_label.sum()),expected_random_positive_count=k*prev,observed_positive_count=hits,unlabeled_count=unlabeled,observed_positive_precision=hits/k,observed_positive_recall=hits/G.pu_label.sum(),enrichment_over_prevalence=(hits/k)/prev,frac_without_mapped_Dmel_ortholog=float((top.has_dmel_ortholog!=1).mean()))
  if k==20: tops=[]
  tops.append(records_tmp)
  if k==100 and method=='M3_XGBoost_mean_OOF':
   top[['gene_id','pu_label','mean_M3_crossfit_PU_score','has_dmel_ortholog']].to_csv(OUT/'R13_top100_M3_OOF_PU_label_status.csv',index=False)
 if method=='M3_XGBoost_mean_OOF':alltops=[]
 alltops+=tops
pd.DataFrame(alltops).to_csv(OUT/'R13_full_universe_topk_utility.csv',index=False)

# Original-targeted top-100, a separate panel of 205 pre-selected genes.
short=read('09_top100_prepopulation_shortlist.csv',cols=['gene_id','pu_label','novel_candidate_rank','top100_novel_prepopulation'])
assert len(short)==100 and int(short.pu_label.sum())==100
# Only 2-dimensional test if there are unlabeled targets among 205; source shortlist isn't a genomewide discovery experiment.


# Evaluate the *different* selected CRISPR shortlist separately from global OOF ranks.
selected=read('09_novel_candidate_ranking.csv')
assert len(selected)==199 and int(selected.pu_label.sum())==150
selected=selected.sort_values('novel_candidate_rank')
selected[['gene_id','novel_candidate_rank','final_prepopulation_rank','pu_label','final_integrated_score','pu_consensus_percentile']].to_csv(OUT/'R13_original_199_novel_targeted_label_provenance.csv',index=False)
selected[selected.pu_label==0][['gene_id','novel_candidate_rank','final_prepopulation_rank','final_integrated_score','pu_consensus_percentile']].to_csv(OUT/'R13_49_unlabeled_genes_in_original_targeted_panel.csv',index=False)
with open(OUT/'R13_original_shortlist_summary.json','w') as f:
 json.dump({'original_targeted_nonbenchmark_genes':199,'observed_positive':150,'unlabeled':49,'top_100_observed_positive':100,'first_unlabeled_rank':int(selected[selected.pu_label==0].novel_candidate_rank.min()),'interpretation':'Top-100 enrichment is conditioned on earlier biological/FlyBase evidence selection and is NOT held-out predictive performance'},f,indent=2)

# Verify 15 original fold PR-AUC metrics from original precomputed output
orig=read('08_pu_cv_metrics.csv')
orig=orig[(orig.model_set==model)&(orig.algorithm==algo)]
assert len(orig)==15
recon=main.merge(orig,on=['repeat_id','fold'],validate='one_to_one')
checks={
 'max_abs_auprc_difference':float((recon.auprc_trapezoidal-recon.apparent_auprc).abs().max()),
 'max_abs_auroc_difference':float((recon.auroc-recon.apparent_auroc).abs().max()),
 'max_abs_recall_at50_difference':float((recon.recall_at_50_x-recon.recall_at_50_y).abs().max()) if 'recall_at_50_y' in recon else None
}
print('RECON',checks)
assert checks['max_abs_auprc_difference']<1e-10
assert checks['max_abs_auroc_difference']<1e-10

# Standalone M1 performance (n=9223, DIFFERENT folds and different prevalence)
m1=read('08_pu_metrics_summary.csv')
m1=m1[(m1.model_set=='M1_expression')&(m1.algorithm=='xgboost')]
assert len(m1)==1
m1.to_csv(OUT/'R13_original_M1_expression_only_xgboost_nonpaired.csv',index=False)

# Labels in full original top 100 are all observed positive; OOF top100 can differ.
summary={
 'source_archive':str(ARCHIVE.name),
 'comparison_population':'13839 eligible M3-labelled genes; 718 observed positive, 13121 unlabeled',
 'n_repeat_folds':15,'n_repeats':3,'n_folds_per_repeat':5,
 'n_random_permutations_per_fold':300,
 'random_total_fold_assessments':4500,
 'score_provenance':'original 08_pu_cv_predictions; no model retrained',
 'metric_choice':'Average precision as primary tie-aware comparison; source-compatible trapezoidal AUPRC retained exclusively for exact model-result verification; expected tie-randomized top-k counts for all baseline comparators',
 'comparator_scope':'predefined outcome-unaware rank-percentiles from 07_pu_M3 raw feature values',
 'naive_composite_components':'equal weight expression ovary, ovary specificity, female-male delta, tau; optional equal-weight Dmel ortholog count',
 'all_expression_missing_imputation':'neutral percentile 0.5',
 'label_dependent_predictors_not_used_in_naive_composite':True,
 'original_top100_targeted_all_PU_positive':True,
 'verification':checks,
 'major_limits':['Observed positives derived from FlyBase/Drosophila orthology; unlabeled are NOT true negatives','The non-PU baseline does not refit supervised PU classifiers and is a diagnostic comparator, not an independent efficacy benchmark','All 15 hold-out folds reuse the same gene universe across three repeats and should not be considered independent replicates','Top-k full-universe scores are the mean of 3 OOF predictions, not an independent held-out biological dataset','Original M1 is on 9223 expression-measured genes, not directly paired with M3 13839 pool','No true recall or true positive yield can be computed for untested genes']}
(OUT/'R13_audit_summary.json').write_text(json.dumps(summary,indent=2))
print('\nFOLD-MATCHED BASELINES')
print(result[['method','auprc_tieaware_trapezoidal_mean','average_precision_tieaware_mean','auroc_mean','precision_at_50_mean','enrichment_at_50_mean','positive_at_50_mean']].to_string(index=False,float_format=lambda v:f'{v:.5f}'))
print('\nFULL UNIVERSE TOP K')
print(pd.DataFrame(alltops).to_string(index=False,float_format=lambda v:f'{v:.4f}'))
print('\nFILES',list(OUT.iterdir()))
