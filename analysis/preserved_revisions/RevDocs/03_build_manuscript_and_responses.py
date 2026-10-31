from pathlib import Path
import pandas as pd
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parent;T=ROOT/'tables'
f=pd.read_csv(T/'R13_fold_matched_comparison_summary.csv').set_index('method')
top=pd.read_csv(T/'R13_full_universe_topk_utility.csv')

def fmt(x,n=3):return f'{x:.{n}f}'

def setup(d,kind):
 sec=d.sections[0];sec.top_margin=Inches(.8);sec.bottom_margin=Inches(.8);sec.left_margin=Inches(.92);sec.right_margin=Inches(.82)
 sty=d.styles['Normal'];sty.font.name='Times New Roman';sty.font.size=Pt(12)
 sty.paragraph_format.line_spacing=2.0;sty.paragraph_format.space_after=Pt(7)
 for h in ['Title','Heading 1','Heading 2']:
  s=d.styles[h];s.font.name='Times New Roman';s.font.color.rgb=RGBColor(27,49,69)
 d.core_properties.author='MosqEdit-R study team'
 d.core_properties.title=kind
 d.add_heading('MosqEdit-R | R1.3 and R2.Q2',0)
 p=d.add_paragraph('Executed performance baselines, top-k proxy utility and interpretation of the original non-benchmark shortlist')
 p.style='Subtitle'
 p=d.add_paragraph('Evidence audit prepared 31 October 2026  |  Manuscript PCOMPBIOL-D-26-02293')
 p.style='Caption'

def para(d,text):d.add_paragraph(text)

def table(d,headers,rows,widths=None):
 t=d.add_table(rows=1, cols=len(headers));t.style='Light Shading Accent 1';t.alignment=WD_TABLE_ALIGNMENT.CENTER
 for j,h in enumerate(headers):t.rows[0].cells[j].text=str(h)
 for row in rows:
  cs=t.add_row().cells
  for j,val in enumerate(row):cs[j].text=str(val)
 for ri,row in enumerate(t.rows):
  for ci,c in enumerate(row.cells):
   c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
   for p in c.paragraphs:
    p.paragraph_format.line_spacing=1.04;p.paragraph_format.space_after=Pt(3)
    for run in p.runs:run.font.name='Times New Roman';run.font.size=Pt(9 if len(headers)>=5 else 10)
  if ri==0:
   for c in row.cells:
    for p in c.paragraphs:
     for run in p.runs:run.bold=True
  if ri==0:
   trPr=row._tr.get_or_add_trPr();hdr=OxmlElement('w:tblHeader');hdr.set(qn('w:val'),'true');trPr.append(hdr)
 d.add_paragraph('')
 return t

DISPLAY={
'M3_XGBoost_PU_out_of_fold':'M3-XGBoost (PU; OOF)',
'Permuted_scores_null':'Random permutation',
'Single_ovary_expression':'Ovary expression',
'Single_ovary_specificity':'Ovary specificity',
'Single_female_male_delta':'Female–male expression',
'Simple_4_expression_mean':'Four-expression composite',
'Simple_4_expression_plus_orthology_mean5':'Expression + orthology composite',
'Single_Dmel_ortholog_count':'Drosophila ortholog count',
'Ortholog_presence_diagnostic':'Ortholog-present indicator (diagnostic)'
}
order=['M3_XGBoost_PU_out_of_fold','Permuted_scores_null','Single_ovary_expression','Single_ovary_specificity','Single_female_male_delta','Simple_4_expression_mean','Simple_4_expression_plus_orthology_mean5','Single_Dmel_ortholog_count']
rows=[]
for name in order:
 r=f.loc[name]
 rows.append((DISPLAY[name],fmt(r.average_precision_tieaware_mean),fmt(r.auroc_mean),fmt(r.positive_at_50_mean,2),fmt(r.enrichment_at_50_mean,2)+'×'))

D=Document();setup(D,'Manuscript insertion — performance baselines and top-k utility')
D.add_heading('Methods: fold-matched performance baselines and label-recovery utility',1)
para(D,'To contextualize PU proxy discrimination, the original M3 expression–orthology–annotation XGBoost out-of-fold predictions were reanalysed without refitting or changing model selection. The 13,839 eligible non-landmark genes comprised 718 FlyBase-derived observed positives and 13,121 unlabeled genes; unlabeled status was not treated as a confirmed negative phenotype. Three repetitions of five-fold cross-fitting yielded 15 held-out fold assessments (approximately 2,768 genes per fold). All diagnostic comparators were evaluated in exactly these held-out gene subsets.')
para(D,'We compared M3 with a permutation null (300 independent within-fold random reorderings of the original prediction scores), single-feature rankings by ovary expression, ovary specificity and female–male expression contrast, and two deliberately simple non-PU heuristics. The four-feature expression heuristic was the arithmetic mean of outcome-unaware rank percentiles for ovary abundance, ovary specificity, female–male contrast and tissue specificity (tau). The five-feature heuristic added the mapped Drosophila ortholog-count percentile at equal weight. Missing feature values received neutral percentile 0.5. These exploratory comparators were specified for the present revision analysis, not selected through outcome-based weight fitting. Variables derived from the FlyBase female-sterility outcome were expressly excluded from the simple baselines. The five-feature heuristic nevertheless retains comparative-orthology ascertainment dependence.')
para(D,'The primary cross-method metric was tie-aware average precision (AP), defined by summing the precision at each distinct score threshold weighted by its corresponding recall increment. The manuscript’s original trapezoidal precision–recall AUC was independently recalculated to verify backward compatibility; because these two integration conventions differ, they are reported with different names. AUROC, expected observed-positive yield among the highest-ranked k genes (k=10, 20, 50, 100 and 200), observed-positive precision-at-k, and label-enrichment relative to each fold’s prevalence were also estimated. Where a heuristic produced scores tied at the kth boundary, expected top-k label counts average uniformly over boundary-tied genes, preventing row-order artefacts. Means and descriptive standard deviations were calculated over 15 repeatedly sampled folds; these assessments reuse genes and are not independent experimental replicates.')
para(D,'A complementary genome-wide descriptive analysis averaged each gene’s three held-out M3 predictions to obtain one cross-fitted gene score and tabulated observed-positive versus unlabeled status at top-k thresholds. This ranked 13,839 modelling-eligible genes, not the distinct 205-gene CRISPR-screened panel. The previously frozen non-benchmark integrated shortlist (199 genes; 100 reported in its top tier) was evaluated separately. Neither comparison can measure recall among truly functional yet unlabeled Anopheles genes, because independent outcome labels for those genes are not available.')

D.add_heading('Results: discriminatory performance relative to simple baselines',1)
para(D,'The original mean apparent precision–recall AUC for M3-XGBoost was reproduced exactly (0.158535; maximum absolute fold-level difference <5×10⁻¹⁶). On a tie-aware AP scale used uniformly for all comparators, M3 achieved mean AP 0.1618 and AUROC 0.7976, compared with AP 0.0545 and AUROC 0.5006 under score permutations. The observed-positive prevalence was approximately 0.0519 per fold. M3 also exceeded ovary abundance (AP 0.0852), ovary specificity (AP 0.0960), the expression-only four-feature mean (AP 0.0984) and the expression-plus-orthology heuristic (AP 0.1015). These diagnostic results indicate improved ranking of the orthology-derived observed-positive labels relative to the tested untuned heuristics, not independently validated prediction of mosquito fertility effects.')
D.add_paragraph('Table X. Fold-matched diagnostic PU-label recovery (means across 15 reused held-out folds).',style='Caption')
table(D,['Ranking method','AP (tie-aware)','AUROC','Positive labels in top 50¹','Top-50 label enrichment'],rows)
para(D,'¹ Top-50 values are mean observed-positive counts within approximately 2,768 held-out genes per fold; when a score cutoff is tied, values denote expected counts under random ordering within the tie. AP and AUROC evaluate observed-positive versus unlabeled label recovery, not true positive versus negative Anopheles phenotypes. The permutation reference averages 4,500 fold-level draws (300 per fold). The simple heuristics were not separately fitted on training labels.')

D.add_heading('Results: top-k label enrichment and novelty boundary',1)
para(D,'On the original held-out folds, the M3 ranking yielded a mean 11.8 observed-positive genes among the top 50, corresponding to proxy precision 23.6%, proxy recall 8.22% of observed positives, and 4.55-fold enrichment over fold-specific prevalence. Random ranking selected 2.59 observed-positive genes on average among 50; the four-feature expression heuristic selected 8.93. These are comparative-label recoveries rather than estimates of newly discovered sterility genes or the number of successful editing targets. At fold-level top 100, M3 captured 21.6 observed positives on average (21.6% proxy precision; 4.16-fold enrichment).')
modeltop=top[top.method=='M3_XGBoost_mean_OOF'].set_index('k')
rows2=[]
for k in [20,50,100,200,500,1000]:
 a=modeltop.loc[k];rows2.append((k,int(a.observed_positive_count),int(a.unlabeled_count),fmt(a.expected_random_positive_count,2),fmt(a.enrichment_over_prevalence,2)+'×'))
D.add_paragraph('Table Y. Genome-wide M3 mean-out-of-fold rankings: comparative-label composition, not functional validation.',style='Caption')
table(D,['Top-k','Observed positives','Unlabeled genes','Random expected positives','Positive-label enrichment'],rows2)
para(D,'The full modelling universe contained 718 observed positives (5.19%) and 13,121 unlabeled genes. Among the top 100 genes ranked solely by mean M3 out-of-fold score, 24 had observed-positive comparative labels and 76 were unlabeled. In contrast, the earlier frozen 199-gene non-benchmark CRISPR-screened integrated ranking comprised 150 observed positives and 49 unlabeled genes; all of its first 100 ranked genes were already observed-positive, with the first unlabeled candidate at non-benchmark rank 112. The two top-100 results therefore arise from different ranking designs and must not be conflated. The latter reflects preceding biological candidate selection and score integration, including phenotype-linked information, and cannot be presented as out-of-sample discovery.')

D.add_heading('Discussion: operational interpretation and limitations',1)
para(D,'The approximately four- to fivefold enrichment among held-out positive annotations quantifies the number of existing comparative labels potentially reviewed per experimental-assessment budget, rather than the yield of genuinely effective target genes. If 50 genes were screened using the cross-fitted proxy ranking, the internally estimated mean observed-positive label count would be about 12 versus approximately 3 from random sampling, but no corresponding estimate of experimentally functional sterility hits can be derived. In addition, the original top-100 non-benchmark shortlist was enriched by an upstream phenotype-informed filter; therefore, its 100% comparative-label-positive composition is not a performance estimate and does not demonstrate novel discovery.')
para(D,'Interpretation remains restricted by non-random Drosophila orthology-dependent positive ascertainment; unlabelled genes are biologically unresolved; the same cross-validation folds influenced model choice; observations are repeatedly assessed across three repeats; and all simple composite weights were selected for a post hoc diagnostic rather than optimized prospectively. The absence of an independent Anopheles-native outcome dataset prevents calculation of genuine sensitivity, specificity, positive predictive value or true recall for newly uncharacterized genes. Claims are accordingly restricted to reproducible refinement of comparative evidence and reference-sequence technical filtering.')
D.add_heading('Figure legends and reproducibility',1)
para(D,'Figure X. Tie-aware average precision of the original M3-XGBoost out-of-fold rankings against randomized-score, single-feature and simple non-PU heuristics. The bar chart uses the 15 original held-out folds with 300 score permutations per fold; no error bars are interpreted as independent-sample confidence intervals. Generate the final publication-resolution PNG using the accompanying base-R script.')
para(D,'Figure Y. Number of observed-positive versus unlabeled genes at increasing top-k thresholds of the genome-wide mean M3 out-of-fold ranking. This is a separate diagnostic from the biologically preselected 199-gene non-benchmark integrated candidate ranking; the plot cannot establish true hits in unlabeled genes.')
para(D,'Reference: Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLoS One. 2015;10(3):e0118432. doi:10.1371/journal.pone.0118432. Full numerical audit, detailed fold-level comparisons, source-compatible metric reproduction and a reproducible base-R script are provided as supporting files. Original PLOS manuscript figure/table numbering must be reconciled before submission.')
D.save(ROOT/'R13_Manuscript_Insert.docx')

R=Document();setup(R,'R1.3 and R2.Q2 point-by-point reviewer response')
R.add_heading('Reviewer 1 — Comment 3: baselines and practical utility',1)
para(R,'We thank the reviewer for requesting formal comparisons that contextualize the reported PU proxy-discrimination results. We reanalysed the exact original M3-XGBoost cross-fitted predictions and verified the reported trapezoidal precision–recall AUC of 0.158535 and AUROC of 0.797601 to numerical precision. The revision includes the same-fold comparison of M3 against 300 within-fold randomized score permutations, individual ovary-expression and ovary-specificity rankings, a four-feature expression composite, and a five-feature expression–orthology composite that does not use FlyBase outcome-derived features. All methods were evaluated over the same 13,839 eligible genes (718 observed positives) and the original three repeats of five held-out folds.')
para(R,'Because simple features generate ties, we applied a single tie-aware average-precision metric across the full comparator panel and retained the legacy trapezoidal metric solely for replication of the initially reported value. M3 achieved mean AP 0.1618 and AUROC 0.7976, versus 0.0545 and 0.5006 under permuted-score ranking, 0.0852 and 0.6140 for ovary-expression alone, 0.0960 and 0.6521 for ovary specificity, and 0.0984 and 0.6183 for the untuned four-expression composite. The corresponding five-feature expression–orthology comparator achieved AP 0.1015. We therefore demonstrate that M3 better recovers available comparative-positive labels than these selected simple baselines, while explicitly avoiding a claim of biological diagnostic accuracy.')
para(R,'At the top-50 held-out threshold, the M3 ranking yielded 11.8 observed-positive labels on average (23.6% proxy precision, 4.55-fold prevalence-relative enrichment), compared with 2.59 under score randomization and 8.93 under the four-expression comparator. This is an operational interpretation of label-enrichment at a finite assay budget, not an estimate of true experimental-hit yield: because the unlabeled group includes untested positives and the labelled positives depend on Drosophila orthology, the probability of successful Anopheles editing or sterility induction cannot be inferred from AUPRC. We have qualified all conclusions accordingly. We also disclose that original model selection used the cross-validation results and that repeated folds do not constitute independent experimental validation.')
R.add_heading('Reviewer 2 — Question 2: novelty, unlabeled genes and Recall@50',1)
para(R,'We agree that “non-benchmark” must not be interpreted as functionally novel or free of reproductive annotations. It indicates exclusion from the six prespecified external Anopheles landmarks. The originally frozen 199 non-benchmark CRISPR-screened genes comprise 150 FlyBase-derived observed positives and 49 unlabeled genes; all 100 genes in the upper integrated shortlist are already observed-positive. The first unlabeled gene appears at non-benchmark rank 112. We have rewritten the interpretation of this upper tier as comparative-evidence refinement and technical filtering, not de novo biological discovery.')
para(R,'To address model utility independently of that selected candidate panel, we report genome-wide mean out-of-fold M3 ranking composition across the 13,839 eligible genes. Its highest 50 genes contain 15 observed positives and 35 unlabeled genes, while its highest 100 contain 24 observed positives and 76 unlabeled genes. These scores demonstrate that the model assigns ranks to unlabeled genes, but do not establish that any unlabeled gene is a true functional reproductive target. A “recall among truly unlabeled positive genes” cannot be calculated from the available data because the true outcome of unlabeled genes is unobserved; treating the entire unlabeled class as confirmed negative or positive would be methodologically misleading. We instead report counts of unlabeled genes and recall restricted explicitly to the observed-positive comparative labels.')
R.add_heading('Evidence status and explicit limitations',1)
para(R,'Status: Baseline calculations, original-metric concordance checks, same-fold top-k diagnostics, complete gene-level provenance tables and point-by-point draft language are completed. R figure-generation and independent R rerun are supplied as a script but have not been executed in the analysis environment. An independently curated set of Anopheles-native fertility outcomes and nested or external model-selection validation remain important for stronger predictive claims; no experiment or unverified outcome is claimed.')
R.add_heading('Files supplied',1)
para(R,'R13_fold_matched_comparison_summary.csv; R13_all_fold_baseline_results_with_random_permutations.csv; R13_full_universe_topk_utility.csv; R13_all_13839_genes_scores_and_comparators.csv; R13_original_199_novel_targeted_label_provenance.csv; R13_49_unlabeled_genes_in_original_targeted_panel.csv; 01_compute_baselines_and_utility.py; 02_reproduce_in_R_and_make_figures.R; R13_Manuscript_Insert.docx; R13_Reviewer_Responses.docx. Future manuscript page, table and figure numbers will be supplied only after integration into the final revised article.')
R.save(ROOT/'R13_Reviewer_Responses.docx')

# Plain text export from full authored docs to bypass download problems.
for stem in ['R13_Manuscript_Insert','R13_Reviewer_Responses']:
 d=Document(ROOT/(stem+'.docx'))
 out=[]
 for p in d.paragraphs:
  if p.text.strip():out.append(p.text)
 for tab in d.tables:
  out.append('\n'+' | '.join(c.text for c in tab.rows[0].cells))
  out += [' | '.join(c.text for c in r.cells) for r in tab.rows[1:]]
 (ROOT/(stem+'.txt')).write_text('\n\n'.join(out),encoding='utf-8')
print('Word documents generated')
