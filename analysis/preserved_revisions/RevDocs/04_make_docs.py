from pathlib import Path
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_LINE_SPACING
import pandas as pd
root=Path(__file__).resolve().parent
out=root/'results'
w=pd.read_csv(out/'R2Q3_step04_domain_weight_sensitivity.csv')
d=pd.read_csv(out/'R2Q3_expression_masking_summary.csv')
c=pd.read_csv(out/'R2Q3_expression_coverage_by_cohort.csv')

def doc(title,subtitle):
    x=Document(); sec=x.sections[0];sec.top_margin=Inches(.8);sec.bottom_margin=Inches(.75);sec.left_margin=Inches(.85);sec.right_margin=Inches(.85)
    normal=x.styles['Normal'];normal.font.name='Times New Roman';normal.font.size=Pt(12); normal.paragraph_format.line_spacing=2;normal.paragraph_format.space_after=Pt(4)
    for nm,sz in [('Title',16),('Heading 1',13),('Heading 2',12)]:
        st=x.styles[nm];st.font.name='Times New Roman';st.font.size=Pt(sz);st.font.bold=True;st.font.color.rgb=RGBColor(26,50,78)
        st.paragraph_format.line_spacing=1.2;st.paragraph_format.space_before=Pt(10);st.paragraph_format.space_after=Pt(5)
    x.add_heading(title,0);p=x.add_paragraph(subtitle);p.style='Subtitle'
    h=sec.header.paragraphs[0];h.text='MosqEdit-R | R2.Q3 expression robustness | Evidence-based manuscript drafting';h.style='Normal';h.paragraph_format.line_spacing=1;h.runs[0].font.size=Pt(9)
    ft=sec.footer.paragraphs[0];ft.text='Analysis date: 31 October 2026 | External RNA-seq comparison pending';ft.paragraph_format.line_spacing=1;ft.alignment=WD_ALIGN_PARAGRAPH.CENTER;ft.runs[0].font.size=Pt(9)
    return x

def para(x,txt): x.add_paragraph(txt)
def table(x,headers,rows,widths=None):
    t=x.add_table(rows=1,cols=len(headers));t.style='Light Shading Accent 1';t.alignment=WD_TABLE_ALIGNMENT.CENTER
    for cell,h in zip(t.rows[0].cells,headers): cell.text=str(h)
    for row in rows:
        cells=t.add_row().cells
        for cell,val in zip(cells,row):cell.text=str(val)
    for row in t.rows:
        for cell in row.cells:
            cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cell.paragraphs:
                p.paragraph_format.line_spacing=1.1;p.paragraph_format.space_after=Pt(1)
                for rr in p.runs:rr.font.name='Times New Roman';rr.font.size=Pt(9)
    x.add_paragraph('')

m=doc('Expression coverage and ranking sensitivity','Insertion-ready Methods, Results, and discussion text | Reviewer 2, Question 3')
m.add_heading('Methods — Original expression provenance and missingness',level=1)
para(m,'Adult tissue-expression evidence originated from the MozAtlas microarray series GSE21689 (Affymetrix GPL1321, first published in 2011), with 60 samples representing 15 adult tissue/sex groups and four replicates per group. The original analysis mapped unique probe annotations to the AgamP4 gene universe and retained expression evidence for 9,228 of 13,845 genes (66.65%). Expression was used in the preliminary Female Reproductive Importance Score (FRIS), whereas functional evidence was derived separately from comparative orthology-linked phenotypes. Genes without an accepted expression mapping were treated as having unavailable expression evidence rather than zero expression.')
m.add_heading('Methods — Executed expression-domain sensitivity analysis',level=1)
para(m,'Using the frozen Step 04 output, the originally reported preliminary score was regenerated from the stored FRIS and functional-evidence components by renormalizing weights over available domains, then multiplying by the original completeness adjustment [0.75 + 0.25 × (available domain weight)]. The primary expression:functional weighting was 0.60:0.40. Expression weights of 0, 0.15, 0.30, 0.45, 0.75, 0.90, and 1.00 were analyzed without modifying the source features, labels or candidate-selection process. Rank concordance was quantified over genes ranked in both conditions, and top-k overlap was measured against the original Step 04 preliminary ranking for k = 20, 50, 100, and 200. The original Step 04 score and rank were reproduced to numerical precision (maximum absolute score difference <1 × 10⁻¹⁵).')
para(m,'An independent missing-expression availability stress test randomly masked 10%, 25%, or 50% of genes originally carrying FRIS values. Each fraction was assessed with 150 reproducible repetitions (fixed pseudorandom seed), preserving other evidence and recalculating the original completeness-adjusted preliminary score. Median, 5th and 95th percentiles of top-k overlap were reported. This is an illustrative missing-data perturbation, not a model refit, a technical-error model, a replacement RNA-seq expression profile or empirical evidence of temporal stability.')
m.add_heading('Results — Expression evidence is incomplete and non-randomly selected',level=1)
para(m,'MozAtlas expression was mapped for 9,228/13,845 genes (66.65%). Among the 10,925 genes carrying at least one of the two primary biological evidence domains, 6,840 had both expression and comparative functional evidence, 2,388 expression alone, and 1,697 comparative functional evidence alone. A further 2,920 genes had neither domain. The original top 100 preliminary-ranking genes all had expression evidence. In contrast, 94/100 of the highest-ranked non-benchmark genes in the downstream integrated shortlist were expression-mapped; all 100 of that shortlist carried comparative functional evidence. The enriched expression coverage within the original shortlist reflects the staged selection design and should not be interpreted as representative of all Anopheles genes.')
m.add_heading('Results — Ranking changes under alternative domain weights',level=1)
para(m,'The preliminary top-100 overlap was 93/100 when expression weight was reduced from 0.60 to 0.45, 87/100 at 0.30, 75/100 at 0.75 and 52/100 at 0.90. Removing the expression domain entirely retained 15/100 genes in the original preliminary top 100; removing comparative functional evidence retained 25/100. Thus, the preliminary shortlist was moderately insensitive to limited changes around the original weights but sensitive to major shifts in expression availability or domain emphasis. The analyses only concern the frozen preliminary ranking, and are not evidence that a newly refitted PU model or final CRISPR-integrated ranking is invariant to RNA-seq replacement.')
m.add_heading('Table 1. Step 04 sensitivity to expression-domain weighting',level=2)
table(m,['Expression weight','Functional weight','Spearman ρ','Top 20 overlap','Top 100 overlap'],[[f'{100*row.expression_weight:.0f}%',f'{100*row.functional_weight:.0f}%',f'{row.spearman_common_with_original:.3f}',f'{int(row.top20_overlap_n)}/20',f'{int(row.top100_overlap_n)}/100'] for row in w.itertuples()])
m.add_heading('Results — Missing-expression availability perturbations',level=1)
para(m,'Across 150 independent masking replicates per fraction, the median top-100 overlap was 90/100 (5th–95th percentile, 85–94) when 10% of mapped expression profiles were hidden, 75/100 (69–82) when 25% were hidden, and 50/100 (44–58) when 50% were hidden. These patterns indicate that broad loss of expression availability can materially affect preliminary candidate ordering, even when orthology-derived functional evidence is retained. Because masks were random, the analysis does not address biological platform bias or structured missingness caused by probe design, developmental stage, infection or tissue context.')
m.add_heading('Table 2. Random masking of originally available expression evidence',level=2)
table(m,['Expression masked','Repetitions','Median top 20','Median top 100','Top 100 p05–p95'],[[f'{100*row.fraction_of_available_expression_masked:.0f}%',int(row.n_replicates),f'{int(row.top20_overlap_n_median)}/20',f'{int(row.top100_overlap_n_median)}/100',f'{int(row.top100_overlap_n_p05)}–{int(round(row.top100_overlap_n_p95))}'] for row in d.itertuples()])
m.add_heading('External RNA-seq comparator: prespecified but not yet executed',level=1)
para(m,'A suitable contemporary partial comparator was identified in the 2023 public RNA-seq series GSE239600 (mosGILT study), comprising four wild-type ovary samples, two wild-type fat-body samples and three wild-type whole-body samples at 24 hours after blood feeding. This series supplies newer ovarian transcriptomic evidence but lacks the matched male, testis and multi-tissue coverage required to reconstruct the original 15-group FRIS or sex-bias predictors. Accordingly, the primary external comparison, once the processed expression file has been obtained and identifier/sample mapping verified, will use only within-study gene abundance ranks and matched-AGAP cross-study correlation for ovarian expression. No RNA-seq cross-platform concordance or ranking robustness result is asserted here, because the GSE239600 processed expression file was not available in the executed analysis environment.')
m.add_heading('Interpretation and limitations',level=1)
para(m,'MozAtlas is a historically important adult tissue atlas, but its vintage, incomplete probe-to-gene coverage and specific physiological conditions may limit portability to other developmental stages, infection states and natural populations. The sensitivity analyses demonstrate quantitative dependence on the original expression domain and its missingness, while highlighting the relative stability under modest changes to domain weighting. They do not remove the need for contemporary matched RNA-seq measurements, recalibration or re-fitting of the PU models, or independent functional validation. Future expression data should be analyzed without redefining the previously frozen genomic population-validation panel or selecting the RNA-seq samples based on their effects on candidate ranks.')
m.add_heading('Figure legends (figures generated separately by base-R script)',level=1)
para(m,'Figure X. Effect of expression-domain weight on retention of the original preliminary top-100 genes. Values represent the unchanged Step 04 FRIS and comparative functional domains under prespecified alternative weightings. Dashed line indicates the original 60% expression weight. Not a PU refit or independent RNA-seq analysis.')
para(m,'Figure Y. Step 04 top-100 retention after random removal of 10%, 25% and 50% of originally observed expression profiles. Each condition was evaluated across 150 fixed-seed permutations; boxplots show the distribution of original top-100 overlap under missing-evidence reweighting. This is a designed availability perturbation, not biological expression variability.')
m.add_heading('References',level=1)
para(m,'1. Baker DA, Nolan T, Fischer B, Pinder A, Crisanti A, Russell S. A comprehensive gene expression atlas of sex- and tissue-specificity in the malaria vector, Anopheles gambiae. BMC Genomics. 2011;12:296. doi:10.1186/1471-2164-12-296.')
para(m,'2. GEO accession GSE21689. MozAtlas: A comprehensive gene expression atlas of sex- and tissue-specificity in the malaria vector, Anopheles gambiae. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21689')
para(m,'3. GEO accession GSE239600. Anopheles gambiae mosGILT regulates innate immune genes and zpg expression. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239600')
m.save(root/'R2Q3_Manuscript_Insert.docx')

r=doc('Reviewer 2, Question 3 — Draft response','Executed original-data sensitivity, with an explicit boundary around unexecuted external RNA-seq validation')
r.add_heading('Comment (substantive summary)',level=1)
para(r,'The 2011 MozAtlas/GSE21689 microarray is dated and retains only a subset of gene probes, with incomplete developmental and infection coverage. The reviewer asks why newer RNA-seq resources were not used and requests sensitivity analysis to demonstrate robustness.')
r.add_heading('Response',level=1)
para(r,'We thank the reviewer for identifying the risk of expression-source dependence. We have now quantified both the evidence coverage and the sensitivity of the original preliminary biological ranking to the weighting and availability of the MozAtlas component. We explicitly distinguish these analyses from independent RNA-seq replication and from any reassessment of experimental editing performance.')
para(r,'The Step 04 scoring formula, including its 60:40 expression:functional weighting and its evidence-domain completeness adjustment, was recovered from the frozen source workflow and numerically reproduced against all 13,845 original genes. The original expression data mapped to 9,228 loci (66.65%); 2,920 loci lacked both expression and comparative functional evidence. The highest-ranked non-benchmark integrated shortlist had a much higher expression-mapping rate (94/100), reflecting the staged candidate-selection design rather than complete genome-wide coverage.')
para(r,'In an executed parameter sensitivity analysis, reducing expression weight from 60% to 45% retained 93 of the original preliminary top 100 genes; 30% expression weight retained 87. At 75% expression weight, 75 were retained. Extreme ablation had much greater effects: only 15 of the original top 100 remained under functional-only scoring, and 25 under expression-only scoring. A separate 150-replicate missingness experiment retained a median of 90, 75 and 50 original top-100 genes after masking 10%, 25% and 50%, respectively, of the available expression evidence. These experiments were performed on the preliminary biological score, not by re-fitting the PU model or recalculating the final integrated CRISPR rankings.')
para(r,'We also identified a newer RNA-seq comparator, GEO GSE239600 (publicly released in 2023), with wild-type ovary, fat-body and whole-body samples taken after blood feeding. We have documented its design and prepared a separate identifier-checked analysis for a cross-platform ovarian-expression concordance check. However, its processed expression table was not available for analysis in the current execution environment, and the study lacks the matched male/testis and complete tissue coverage required to reconstruct the full MozAtlas-derived predictors. We therefore do not claim that independent RNA-seq robustness has already been demonstrated. If the processed dataset can be acquired and matched, this external concordance check will be added as an additional limited sensitivity analysis before resubmission.')
para(r,'The manuscript text now quantifies incomplete expression coverage, describes the executed rank sensitivity under weight and availability perturbations, and restricts conclusions to the conditional stability of the original preliminary ranking. It explicitly states that perturbations of summary evidence scores do not replace contemporary molecular measurements or refitting of the supervised PU learners.')
r.add_heading('Material provided',level=1)
table(r,['Deliverable','Contents'],[
 ('Analysis script','01_audit_sensitivity.py: exact score replication, eight weight scenarios, 450 missingness replicates'),
 ('Audited numerical outputs','Coverage, gene-rank sensitivity, masking replicates, effect summaries, source-data checksums'),
 ('Publication figures','02_reproduce_tables_figures.R: two 600-dpi R figures'),
 ('External validation protocol','03_external_GSE239600_RNAseq_comparison.R: source sample and gene-ID guards'),
 ('Manuscript insertion','Methods, Results, two tables, figure legends, references and limitations')])
r.add_heading('Outstanding items before final response',level=1)
para(r,'The GSE239600 processed expression dataset must be obtained and the external cross-platform comparison executed, or its omission explicitly retained as a limitation. If stronger claims of full-model robustness are desired, the original PU models must be refitted with alternative or updated expression features using the same label provenance, nested selection and held-out assessment. The final manuscript and supplementary cross-references remain to be reconciled.')
r.save(root/'R2Q3_Reviewer_Response.docx')

# plain-text mirrors for download fallback
for source in ('R2Q3_Manuscript_Insert.docx','R2Q3_Reviewer_Response.docx'):
    x=Document(root/source)
    lines=[]
    for child in x.element.body:
        from docx.oxml.ns import qn
        if child.tag==qn('w:p'):
            from docx.text.paragraph import Paragraph
            p=Paragraph(child,x);lines.append(p.text)
        elif child.tag==qn('w:tbl'):
            from docx.table import Table
            t=Table(child,x); lines += ['\t'.join(c.text for c in row.cells) for row in t.rows]
    (root/source.replace('.docx','.txt')).write_text('\n'.join(lines),encoding='utf-8')
print('Wrote two Word documents and plain-text mirrors')
