"""Gene-level Ag1000G Phase 3 variation summaries for MosqEdit-R.

Uses the maintained malariagen_data Ag3.snp_allele_frequencies() API, which accepts
a gene or region identifier and returns allele frequencies without materialising a
chromosome-scale genotype tensor in memory. The routine is resumable via a cache CSV.
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import malariagen_data


def _feature_gene_table(ag3):
    g = ag3.genome_features(attributes=["ID", "Name", "description"])
    g = g[g["type"].eq("gene")].copy()
    g["gene_id"] = g["ID"].astype(str).str.replace("gene:", "", regex=False)
    g = g[g["gene_id"].str.startswith("AGAP")]
    return g[["gene_id", "contig", "start", "end"]].drop_duplicates("gene_id")


def compute_gene_variation_ag3(annotation_csv, species=("gambiae", "coluzzii"),
                               min_minor_allele_frequency=0.01,
                               cache_csv="data_raw/ag3_gene_variation_cache.csv"):
    ag3 = malariagen_data.Ag3()
    ann = pd.read_csv(annotation_csv, dtype={"gene_id": str})
    genes = _feature_gene_table(ag3)
    work = ann[["gene_id"]].drop_duplicates().merge(genes, on="gene_id", how="left")
    work = work.dropna(subset=["contig", "start", "end"]).copy()

    # Reproducible test mode is useful before a full multi-hour genome-wide run.
    max_genes = int(os.environ.get("MOSQEDIT_MAX_GENES", "0") or 0)
    if max_genes > 0:
        work = work.head(max_genes)

    cache_path = Path(cache_csv)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        cached = pd.read_csv(cache_path, dtype={"gene_id": str})
        done = set(cached["gene_id"].astype(str))
        rows = cached.to_dict("records")
    else:
        done = set()
        rows = []

    taxon_query = "taxon in " + repr(list(species))
    masks = set(getattr(ag3, "site_mask_ids", []) or [])
    site_mask = "gamb_colu" if "gamb_colu" in masks else None

    for i, r in enumerate(work.itertuples(index=False), start=1):
        if r.gene_id in done:
            continue
        length_kb = max((float(r.end) - float(r.start) + 1.0) / 1000.0, 1e-9)
        try:
            af = ag3.snp_allele_frequencies(
                region=r.gene_id,
                cohorts={"combined": taxon_query},
                min_cohort_size=10,
                site_mask=site_mask,
                drop_invariant=True,
                effects=False,
                include_counts=True,
            )
            af = af.reset_index()
            if af.empty:
                row = dict(gene_id=r.gene_id, nucleotide_diversity_pi_proxy=np.nan,
                           snp_density=0.0, high_freq_variant_density=0.0,
                           high_freq_variant_count=0, n_variant_alleles=0,
                           ag3_status="OK_NO_VARIATION")
            else:
                # max_af is documented by malariagen_data and equals the largest cohort
                # alternate-allele frequency. With one custom cohort it is that cohort's AF.
                freq = pd.to_numeric(af.get("max_af"), errors="coerce")
                # Count unique SNP positions where possible; fall back to allele rows.
                pos_col = next((c for c in af.columns if str(c).lower() in {"position", "pos"}), None)
                contig_col = next((c for c in af.columns if str(c).lower() in {"contig", "chrom", "chromosome"}), None)
                if pos_col is not None:
                    site_key = af[pos_col].astype(str)
                    if contig_col is not None:
                        site_key = af[contig_col].astype(str) + ":" + site_key
                    n_sites = int(site_key.nunique())
                    hi_sites = int(site_key[freq.fillna(-1) >= min_minor_allele_frequency].nunique())
                else:
                    n_sites = int(freq.notna().sum())
                    hi_sites = int((freq >= min_minor_allele_frequency).sum())
                # Explicitly a variant-allele heterozygosity proxy, not Nei pi over callable bases.
                pi_proxy = float(np.nanmean(2.0 * freq.to_numpy(float) * (1.0 - freq.to_numpy(float)))) if freq.notna().any() else np.nan
                row = dict(gene_id=r.gene_id, nucleotide_diversity_pi_proxy=pi_proxy,
                           snp_density=n_sites / length_kb,
                           high_freq_variant_density=hi_sites / length_kb,
                           high_freq_variant_count=hi_sites,
                           n_variant_alleles=int(freq.notna().sum()), ag3_status="OK")
        except Exception as exc:
            row = dict(gene_id=r.gene_id, nucleotide_diversity_pi_proxy=np.nan,
                       snp_density=np.nan, high_freq_variant_density=np.nan,
                       high_freq_variant_count=np.nan, n_variant_alleles=np.nan,
                       ag3_status=f"ERROR:{type(exc).__name__}:{str(exc)[:180]}")
        rows.append(row)

        # Frequent checkpoints make the full run resumable.
        if i % 100 == 0:
            pd.DataFrame(rows).drop_duplicates("gene_id", keep="last").to_csv(cache_path, index=False)

    out = pd.DataFrame(rows).drop_duplicates("gene_id", keep="last")
    out.to_csv(cache_path, index=False)
    return out
