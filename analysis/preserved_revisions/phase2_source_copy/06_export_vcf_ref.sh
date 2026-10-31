#!/usr/bin/env bash
# Export REF allele at Phase 2 PASS SNP positions. No genotype files downloaded.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output
out="output/AR1_PASS_SNP_REF_5arms.tsv"
if ! command -v bcftools >/dev/null 2>&1; then echo 'ERROR: bcftools is not installed' >&2; exit 2; fi
printf 'chromosome\tpos\tref\talt\n' > "${out}.tmp"
expected=(727 395 239 576 69)
arms=(2L 2R 3L 3R X)
for i in "${!arms[@]}"; do
  arm="${arms[$i]}"
  f="output/AR1_frozen_${arm}_genotypes.vcf.gz"
  if [[ ! -s "$f" ]]; then echo "STOP: missing $f" >&2; rm -f "${out}.tmp"; exit 2; fi
  echo "Checking $arm ..."
  bcftools view -h "$f" | grep -Ei '^##(reference|source|contig)=' | head -n 5 || true
  # Query all REF alleles, even if ALT is multiallelic. Comparison does not use genotypes.
  n=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$f" | tee -a "${out}.tmp" | wc -l)
  echo "  $arm: $n variants"
  if [[ "$n" != "${expected[$i]}" ]]; then
      echo "STOP: expected ${expected[$i]} rows for $arm; obtained $n" >&2
      rm -f "${out}.tmp"; exit 2
  fi
done
n=$(($(wc -l < "${out}.tmp") - 1))
[[ "$n" == "2006" ]] || { echo "STOP: expected 2006 variant rows, got $n" >&2; exit 2; }
mv "${out}.tmp" "$out"
echo "Exported $n PASS variant REF alleles: $out"
echo 'This is a partial reference check at variant positions only; full FASTA validation remains necessary.'
