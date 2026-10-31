#!/usr/bin/env bash
# Usage: bash 02_extract_region.sh 2R pilot
#        bash 02_extract_region.sh 2R full
# Requirements: WSL Ubuntu with bcftools + tabix installed; HTTPS Range support.
set -euo pipefail
cd "$(dirname "$0")"
chrom="${1:-2R}"
mode="${2:-pilot}"
case "$chrom" in 2L|2R|3L|3R|X) ;; *) echo 'Allowed chromosome arms: 2L 2R 3L 3R X' >&2; exit 2;; esac
case "$mode" in pilot|full) ;; *) echo 'Choose pilot or full' >&2; exit 2;; esac
command -v bcftools >/dev/null || { echo 'Install bcftools first: sudo apt-get update && sudo apt-get install -y bcftools tabix' >&2; exit 3; }
python3 01_make_regions.py
base='https://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/variation/main/vcf/pass'
url="$base/ag1000g.phase2.ar1.pass.${chrom}.vcf.gz"
mkdir -p output
if [ "$mode" = pilot ]; then
    IFS=$'\t' read -r c start0 end0 < "regions/${chrom}_merged_0based.bed"
    start=$((start0+1)); end=$((end0+10000))
    # Extend pilot around first selected site for likelihood of a variant, rather than assuming one.
    region="${chrom}:${start}-${end}"
    outfile="output/PILOT_${chrom}.vcf.gz"
    echo "PILOT: querying remote indexed VCF in $region (not downloading chromosome)."
    bcftools view -r "$region" -Oz -o "$outfile" "$url"
else
    bed="regions/${chrom}_merged_0based.bed"
    outfile="output/AR1_frozen_${chrom}_genotypes.vcf.gz"
    echo "FULL TARGETED: $(wc -l < "$bed") merged intervals from the frozen panel on $chrom"
    bcftools view -R "$bed" -Oz -o "$outfile" "$url"
fi
bcftools index -t -f "$outfile"
echo "OUTPUT: $outfile"
echo "VCF chromosome names:"; bcftools index -s "$outfile" || true
echo "Individual VCF sample columns: $(bcftools query -l "$outfile" | wc -l)"
echo "Variant records extracted: $(bcftools view -H "$outfile" | wc -l)"
echo 'IMPORTANT: Missing variants do not prove invariant target sequences. Download and assess accessibility separately.'
