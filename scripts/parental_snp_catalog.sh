#!/usr/bin/env bash
set -euo pipefail

# Concatenate per-chromosome VCF files into genome-wide VCFs.
# Usage: bash 05_snp_catalog.sh <listfile> <prefix> <output_dir> <gene_gff> <partition>

if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
    echo "Usage: $0 listfile prefix output_dir gene_gff partition [mem] [time]"
    exit 1
fi

listfile="$1"
prefix="$2"
output_dir="$3"
gene_gff="$4"
partition="$5"
mem="${6:-14G}"
time_limit="${7:-24:00:00}"


if [ ! -f "$listfile" ]; then
    echo "Error: File '$listfile' not found."
    exit 2
fi

# Validate list file format: should have exactly 1 column (region), no headers
if ! awk 'NF == 0 || $1 ~ /^#/ {next} NF != 1 {exit 1}' "$listfile"; then
    echo "Error: List file '$listfile' must have exactly 1 column (region) with no headers."
    exit 2
fi


# Build file lists from the chromosome name before the ':' in each region
filtered_list="${output_dir}/${prefix}.filtered.list"
coding_list="${output_dir}/${prefix}.coding_snps.list"
final_list="${output_dir}/${prefix}.final_snps.list"

while IFS= read -r region; do
    [[ "$region" =~ ^# ]] && continue
    chrname="${region%%:*}"
    echo "${output_dir}/${chrname}.${prefix}.filtered.vcf.gz"
done < "$listfile" > "$filtered_list"

while IFS= read -r region; do
    [[ "$region" =~ ^# ]] && continue
    chrname="${region%%:*}"
    echo "${output_dir}/${chrname}.${prefix}.coding_snps.vcf.gz"
done < "$listfile" > "$coding_list"

while IFS= read -r region; do
    [[ "$region" =~ ^# ]] && continue
    chrname="${region%%:*}"
    echo "${output_dir}/${chrname}.${prefix}.final_snps.vcf.gz"
done < "$listfile" > "$final_list"

# Generate and submit concat job
script_file="${output_dir}/${prefix}.concat.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=${prefix}.concat"
echo "#SBATCH --output=${output_dir}/${prefix}.concat.out"
echo "#SBATCH --error=${output_dir}/${prefix}.concat.err"
echo '#SBATCH --cpus-per-task=1'
echo "#SBATCH -p $partition"
echo "#SBATCH --chdir=${output_dir}"
echo "#SBATCH --mem=$mem"
echo "#SBATCH --time=$time_limit"
echo ''
echo "source /hpc/group/willislab/eo22/miniforge3/etc/profile.d/conda.sh"
echo "conda activate imprinting-align"
echo ''
echo "bcftools concat -f ${filtered_list} -Oz -o ${output_dir}/${prefix}.filtered.vcf.gz"
echo "tabix ${output_dir}/${prefix}.filtered.vcf.gz"
echo ''
echo "bcftools concat -f ${coding_list} -Oz -o ${output_dir}/${prefix}.coding_snps.vcf.gz"
echo "tabix ${output_dir}/${prefix}.coding_snps.vcf.gz"
echo ''
echo "bcftools concat -f ${final_list} -Oz -o ${output_dir}/${prefix}.final_snps.vcf.gz"
echo "tabix ${output_dir}/${prefix}.final_snps.vcf.gz"
echo ''
echo "echo 'Concatenation complete.'"
echo "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MQ[\t%GT:%DP]' ${output_dir}/${prefix}.final_snps.vcf.gz > ${output_dir}/${prefix}.genic_snps.txt"
echo "python ${output_dir}/call_pairwise_snps_w_genes.py ${output_dir}/${prefix}.genic_snps.txt $prefix 40 ${gene_gff}"
chmod +x "$script_file"
} > "$script_file"

echo "Created script: $script_file"
echo "To submit jobs sbatch $script_file"