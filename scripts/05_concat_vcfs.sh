#!/usr/bin/env bash
set -euo pipefail

# Concatenate per-chromosome VCF files into genome-wide VCFs.
# Usage: bash 05_concat_vcfs.sh <listfile> <prefix> <output_dir>

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 listfile prefix output_dir"
    exit 1
fi

listfile="$1"
prefix="$2"
output_dir="$3"

if [ ! -f "$listfile" ]; then
    echo "Error: File '$listfile' not found."
    exit 2
fi

# Build file lists from the second column (chrname)
filtered_list="${output_dir}/${prefix}.filtered.list"
coding_list="${output_dir}/${prefix}.coding_snps.list"
final_list="${output_dir}/${prefix}.final_snps.list"

while IFS=$'\t' read -r region chrname; do
    [[ "$region" =~ ^# ]] && continue
    echo "${output_dir}/Chr_${chrname}.${prefix}.filtered.vcf.gz"
done < "$listfile" > "$filtered_list"

while IFS=$'\t' read -r region chrname; do
    [[ "$region" =~ ^# ]] && continue
    echo "${output_dir}/Chr_${chrname}.${prefix}.coding_snps.vcf.gz"
done < "$listfile" > "$coding_list"

while IFS=$'\t' read -r region chrname; do
    [[ "$region" =~ ^# ]] && continue
    echo "${output_dir}/Chr_${chrname}.${prefix}.final_snps.vcf.gz"
done < "$listfile" > "$final_list"

# Generate and submit concat job
script_file="${output_dir}/${prefix}.concat.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=${prefix}.concat"
echo "#SBATCH --output=${output_dir}/${prefix}.concat.out"
echo "#SBATCH --error=${output_dir}/${prefix}.concat.err"
echo '#SBATCH --cpus-per-task=1'
echo '#SBATCH -p common,scavenger'
echo "#SBATCH --chdir=${output_dir}"
echo '#SBATCH --mem=14G'
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
echo "rm *.snpcluster.vcf *.snpcluster.vcf.idx"
echo "python ${output_dir}/call_pairwise_snps_w_genes.py ${output_dir}/${prefix}.genic_snps.txt $prefix 40 ${output_dir}/Mguttatusvar_IM767_887_v2.1.gene_exons.gff3"
} > "$script_file"

echo "Submitting ${script_file}..."
sbatch "$script_file"