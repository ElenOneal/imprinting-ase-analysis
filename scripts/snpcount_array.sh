#!/usr/bin/env bash
set -euo pipefail

# Count parental allele-specific reads per gene using a SLURM array job.
# Parallelizes by chromosome. For each sample, generates and submits one
# array job script with one task per chromosome in the SNP catalog.
#
# Usage: bash 11_snpcount_array.sh <samples.tsv> <snp_catalog> <bam_dir> <output_dir> <scripts_dir> <partition> <conda_env>
#
# samples.tsv columns: r1  r2  sample_id  job_id  cross
# snp_catalog: output of call_pairwise_snps_w_genes.py (_coding_sequence_snps.bed)
# bam_dir: directory containing filtered BAM files (*.f2.bam)
# output_dir: directory for output files
# scripts_dir: directory containing classify_parental_alleles.py
# partition: SLURM partition (e.g., common,scavenger)
# conda_env: Conda environment name (e.g., imprinting-align)

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 samples snp_catalog bam_dir output_dir scripts_dir partition conda_env"
    exit 1
fi

samples="$1"
snp_catalog="$2"
bamdir="$3"
output_dir="$4"
scripts_dir="$5"
partition="$6"
conda_env="$7"

# Validate inputs
if [ ! -f "$samples" ]; then
    echo "Error: Sample file '$samples' not found."
    exit 2
fi

if [ ! -f "$snp_catalog" ]; then
    echo "Error: SNP catalog '$snp_catalog' not found."
    exit 2
fi

if [ ! -d "$bamdir" ]; then
    echo "Error: BAM directory '$bamdir' not found."
    exit 2
fi

if [ ! -f "$scripts_dir/classify_parental_alleles.py" ]; then
    echo "Error: classify_parental_alleles.py not found in '$scripts_dir'."
    exit 2
fi

mkdir -p "$output_dir"

# Get unique chromosome list from SNP catalog (skip header)
chr_list=($(tail -n +2 "$snp_catalog" | cut -f1 | sort -u))
nchr=${#chr_list[@]}

if [ "$nchr" -eq 0 ]; then
    echo "Error: No chromosomes found in '$snp_catalog'."
    exit 2
fi

echo "Found ${nchr} chromosomes in SNP catalog."

# Main loop — one array job script per sample
while IFS=$'\t' read -r r1 r2 sample_id job_id cross; do
    [[ "$r1" =~ ^# ]] && continue

    script_file="${output_dir}/${sample_id}.snpcounts.sh"

    {
        echo '#!/bin/bash'
        echo "#SBATCH --job-name=${job_id}"
        echo "#SBATCH --error=${output_dir}/${sample_id}.%a.counts.out"
        echo "#SBATCH --error=${output_dir}/${sample_id}.%a.counts.err"
        echo '#SBATCH --cpus-per-task=1'
        echo "#SBATCH -p $partition"
        echo "#SBATCH --chdir=${output_dir}"
        echo "#SBATCH --array=1-${nchr}"
        echo '#SBATCH --mem=6G'
        echo ''
        echo 'source "$(conda info --base)/etc/profile.d/conda.sh"'
        echo "conda activate $conda_env"
        echo "# Get chromosome name for this array task"
        echo "chr_list=(${chr_list[@]})"
        echo 'chr_name="${chr_list[$((SLURM_ARRAY_TASK_ID-1))]}"'
        echo ''
        echo "# Extract SNPs for this chromosome (3-col BED for samtools, 6-col for python)"
        echo "awk -v chr=\"\$chr_name\" '\$1==chr {print \$1\"\t\"\$2\"\t\"\$3}' ${snp_catalog} > ${output_dir}/${sample_id}.\${chr_name}.bed"
        echo "awk -v chr=\"\$chr_name\" '\$1==chr' ${snp_catalog} > ${output_dir}/${sample_id}.\${chr_name}.snpbed"
        echo ''
        echo "# Subset BAM to SNP regions"
        echo "samtools view -bh ${bamdir}/${sample_id}.f2.bam -L ${output_dir}/${sample_id}.\${chr_name}.bed -o ${output_dir}/${sample_id}.\${chr_name}.coding.bam"
        echo "samtools index ${output_dir}/${sample_id}.\${chr_name}.coding.bam"
        echo ''
        echo "# Count parental allele-specific reads"
        echo "python ${scripts_dir}/classify_parental_alleles.py ${output_dir}/${sample_id}.\${chr_name}.coding.bam ${output_dir}/${sample_id}.\${chr_name}.snpbed ${output_dir}/${sample_id}.\${chr_name}"
        echo ''
        echo "# Remove intermediate files"
        echo "rm ${output_dir}/${sample_id}.\${chr_name}.bed"
        echo "rm ${output_dir}/${sample_id}.\${chr_name}.snpbed"
        echo "rm ${output_dir}/${sample_id}.\${chr_name}.coding.bam"
        echo "rm ${output_dir}/${sample_id}.\${chr_name}.coding.bam.bai"
    } > "$script_file"

    chmod 755 "$script_file"
    echo "Created filtering script: $script_file"
done < "$samples"

echo "All counting scripts created. To submit jobs, run:"
echo "for script in $output_dir/*.snpcounts.sh; do sbatch \"\$script\"; done"