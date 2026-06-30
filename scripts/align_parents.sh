#!/usr/bin/env bash
set -euo pipefail

# Align parental DNA-seq reads to reference genome with BWA.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 02_align_parents.sh <samples.tsv> <execute_dir> <reads_dir> <genome_dir> <genome_fasta> <partition> [time]
#
# samples.tsv columns: r1  r2  sample_id  genotype  barcode population species

# Check if the correct number of arguments is provided
if [ "$#" -lt 6 ] || [ "$#" -gt 7 ]; then
  echo "Usage: $0 samples.tsv execute_dir reads_dir genome_dir genome_fasta partition [time]"
    exit 1
fi


# Assign arguments to variables
samples="$1"
execute_dir="$2"
reads_dir="$3"
genome_dir="$4"
genome="$5"
partition="$6"
time_limit="${7:-24:00:00}"

# Check if the input file exists
if [ ! -f "$samples" ]; then
    echo "Error: sample sheet '$samples' not found."
    exit 2
fi

# Create execute directory if it doesn't exist
mkdir -p "$execute_dir"

# Main loop
while IFS=$'\t' read -r r1 r2 sample_id barcode _; do
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue
  script_file="${sample_id}.align.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=${sample_id}.align"
    echo "#SBATCH --output=${sample_id}.align.out"
    echo "#SBATCH --error=${sample_id}.align.err"
    echo '#SBATCH --cpus-per-task=6'
    echo "#SBATCH -p $partition"
    echo "#SBATCH --time=$time_limit"
    echo "#SBATCH --chdir=$execute_dir"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "source $(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "bwa mem -t 6 -M $genome_dir/$genome $reads_dir/${sample_id}.PE.R1.fq.gz $reads_dir/${sample_id}.PE.R2.fq.gz | samtools view -Shb | samtools sort -T ${sample_id}.sort -o ${sample_id}.sort.bam"
    echo "samtools index ${sample_id}.sort.bam"
  } > "$script_file"
done < "$samples"

echo "All alignment scripts created. To submit jobs, run:"
echo "for script in *.align.sh; do sbatch \"\$script\"; done"