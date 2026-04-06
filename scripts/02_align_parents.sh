#!/usr/bin/env bash
set -euo pipefail

# Align parental DNA-seq reads to reference genome with BWA.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 02_align_parents.sh <samples.tsv> <output_dir> <clean_dir> <genome_dir> <genome_fasta>
#
# samples.tsv columns: r1  r2  sample_id  genotype  barcode

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 samples.tsv output_dir clean_dir genome_dir genome_fasta"
    exit 1
fi


# Assign arguments to variables
samples="$1"
output_dir="$2"
clean_dir="$3"
genome_dir="$4"
genome="$5"

# Check if the input file exists
if [ ! -f "$samples" ]; then
    echo "Error: sample sheet '$samples' not found."
    exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Main loop
while IFS=$'\t' read -r f c a b d; do
  [[ "$f" =~ ^#  ]] && continue
  script_file="$output_dir/$a.align.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.align.out"
    echo "#SBATCH --error=${a}.align.err"
    echo '#SBATCH --cpus-per-task=6'
    echo '#SBATCH -p common,scavenger'
    echo "#SBATCH --chdir=$output_dir"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "source $(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "bwa mem -t 6 -M $genome_dir/$genome $clean_dir/$a.PE.R1.fq.gz $clean_dir/$a.PE.R2.fq.gz | samtools view -Shb | samtools sort -T $a.sort -o $a.sort.bam"
    echo "samtools index $a.sort.bam"
    echo "sbatch $output_dir/${a}.filter.sh"
  } > "$script_file"
done < "$samples"

echo "All alignment scripts created. To submit jobs, run:"
echo "for script in $output_dir/*.align.sh; do sbatch \"\$script\"; done"