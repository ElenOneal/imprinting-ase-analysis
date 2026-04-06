#!/usr/bin/env bash
set -euo pipefail

# Filter hybrid rna alignments to retain only uniquely mapping, properly paired reads with no secondary alignments. These reads will be used to count allele-specific expression.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 10_filter_rna.sh <samples.tsv> <output_dir>
#
# samples.tsv columns: r1  r2  sample_id  genotype  cross

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 samples.tsv output_dir"
    exit 1
fi

# Assign arguments to variables
samples="$1"
output_dir="$2"

# Validate inputs
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

if [ ! -d "$output_dir" ]; then
    echo "Error: Output directory '$output_dir' does not exist."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r r1 r2 sample_id job_id cross; do
  # Skip comment lines
  [[ "$r1" =~ ^# ]] && continue
  
  # Validate BAM file exists
  bam_file="$output_dir/${sample_id}.2.Aligned.sortedByCoord.out.bam"
  if [ ! -f "$bam_file" ]; then
    echo "Error: BAM file '$bam_file' not found. Skipping $sample_id."
    continue
  fi
  
  script_file="$output_dir/${sample_id}.filter.sh"
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name='"${job_id}"
    echo '#SBATCH --output='"${output_dir}/${sample_id}.filter.out"
    echo '#SBATCH --error='"${output_dir}/${sample_id}.filter.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH -p common,scavenger'
    echo "#SBATCH --chdir=$output_dir"
    echo '#SBATCH --mem=6G'
    echo ''
    echo 'source "$(conda info --base)/etc/profile.d/conda.sh"'
    echo 'conda activate imprinting-align'
    echo ''
    echo "# Index original BAM file"
    echo "samtools index ${sample_id}.2.Aligned.sortedByCoord.out.bam"
    echo ''
    echo "# Filter for high-quality, properly paired reads (MAPQ >= 60)"
    echo "samtools view -q 60 -f 0x2 -b ${sample_id}.2.Aligned.sortedByCoord.out.bam > ${sample_id}.f1.bam"
    echo "samtools index ${sample_id}.f1.bam"
    echo ''
    echo "# Filter for uniquely mapping reads (NH:i:1)"
    echo "samtools view -h ${sample_id}.f1.bam | grep -E '(NH:i:1|^@)' | samtools view -Shb > ${sample_id}.f2.bam"
    echo "samtools index ${sample_id}.f2.bam"
  } > "$script_file"

  chmod 755 "$script_file"
  echo "Created filtering script: $script_file"
done < "$samples"

echo "All filtering scripts created. To submit jobs, run:"
echo "for script in $output_dir/*.filter.sh; do sbatch \"\$script\"; done"