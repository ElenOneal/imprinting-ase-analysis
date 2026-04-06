#!/usr/bin/env bash
set -euo pipefail

# Align hybrid endosperm RNA-seq reads to genome with STAR. This is a basic alignment. For pseudogenome alignments, see pseudogenome.sh.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 07_align_star.sh <samples.tsv> <output_dir> <read_dir> <genome_dir> <genome_fasta> <genome_gtf>
#
# samples.tsv columns: r1  r2  sample_id  genotype  cross


# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 samples output_dir read_dir genome_dir genome gtf"
    exit 1
fi

# Assign arguments to variables
samples="$1"
output_dir="$2"
read_dir="$3"
genome_dir="$4"
genome="$5"
gtf="$6"

# Validate inputs
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r r1 r2 sample_id genotype cross; do
  # Skip comment lines
  [[ "$r1" =~ ^# ]] && continue

  script_file="$output_dir/${sample_id}.star1.sh"
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name=star_align_'"${sample_id}"
    echo '#SBATCH --output='"${output_dir}/${sample_id}.star1.out"
    echo '#SBATCH --error='"${output_dir}/${sample_id}.star1.err"
    echo '#SBATCH --cpus-per-task=6'
    echo '#SBATCH -p common,scavenger'
    echo '#SBATCH --mem=48G'
    echo "#SBATCH --chdir=$output_dir"
    echo ''
    echo "STAR --runThreadN 6 \\"
    echo "     --genomeDir $genome_dir \\"
    echo "     --readFilesCommand gunzip -c \\"
    echo "     --readFilesIn $read_dir/${sample_id}.PE.R1.fq.gz $read_dir/${sample_id}.PE.R2.fq.gz \\"
    echo "     --outSAMtype BAM SortedByCoordinate \\"
    echo "     --outSAMmapqUnique 60 \\"
    echo "     --outFileNamePrefix $output_dir/${sample_id}. \\"
    echo "     --outFilterScoreMinOverLread 0.33 \\"
    echo "     --outFilterMatchNminOverLread 0.33 \\"
    echo "     --outSAMattributes NH HI AS NM nM MD jM jI \\"
    echo "     --sjdbGTFfile $genome_dir/$gtf \\"
    echo "     --outSAMunmapped Within KeepPairs \\"
    echo "     --alignEndsProtrude 10 ConcordantPair \\"
    echo "     --outFilterMismatchNoverReadLmax 0.06"
  } > "$script_file"

  chmod 755 "$script_file"
  echo "Created alignment script: $script_file"
done < "$samples"

echo "All STAR alignment scripts created. To submit jobs, run:"
echo "for script in $output_dir/*.star1.sh; do sbatch \"\$script\"; done"
