#!/usr/bin/env bash
set -euo pipefail

# Realign hybrid endosperm RNA-seq reads to genome with STAR. This will incorporate novel junctions discovered by STAR in the first alignment.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 09_realign_star.sh <samples.tsv> <output_dir> <read_dir> <genome_dir> <genome_fasta> <genome_gtf>
#
# samples.tsv columns: r1  r2  sample_id  job_id  cross

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 samples.tsv output_dir read_dir genome_dir genome gtf"
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

if [ ! -d "$output_dir" ]; then
    echo "Error: Output directory '$output_dir' does not exist."
    exit 2
fi

if [ ! -d "$read_dir" ]; then
    echo "Error: Read directory '$read_dir' does not exist."
    exit 2
fi

if [ ! -d "$genome_dir" ]; then
    echo "Error: Genome directory '$genome_dir' does not exist."
    exit 2
fi

if [ ! -f "$genome_dir/$genome" ]; then
    echo "Error: Genome file '$genome_dir/$genome' does not exist."
    exit 2
fi

if [ ! -f "$genome_dir/$gtf" ]; then
    echo "Error: GTF file '$genome_dir/$gtf' does not exist."
    exit 2
fi

if [ ! -f "$genome_dir/spliced.tab" ]; then
    echo "Error: Splice junction file '$genome_dir/spliced.tab' does not exist."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r r1 r2 sample_id job_id cross; do
  # Skip comment lines
  [[ "$r1" =~ ^# ]] && continue

  script_file="$output_dir/${sample_id}.star2.sh"
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name='"${job_id}"
    echo '#SBATCH --output='"${output_dir}/${sample_id}.star2.out"
    echo '#SBATCH --error='"${output_dir}/${sample_id}.star2.err"
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
    echo "     --outFileNamePrefix $output_dir/${sample_id}.2. \\"
    echo "     --outFilterScoreMinOverLread 0.33 \\"
    echo "     --outFilterMatchNminOverLread 0.33 \\"
    echo "     --alignEndsProtrude 10 ConcordantPair \\"
    echo "     --outSAMattributes NH HI AS NM nM MD jM jI \\"
    echo "     --sjdbGTFfile $genome_dir/$gtf \\"
    echo "     --sjdbFileChrStartEnd $genome_dir/spliced.tab \\"
    echo "     --outSAMunmapped Within KeepPairs \\"
    echo "     --outFilterMismatchNoverReadLmax 0.06"
  } > "$script_file"

  chmod 755 "$script_file"
  echo "Created realignment script: $script_file"
done < "$samples"

echo "All STAR realignment scripts created. To submit jobs, run:"
echo "for script in $output_dir/*.star2.sh; do sbatch \"\$script\"; done"
