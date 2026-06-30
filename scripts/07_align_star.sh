#!/usr/bin/env bash
set -euo pipefail

# Align hybrid endosperm RNA-seq reads to genome with STAR. This is a basic alignment. For pseudogenome alignments, see pseudogenome.sh.
# Generates one SLURM job script per SAMPLE (not per lane) -- all of a sample's lanes
# (one, two, or more) are merged via comma-separated --readFilesIn, so samples with
# different lane counts (e.g. a resequenced sample with only one usable lane) are
# handled by the same logic, no special-casing required.
#
# Usage: bash 07_align_star.sh <samples.tsv> <output_dir> <read_dir> <genome_dir> <genome> <gtf> <partition>
# samples.tsv columns: r1  r2  sample_id  barcode  lane  replicate  cross_direction  species  notes
#
# Expects trimmed files at: $read_dir/${sample_id}.${lane}_val_1.fq.gz / _val_2.fq.gz
# (matches the --basename ${sample_id}.${lane} convention used by the Trim Galore script)

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 samples output_dir read_dir genome_dir genome gtf partition"
    exit 1
fi

# Assign arguments to variables
samples="$1"
output_dir="$2"
read_dir="$3"
genome_dir="$4"
genome="$5"
gtf="$6"
partition="$7"

# Validate inputs
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

declare -A r1_list
declare -A r2_list

# Group lanes by sample_id (equivalently, by replicate -- the two are 1:1 in this design)
while IFS=$'\t' read -r r1 r2 sample_id barcode lane replicate cross_direction species notes; do
  r1="${r1#$'\ufeff'}"
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue

  lane_id="${sample_id}.${lane}"
  trimmed_r1="$read_dir/${lane_id}_val_1.fq.gz"
  trimmed_r2="$read_dir/${lane_id}_val_2.fq.gz"

  if [[ -z "${r1_list[$sample_id]:-}" ]]; then
    r1_list[$sample_id]="$trimmed_r1"
    r2_list[$sample_id]="$trimmed_r2"
  else
    r1_list[$sample_id]+=",$trimmed_r1"
    r2_list[$sample_id]+=",$trimmed_r2"
  fi
done < "$samples"

# One alignment job per sample, using all of that sample's lanes
for sample_id in "${!r1_list[@]}"; do
  script_file="$output_dir/${sample_id}.star1.sh"
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name=star_'"${sample_id}"
    echo '#SBATCH --output='"${output_dir}/${sample_id}.star1.out"
    echo '#SBATCH --error='"${output_dir}/${sample_id}.star1.err"
    echo '#SBATCH --cpus-per-task=6'
    echo "#SBATCH -p $partition"
    echo '#SBATCH --mem=48G'
    echo "#SBATCH --chdir=$output_dir"
    echo ''
    echo "source $(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "STAR --runThreadN 6 \\"
    echo "     --genomeDir $genome_dir \\"
    echo "     --readFilesCommand gunzip -c \\"
    echo "     --readFilesIn ${r1_list[$sample_id]} ${r2_list[$sample_id]} \\"
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
done

echo "All STAR alignment scripts created (one per sample). To submit jobs, run:"
echo "for script in $output_dir/*.star1.sh; do sbatch \"\$script\"; done"
