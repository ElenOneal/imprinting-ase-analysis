#!/usr/bin/env bash
set -euo pipefail

# Filter parental bam files with a combination of Picard and samtools.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 03_filter_parent_bams.sh <samples.tsv> <output_dir> <picard> <genome_dir> <genome_fasta>
#
# samples.tsv columns: r1  r2  sample_id  genotype  barcode


# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <samples.tsv> <output_dir> <picard> <genome_dir> <genome_fasta>"
    exit 1
fi


# Assign arguments to variables
samples="$1"
output_dir="$2"
picard="$3"
genome_dir="$4"
genome="$5"

# Check if the input file exists
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Main loop
while IFS=$'\t' read -r f c a b d; do
  [[ "$f" =~ ^#  ]] && continue
  script_file="$a.filter.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.filter.out"
    echo "#SBATCH --error=${a}.filter.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH -p common,scavenger'
    echo "#SBATCH --chdir=$output_dir"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "source $(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "mkdir -p $output_dir/temp/$a"
    echo "$picard MarkDuplicates INPUT=$a.sort.bam OUTPUT=$a.MD.bam M=$a.metrics_file VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true"
    echo "samtools index $a.MD.bam" 
    echo "samtools view -h $a.MD.bam | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o $a.filtered.bam"
    echo "samtools index $a.filtered.bam"
    echo "$picard FixMateInformation INPUT=$a.filtered.bam OUTPUT=$a.FM.bam SORT_ORDER=coordinate TMP_DIR=$output_dir/temp/$a VALIDATION_STRINGENCY=LENIENT"
    echo "samtools index $a.FM.bam"
    echo "$picard AddOrReplaceReadGroups RGLB=$a RGPL=illumina RGPU=run RGSM=$b I=$a.FM.bam O=$a.RG.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT TMP_DIR=$output_dir/temp/$a"
    echo "samtools index $a.RG.bam"
  } > "$script_file"
done < "$samples"

