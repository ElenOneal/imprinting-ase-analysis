#!/usr/bin/env bash
set -euo pipefail

# Filter parental bam files with a combination of Picard and samtools.
# Generates and submits a SLURM job script per sample.
#
# Usage: bash 03_filter_parent_bams.sh <samples.tsv> <output_dir> <picard> <genome_dir> <genome_fasta> <partition> <mem> [time]
#
# samples.tsv columns: r1  r2  sample_id  genotype  barcode
# Optional: Set CONDA_ENV environment variable to override conda environment name (default: imprinting-align)


# Check if the correct number of arguments is provided
if [ "$#" -lt 7 ] || [ "$#" -gt 8 ]; then
    echo "Usage: $0 <samples.tsv> <output_dir> <picard> <genome_dir> <genome_fasta> <partition> <mem> [time]"
    exit 1
fi

# Check if conda is available
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is not available in PATH."
    exit 2
fi

# Assign arguments to variables
samples="$1"
output_dir="$2"
picard="$3"
genome_dir="$4"
genome="$5"
partition="$6"
mem="$7"
time="${8:-01:00:00}"
java_mem=$(( ${mem%G} - 2 ))g
conda_env="${CONDA_ENV:-imprinting-align}"

# Validation
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

if [ ! -f "$picard" ]; then
    echo "Error: Picard jar '$picard' not found."
    exit 2
fi

if [ ! -d "$genome_dir" ]; then
    echo "Error: Genome directory '$genome_dir' not found."
    exit 2
fi

if [ ! -f "$genome" ]; then
    echo "Error: Genome fasta '$genome' not found."
    exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

script_count=0

# Main loop
while IFS=$'\t' read -r r1 r2 sample_id genotype barcode_; do
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue
  [[ -z "${sample_id:-}" ]] && continue
  script_file="$output_dir/$sample_id.filter.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --job-name=$sample_id"
    echo "#SBATCH --output=$output_dir/${sample_id}.filter.out"
    echo "#SBATCH --error=$output_dir/${sample_id}.filter.err"
    echo '#SBATCH --cpus-per-task=1'
    echo "#SBATCH -p $partition"
    echo "#SBATCH --chdir=$output_dir"
    echo "#SBATCH --mem=$mem"
    echo "#SBATCH --time=$time"
    echo ''
    echo "source \$(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate $conda_env"
    echo "mkdir -p $output_dir/temp/$sample_id"
    echo "java -Xmx${java_mem} -jar $picard MarkDuplicates INPUT=$sample_id.sort.bam OUTPUT=$sample_id.MD.bam M=$sample_id.metrics_file VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true"
    echo "samtools index $sample_id.MD.bam" 
    echo "samtools view -h $sample_id.MD.bam | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o $sample_id.filtered.bam"
    echo "samtools index $sample_id.filtered.bam"
    echo "java -Xmx${java_mem} -jar $picard FixMateInformation INPUT=$sample_id.filtered.bam OUTPUT=$sample_id.FM.bam SORT_ORDER=coordinate TMP_DIR=$output_dir/temp/$sample_id VALIDATION_STRINGENCY=LENIENT"
    echo "java -Xmx${java_mem} -jar $picard AddOrReplaceReadGroups RGLB=$sample_id RGPL=illumina RGPU=run RGSM=$genotype I=$sample_id.FM.bam O=$sample_id.RG.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT TMP_DIR=$output_dir/temp/$sample_id"
  } > "$script_file"
  chmod +x "$script_file"
  ((script_count++))
done < "$samples"

echo "Created $script_count filtering job scripts in $output_dir"
if [ "$script_count" -eq 0 ]; then
    echo "Warning: No scripts were created. Check that $samples has valid data."
    exit 1
fi
echo "To submit jobs, run:"
echo "for script in $output_dir/*.filter.sh; do sbatch \"\$script\"; done"