#!/usr/bin/env bash
set -euo pipefail

# Trim reads with Trimmomatic. 
# Generates a separate SLURM job script for each sample, which can be submitted to the cluster for parallel processing. Each job will trim the paired-end reads for a single sample and save the cleaned reads to the specified output directory. The trimming parameters can be adjusted as needed, and the script assumes that the input file is a tab-delimited text file with the required columns.
# Usage: bash 01_trim_reads.sh <samples.tsv> <read_dir> <clean_dir> <execute_dir> <partition>
#
# samples.tsv columns: r1  r2  sample_id  barcode population
# Check if the correct number of arguments is provided

if [ "$#" -ne 5 ]; then
echo "Usage: $0 samples.tsv read_dir clean_dir execute_dir partition"
exit 1
fi

# Assign arguments to variables
samples="$1"
read_dir="$2"
clean_dir="$3"
execute_dir="$4"
partition="$5"

# Check if the input file exists
if [ ! -f "$samples" ]; then
echo "Error: File '$samples' not found."
exit 2
fi

# Create clean directory if it doesn't exist
mkdir -p "$clean_dir" "$execute_dir"
# Main loop
while IFS=$'\t' read -r r1 r2 sample_id barcode _; do
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue
script_file="$execute_dir/${sample_id}.trim.sh"
  {
echo '#!/bin/bash'
echo '#'
echo "#SBATCH --get-user-env"
echo "#SBATCH --job-name=${sample_id}.trim"
echo "#SBATCH --output=${sample_id}.trim.out"
echo "#SBATCH --error=${sample_id}.trim.err"
echo '#SBATCH --cpus-per-task=6'
echo "#SBATCH --chdir=$execute_dir"
echo "#SBATCH -p $partition"
echo '#SBATCH --mem=48G'
echo ''
echo "source $(conda info --base)/etc/profile.d/conda.sh"
echo "conda activate imprinting-align"
echo "trimmomatic PE -threads 6 -phred33 -quiet -validatePairs $read_dir/$r1 $read_dir/$r2 $clean_dir/${sample_id}.PE.R1.fq.gz $clean_dir/${sample_id}.U.R1.fq.gz $clean_dir/${sample_id}.PE.R2.fq.gz $clean_dir/${sample_id}.U.R2.fq.gz ILLUMINACLIP:$HOME/TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:25 MINLEN:75 SLIDINGWINDOW:4:15 2>&1"
  } > "$script_file"
done < "$samples"
echo "All trimming scripts created. To submit jobs, run:"
echo "for script in *.trim.sh; do sbatch \"\$script\"; done"