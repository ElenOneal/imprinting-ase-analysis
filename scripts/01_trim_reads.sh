#!/usr/bin/env bash
set -euo pipefail

# Trim reads with Trimmomatic. 
# Generates a separate SLURM job script for each sample, which can be submitted to the cluster for parallel processing. Each job will trim the paired-end reads for a single sample and save the cleaned reads to the specified output directory. The trimming parameters can be adjusted as needed, and the script assumes that the input file is a tab-delimited text file with the required columns.
# Usage: bash 01_trim_reads.sh <samples.tsv> <read_dir> <clean_dir> <execute_dir>
#
# samples.tsv columns: r1  r2  sample_id  genotype

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 samples.tsv read_dir clean_dir execute_dir"
    exit 1
fi

# Assign arguments to variables
samples="$1"
read_dir="$2"
clean_dir="$3"
execute_dir="$4"

# Check if the input file exists
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

# Create clean directory if it doesn't exist
mkdir -p "$clean_dir"

# Main loop
while IFS=$'\t' read -r f c a b; do
  [[ "$f" =~ ^#  ]] && continue
  script_file="$execute_dir/$a.trim.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.out"
    echo "#SBATCH --error=${a}.err"
    echo '#SBATCH --cpus-per-task=6'
    echo "#SBATCH --chdir=$execute_dir"
    echo '#SBATCH --mem=48G'
    echo ''
    echo "source $(conda info --base)/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "trimmomatic PE -threads 6 -phred33 -trimlog $clean_dir/$a.trimlog -quiet -validatePairs $read_dir/$f $read_dir/$c $clean_dir/$a.PE.R1.fq.gz $clean_dir/$a.U.R1.fq.gz $clean_dir/$a.PE.R2.fq.gz $clean_dir/$a.U.R2.fq.gz ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 AVGQUAL:25"
    echo "sbatch $execute_dir/${a}.align.sh"
  } > "$script_file"
done < "$samples"

