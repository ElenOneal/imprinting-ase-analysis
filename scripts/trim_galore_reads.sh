#!/usr/bin/env bash
set -euo pipefail

# Usage: bash 01_trim_galore_reads.sh <samples.tsv> <read_dir> <clean_dir> <execute_dir> <partition>
# samples.tsv columns: r1  r2  sample_id  barcode  lane  replicate  cross_direction  species  notes
#
# One job is generated per LANE (per row). Samples with multiple lanes get multiple
# job scripts, named ${sample_id}.${lane}.trim_rna.sh so they don't collide.
#
# Optional: Set TRIM_GALORE_SIF environment variable to override the container path:
#   export TRIM_GALORE_SIF=/path/to/trim-galore.sif

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 samples.tsv read_dir clean_dir execute_dir partition"
    exit 1
fi

samples="$1"
read_dir="$2"
clean_dir="$3"
execute_dir="$4"
partition="$5"

# Path to Apptainer image (override with TRIM_GALORE_SIF environment variable)
trim_galore_sif="${TRIM_GALORE_SIF:-$HOME/containers/trim-galore.sif}"

if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

for dir in "$read_dir" "$clean_dir" "$execute_dir"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' not found."
        exit 2
    fi
done

if [ ! -f "$trim_galore_sif" ]; then
    echo "Error: Apptainer image '$trim_galore_sif' not found."
    exit 2
fi

mkdir -p "$clean_dir" "$execute_dir"

script_count=0

while IFS=$'\t' read -r r1 r2 sample_id barcode lane _; do
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue
  [[ -z "${lane:-}" ]] && continue

  lane_id="${sample_id}.${lane}"
  script_file="$execute_dir/${lane_id}.trim_rna.sh"

  {
    echo '#!/bin/bash'
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=${lane_id}.trim"
    echo "#SBATCH --output=\"$execute_dir/${lane_id}.trim.out\""
    echo "#SBATCH --error=\"$execute_dir/${lane_id}.trim.err\""
    echo "#SBATCH --cpus-per-task=6"
    echo "#SBATCH --chdir=\"$execute_dir\""
    echo "#SBATCH -p $partition"
    echo "#SBATCH --mem=24G"
    echo ''
    echo 'set -euo pipefail'
    echo ''
    echo "apptainer exec --bind \"$read_dir\":/reads,\"$clean_dir\":/clean \"$trim_galore_sif\" \\"
    echo "  trim_galore \\"
    echo "    --paired \\"
    echo "    --cores 6 \\"
    echo "    --quality 25 \\"
    echo "    --length 36 \\"
    echo "    --fastqc \\"
    echo "    --basename ${lane_id} \\"
    echo "    --output_dir /clean \\"
    echo "    /reads/$r1 /reads/$r2"
  } > "$script_file"

  chmod +x "$script_file"
  ((script_count++))

done < "$samples"

echo "Created $script_count Trim Galore job scripts in $execute_dir"
if [ "$script_count" -eq 0 ]; then
    echo "Warning: No scripts were created. Check that $samples has valid data."
    exit 1
fi
echo "To submit jobs, run:"
echo "for script in $execute_dir/*.trim_rna.sh; do sbatch \"\$script\"; done"
