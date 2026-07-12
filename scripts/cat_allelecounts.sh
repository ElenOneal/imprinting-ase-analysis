#!/usr/bin/env bash
set -euo pipefail

# Concatenate per-chromosome allele count files into per-sample genome-wide files.
# Usage: bash 12_concat_counts.sh <samples.tsv> <count_dir> <output_suffix>
#
# samples.tsv columns: r1  r2  sample_id  job_id  cross

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 samples.tsv count_dir"
    exit 1
fi

samples="$1"
count_dir="$2"

if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

while IFS=$'\t' read -r r1 r2 sample_id job_id cross; do
    [[ "$r1" =~ ^# ]] && continue

    tmp_file="${count_dir}/${sample_id}.genecounts.tmp"
    out_file="${count_dir}/${sample_id}.allelecounts.txt"

    # Concatenate all per-chromosome count files for this sample
    cat "${count_dir}/${sample_id}".Chr_*_genecounts.txt > "$tmp_file"

    # Remove duplicate headers and write final output
    grep -v "^Gene" "$tmp_file" > "$out_file"

    rm "$tmp_file"
    echo "Written: $out_file"

done < "$samples"