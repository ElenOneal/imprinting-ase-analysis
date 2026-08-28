#!/usr/bin/env bash
set -euo pipefail

# Filter parental bam files with a combination of Picard and samtools.
# Generates and submits a SLURM job script per BIOLOGICAL SAMPLE
# (grouped by the "genotype" column), not per fastq/lane row.
#
# Usage: bash filter_parent_bams.sh <samples.tsv> <output_dir> <partition> <mem> [time]
#
# samples.tsv columns: r1  r2  sample_id  genotype  barcode  parent  species  lane
#
# - sample_id  : unique per-lane identifier. This is the basename of the
#                ALREADY-ALIGNED, coordinate-sorted input bam, i.e. this
#                script expects "${sample_id}.sort.bam" to exist.
# - genotype   : the true biological sample identity. Rows that share the
#                same genotype are treated as multiple lanes of the SAME
#                sample and are merged together (see below). This becomes
#                the basename for all downstream output files, and the
#                RGSM/RGLB tag.
# - lane       : lane number for this row (used to build RGPU so multi-lane
#                reads are tagged with which physical lane they came from).
#
# Why grouping is done on genotype and not sample_id: sequencers/demux
# pipelines often assign a *different* S-number (and therefore a different
# sample_id, e.g. CSH10_S6 vs CSH10_S21) to the same physical sample when
# it's split across lanes. sample_id is only guaranteed unique per fastq
# pair/lane, not per biological sample -- genotype is the column that
# actually identifies "this is the same sample".
#
# Multi-lane handling:
#   For a genotype with N>1 rows (lanes), each lane's sort.bam is first
#   tagged individually with AddOrReplaceReadGroups (distinct --RGID per
#   lane, RGPU encoding the lane number, but a shared RGSM/RGLB = genotype).
#   All of that genotype's tagged, per-lane bams are then passed to a
#   SINGLE MarkDuplicates call as multiple --INPUT arguments (Picard
#   accepts repeated --INPUT and treats them as one combined set). This is
#   deliberate: PCR duplicates from the same library can appear on more
#   than one lane, and MarkDuplicates can only catch duplicates that are
#   compared against each other in the same run -- it doesn't matter that
#   the reads share RGSM/RGLB tags if they're never fed to the same
#   invocation. Passing multiple --INPUT files achieves this without the
#   extra step (and extra intermediate bam) of physically merging with
#   samtools merge first. Optical duplicate detection still correctly
#   scopes to within a single lane/tile, since Picard's READ_NAME_REGEX
#   parses lane/tile/x/y straight from the read name regardless of how
#   many files were given.
#   For a genotype with a single row, MarkDuplicates just gets one --INPUT.


# Picard is expected to be installed in the f1mapping conda environment
# (see environment.yml) and invoked as the `picard` command rather than
# `java -jar <path>`.


# Check if the correct number of arguments is provided
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
  echo "Usage: $0 <samples.tsv> <output_dir> <partition> <mem> [time]"
    exit 1
fi


# Assign arguments to variables
samples="$1"
output_dir="$2"
partition="$3"
mem="$4"
time_limit="${5:-24:00:00}"

# Require an explicit G suffix (e.g. 48G). SLURM interprets --mem values
# with no unit suffix as megabytes, so a bare "48" silently becomes 48MB
# instead of 48GB -- the job then gets OOM-killed almost immediately with
# no obvious clue why, since the java_mem calculation below still "works"
# whether or not G is present.
if ! [[ "$mem" =~ ^[0-9]+G$ ]]; then
    echo "Error: --mem must be a whole number followed by G (e.g. 48G). Got: '$mem'"
    exit 4
fi

# Leave a 6G buffer between the JVM heap ceiling and the SLURM --mem hard
# limit. Picard's off-heap memory use (native codecs, metaspace, thread
# stacks) sits outside -Xmx, so a small gap (e.g. 2G) can still let a job
# get OOM-killed even at large --mem.
java_mem=$(( ${mem%G} - 6 ))g

if [ "${mem%G}" -le 6 ]; then
    echo "Error: --mem ($mem) must be greater than 6G to leave room for JVM overhead."
    exit 3
fi

# Check if the input file exists
if [ ! -f "$samples" ]; then
    echo "Error: File '$samples' not found."
    exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Compute once, outside the loop -- conda info --base is slow (~1-5s) and
# calling it once per sample inside the loop adds up for large sample sheets.
conda_base="$(conda info --base)"

# ---------------------------------------------------------------------------
# Pass 1: read the tsv and group rows by genotype, preserving first-seen
# order of genotypes. group_rows[$genotype] accumulates one line per row
# as "sample_id<TAB>lane".
# ---------------------------------------------------------------------------
declare -A group_rows
declare -A seen
genotype_order=()

while IFS=$'\t' read -r r1 r2 sample_id genotype barcode parent species lane; do
  [[ "$r1" =~ ^# ]] && continue
  [[ -z "${r1:-}" ]] && continue
  if [[ -z "${seen[$genotype]:-}" ]]; then
    seen[$genotype]=1
    genotype_order+=("$genotype")
  fi
  group_rows["$genotype"]+="${sample_id}"$'\t'"${lane}"$'\n'
done < "$samples"

# ---------------------------------------------------------------------------
# Pass 2: one SLURM script per genotype (biological sample).
# ---------------------------------------------------------------------------
for genotype in "${genotype_order[@]}"; do
  # Split this genotype's accumulated rows into an array, one entry per lane.
  mapfile -t lane_rows <<< "${group_rows[$genotype]}"
  # mapfile picks up a trailing empty element from the final \n -- drop it.
  if [[ -z "${lane_rows[-1]:-}" ]]; then
    unset 'lane_rows[-1]'
  fi
  n_lanes="${#lane_rows[@]}"

  script_file="${genotype}.filter.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --job-name=${genotype}"
    echo "#SBATCH --output=${genotype}.filter.out"
    echo "#SBATCH --error=${genotype}.filter.err"
    echo '#SBATCH --cpus-per-task=1'
    echo "#SBATCH -p $partition"
    echo "#SBATCH --time=$time_limit"
    echo "#SBATCH --chdir=$output_dir"
    echo "#SBATCH --mem=$mem"
    echo ''
    echo "source ${conda_base}/etc/profile.d/conda.sh"
    echo "conda activate imprinting-align"
    echo "mkdir -p $output_dir/temp/${genotype}"

    md_inputs=()
    for row in "${lane_rows[@]}"; do
      IFS=$'\t' read -r sample_id lane <<< "$row"
      rgpu="flowcell.lane${lane}"
      tagged_bam="${sample_id}.tagged.bam"
      md_inputs+=(--INPUT "$tagged_bam")
      # --RGID must be unique per lane so the combined header ends up with
      # distinct @RG lines instead of Picard's collapsing them all under
      # the default RGID of "1".
      echo "picard -Xmx${java_mem} AddOrReplaceReadGroups --RGID ${sample_id} --RGLB ${genotype} --RGPL illumina --RGPU ${rgpu} --RGSM ${genotype} --INPUT ${sample_id}.sort.bam --OUTPUT ${tagged_bam} --SORT_ORDER coordinate --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --TMP_DIR $output_dir/temp/${genotype}"
    done

    # One --INPUT per lane (Picard accepts repeated --INPUT and dedups
    # across all of them jointly); no separate merge step needed.
    echo "picard -Xmx${java_mem} MarkDuplicates ${md_inputs[*]} --OUTPUT ${genotype}.MD.bam --METRICS_FILE ${genotype}.metrics_file --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true"
    echo "samtools index ${genotype}.MD.bam"
    echo "samtools view -h ${genotype}.MD.bam | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o ${genotype}.filtered.bam"
    echo "samtools index ${genotype}.filtered.bam"
    echo "picard -Xmx${java_mem} FixMateInformation --INPUT ${genotype}.filtered.bam --OUTPUT ${genotype}.RG.bam --SORT_ORDER coordinate --TMP_DIR $output_dir/temp/${genotype} --VALIDATION_STRINGENCY LENIENT"
    echo "samtools index ${genotype}.RG.bam"
  } > "$script_file"
  chmod +x "$script_file"
done

echo "All filtering scripts created. To submit jobs, run:"
echo "for script in *.filter.sh; do sbatch \"\$script\"; done"
