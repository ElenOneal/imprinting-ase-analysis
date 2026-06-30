#!/usr/bin/env bash
set -euo pipefail

# Reindex STAR genome index with junctions detected by first round of alignment
# Usage: bash 08_reindex_star.sh <genome_dir> <genome> <gtf> <tab_dir>

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 genome_dir genome gtf tab_dir"
    exit 1
fi

# Assign arguments to variables
genome_dir="$1"
genome="$2"
gtf="$3"
tab_dir="$4"

# Validate inputs
if [ ! -d "$genome_dir" ]; then
    echo "Error: Genome directory '$genome_dir' does not exist."
    exit 1
fi

if [ ! -f "$genome_dir/$genome" ]; then
    echo "Error: Genome file '$genome_dir/$genome' does not exist."
    exit 1
fi

if [ ! -f "$genome_dir/$gtf" ]; then
    echo "Error: GTF file '$genome_dir/$gtf' does not exist."
    exit 1
fi

if [ ! -d "$tab_dir" ]; then
    echo "Error: Tab directory '$tab_dir' does not exist."
    exit 1
fi

reindex_script="$genome_dir/reindex.sh"

{
echo '#!/bin/bash'
echo '#SBATCH --job-name=star_reindex'
echo '#SBATCH --output=star_reindex.out'
echo '#SBATCH --error=star_reindex.err'
echo '#SBATCH -p common,scavenger'
echo '#SBATCH --cpus-per-task=6'
echo "#SBATCH --chdir=$genome_dir"
echo '#SBATCH --mem=24G'
echo ''
echo "source $(conda info --base)/etc/profile.d/conda.sh"
echo "conda activate imprinting-align"
echo "# Copy splice junction files"
echo "cp $tab_dir/*.tab ."
echo ''
echo "# Filter and prepare splice junctions"
echo "cat *.tab | awk '(\$5 > 0 && \$7 > 2 && \$6==0)' | cut -f1-6 | sort | uniq > spliced.tab"
echo ''
echo "# Re-index genome with detected splice junctions"
echo "STAR --runMode genomeGenerate \\"
echo "     --genomeDir . \\"
echo "     --genomeFastaFiles $genome \\"
echo "     --runThreadN 6 \\"
echo "     --genomeSAindexNbases 13 \\"
echo "     --sjdbGTFfile $gtf \\"
echo "     --sjdbOverhang 149 \\"
echo "     --sjdbFileChrStartEnd spliced.tab \\"
echo "     --limitSjdbInsertNsj 10000000"
} > "$reindex_script"

chmod 755 "$reindex_script"

echo "STAR re-indexing script created: $reindex_script"
echo "To submit the job, run: sbatch $reindex_script"