
#!/usr/bin/env bash
set -euo pipefail

# Create STAR index
# Usage: bash 07_index_star.sh <genome_dir> <genome> <gtf> <partition>

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 genome_dir genome gtf partition"
    exit 1
fi

# Assign arguments to variables
genome_dir="$1"
genome="$2"
gtf="$3"
partition="$4"

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


echo '#!/bin/bash' > index.sh
echo '#SBATCH --job-name=star_index' >> index.sh
echo '#SBATCH --output=star_index.out' >> index.sh
echo '#SBATCH --error=star_index.err' >> index.sh
echo "#SBATCH -p $partition" >> index.sh
echo '#SBATCH --cpus-per-task=4' >> index.sh
echo "#SBATCH --chdir=$genome_dir" >> index.sh
echo '#SBATCH --mem=24G' >> index.sh
echo '' >> index.sh
echo "source $(conda info --base)/etc/profile.d/conda.sh" >> index.sh
echo "conda activate imprinting-align" >> index.sh
echo "STAR --runMode genomeGenerate \\" >> index.sh
echo "     --genomeDir . \\" >> index.sh
echo "     --genomeFastaFiles $genome \\" >> index.sh
echo "     --runThreadN 4 \\" >> index.sh
echo "     --genomeSAindexNbases 13 \\" >> index.sh
echo "     --sjdbGTFfile $gtf \\" >> index.sh
echo "     --sjdbOverhang 149 \\" >> index.sh
echo "     --limitSjdbInsertNsj 10000000" >> index.sh

chmod 755 index.sh

echo "STAR index script created: index.sh"
echo "To submit the job, run: sbatch index.sh"