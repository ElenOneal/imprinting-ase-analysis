#!/usr/bin/env bash
set -euo pipefail

# Get filtered variant files for each parent, generate pseudogenomes with snps only. Then use STAR to create genomes for alignment of RNA-seq reads.
# Usage: bash scripts/generate_pseudogenome.sh <vcffile> <output_dir> <genome_dir> <genome> <gtf> <parent1> <parent2> <gene_bed> <mindepth> <maxdepth> <partition> <mem> [time]

if [ "$#" -lt 12 ] || [ "$#" -gt 13 ]; then
    echo "Usage: $0 vcffile output_dir genome_dir genome gtf parent1 parent2 gene_bed mindepth maxdepth partition mem [time]"
    exit 1
fi

vcffile="$1"
output_dir="$2"
genome_dir="$3"
genome="$4"
gtf="$5"
parent1="$6"
parent2="$7"
gene_bed="$8"
mindepth="$9"
maxdepth="${10}"
partition="${11}"
mem="${12}"
time="${13:-01:00:00}"

if [ ! -f "$vcffile" ]; then
    echo "Error: File '$vcffile' not found."
    exit 2
fi

if [ ! -f "$gene_bed" ]; then
    echo "Error: File '$gene_bed' not found."
    exit 2
fi

if [ ! -f "${genome_dir}/${genome}" ]; then
    echo "Error: File '${genome_dir}/${genome}' not found."
    exit 2
fi

if [ ! -f "${genome_dir}/${gtf}" ]; then
    echo "Error: File '${genome_dir}/${gtf}' not found."
    exit 2
fi

if ! [[ "$mindepth" =~ ^[0-9]+$ ]]; then
    echo "Error: mindepth must be an integer, got '$mindepth'."
    exit 2
fi

if ! [[ "$maxdepth" =~ ^[0-9]+$ ]]; then
    echo "Error: maxdepth must be an integer, got '$maxdepth'."
    exit 2
fi

if (( mindepth > maxdepth )); then
    echo "Error: mindepth ($mindepth) cannot be greater than maxdepth ($maxdepth)."
    exit 2
fi

if [ -z "$partition" ]; then
    echo "Error: partition must be a non-empty string."
    exit 2
fi

if ! [[ "$mem" =~ ^[0-9]+[KMGTP]$ ]]; then
    echo "Error: mem must look like 14G, 30000M, etc. Got '$mem'."
    exit 2
fi

if ! [[ "$time" =~ ^[0-9]+:[0-9]{2}:[0-9]{2}$ ]]; then
    echo "Error: time must be in HH:MM:SS format. Got '$time'."
    exit 2
fi

mkdir -p "$output_dir"

# Resolve conda setup; CONDA_SH can override auto-detection.
if [ -n "${CONDA_SH:-}" ]; then
    conda_sh="$CONDA_SH"
else
    if ! command -v conda >/dev/null 2>&1; then
        echo "Error: conda was not found in PATH and CONDA_SH is not set."
        echo "Set CONDA_SH to your conda.sh path, then rerun."
        exit 2
    fi
    conda_base="$(conda info --base 2>/dev/null || true)"
    if [ -z "$conda_base" ]; then
        echo "Error: unable to determine conda base path from 'conda info --base'."
        echo "Set CONDA_SH to your conda.sh path, then rerun."
        exit 2
    fi
    conda_sh="${conda_base}/etc/profile.d/conda.sh"
fi

conda_env="${CONDA_ENV:-imprinting-align}"

if [ ! -f "$conda_sh" ]; then
    echo "Error: conda init script '$conda_sh' not found."
    echo "Set CONDA_SH to the correct path, then rerun."
    exit 2
fi

for tool in sbatch awk sed cp cat mkdir chmod; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: required command '$tool' not found in PATH."
        exit 2
    fi
done

# Generate pseudogenome and STAR index job scripts
script_file="${output_dir}/pseudogenome.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=pgen"
echo "#SBATCH --output=${output_dir}/pgen.out"
echo "#SBATCH --error=${output_dir}/pgen.err"
echo '#SBATCH --cpus-per-task=1'
echo "#SBATCH -p ${partition}"
echo "#SBATCH --chdir=${output_dir}"
echo "#SBATCH --mem=$mem"
echo "#SBATCH --time=$time"
echo ''
echo 'set -euo pipefail'
echo ''
echo "source ${conda_sh}"
echo "conda activate ${conda_env}"
echo "command -v bcftools >/dev/null 2>&1 || { echo 'Error: bcftools not found in PATH.'; exit 2; }"
echo "command -v tabix >/dev/null 2>&1 || { echo 'Error: tabix not found in PATH.'; exit 2; }"
echo "command -v sbatch >/dev/null 2>&1 || { echo 'Error: sbatch not found in PATH.'; exit 2; }"
echo ''
echo "bcftools view -s \"$parent1\" -v snps -T \"$gene_bed\" \"$vcffile\" | bcftools filter -e 'GT==\"mis\" || GT==\"het\" || FORMAT/DP<${mindepth} || FORMAT/DP>${maxdepth}' -Oz -o \"${output_dir}/${parent1}.genic_snps.vcf.gz\""
echo "tabix \"${output_dir}/${parent1}.genic_snps.vcf.gz\""
echo "bcftools consensus -f \"${genome_dir}/${genome}\" -s \"$parent1\" -o \"${output_dir}/${parent1}_pseudogenome.fa\" \"${output_dir}/${parent1}.genic_snps.vcf.gz\""
echo "sed -i 's/>/>${parent1}_/g' \"${output_dir}/${parent1}_pseudogenome.fa\""
echo ''
echo "bcftools view -s \"$parent2\" -v snps -T \"$gene_bed\" \"$vcffile\" | bcftools filter -e 'GT==\"mis\" || GT==\"het\" || FORMAT/DP<${mindepth} || FORMAT/DP>${maxdepth}' -Oz -o \"${output_dir}/${parent2}.genic_snps.vcf.gz\""
echo "tabix \"${output_dir}/${parent2}.genic_snps.vcf.gz\""
echo "bcftools consensus -f \"${genome_dir}/${genome}\" -s \"$parent2\" -o \"${output_dir}/${parent2}_pseudogenome.fa\" \"${output_dir}/${parent2}.genic_snps.vcf.gz\""
echo "sed -i 's/>/>${parent2}_/g' \"${output_dir}/${parent2}_pseudogenome.fa\""
echo "cat \"${output_dir}/${parent1}_pseudogenome.fa\" \"${output_dir}/${parent2}_pseudogenome.fa\" > \"${output_dir}/combined_pseudogenome.fa\""
echo "awk -v prefix=\"${parent1}_\" '{ \$1 = prefix \$1; print }' OFS='\t' \"${genome_dir}/${gtf}\" > \"${output_dir}/${parent1}.gtf\""
echo "awk -v prefix=\"${parent2}_\" '{ \$1 = prefix \$1; print }' OFS='\t' \"${genome_dir}/${gtf}\" > \"${output_dir}/${parent2}.gtf\""
echo "cat \"${output_dir}/${parent1}.gtf\" \"${output_dir}/${parent2}.gtf\" > \"${output_dir}/combined.gtf\""
echo "mkdir -p \"${output_dir}/star_index\""
echo "cp \"${output_dir}/combined_pseudogenome.fa\" \"${output_dir}/star_index/combined_pseudogenome.fa\""
echo "cp \"${output_dir}/combined.gtf\" \"${output_dir}/star_index/combined.gtf\""
echo "sbatch \"${output_dir}/star_index.sh\""
} > "$script_file"
chmod +x "$script_file"

# Create STAR index script
star_script_file="${output_dir}/star_index.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=star_index"
echo "#SBATCH --output=${output_dir}/star_index.out"
echo "#SBATCH --error=${output_dir}/star_index.err"
echo '#SBATCH --cpus-per-task=6'
echo "#SBATCH -p ${partition}"
echo "#SBATCH --chdir=${output_dir}/star_index"
echo "#SBATCH --mem=$mem"
echo "#SBATCH --time=$time"
echo ''
echo 'set -euo pipefail'
echo ''
echo "source ${conda_sh}"
echo "conda activate ${conda_env}"
echo "command -v STAR >/dev/null 2>&1 || { echo 'Error: STAR not found in PATH.'; exit 2; }"
echo ''
echo "STAR --runThreadN 6 --runMode genomeGenerate \\"
echo "     --genomeDir . \\"
echo "     --genomeFastaFiles combined_pseudogenome.fa \\"
echo "     --sjdbGTFfile combined.gtf \\"
echo "     --genomeSAindexNbases 13"
} > "$star_script_file"
chmod +x "$star_script_file"

echo "Pseudogenome script created. To submit jobs, run:"
echo "sbatch ${script_file}"