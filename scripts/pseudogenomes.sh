#!/usr/bin/env bash
set -euo pipefail

# Get filtered variant files for each parent, generate pseudogenomes with snps only. Then use STAR to create genomes for alignment of RNA-seq reads.
# Usage: bash 06_pseudogenomes.sh <vcffile> <output_dir> <genome_dir> <genome> <gtf> <parent1> <parent2> <gene_bed> <mindepth> <maxdepth>

if [ "$#" -ne 10 ]; then
    echo "Usage: $0 vcffile output_dir genome_dir genome gtf parent1 parent2 gene_bed mindepth maxdepth"
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

if [ ! -f "$vcffile" ]; then
    echo "Error: File '$vcffile' not found."
    exit 2
fi

# Generate and submit pseudogenome job
script_file="${output_dir}/pseudogenome.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=pgen"
echo "#SBATCH --output=${output_dir}/pgen.out"
echo "#SBATCH --error=${output_dir}/pgen.err"
echo '#SBATCH --cpus-per-task=1'
echo '#SBATCH -p common,scavenger'
echo "#SBATCH --chdir=${output_dir}"
echo '#SBATCH --mem=14G'
echo ''
echo "source /hpc/group/willislab/eo22/miniforge3/etc/profile.d/conda.sh"
echo "conda activate imprinting-align"
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
echo "cat ${output_dir}/${parent1}_pseudogenome.fa ${output_dir}/${parent2}_pseudogenome.fa > ${output_dir}/combined_pseudogenome.fa"
echo "awk -v prefix=\"${parent1}_\" '{ \$1 = prefix \$1; print }' OFS='\t' \"${genome_dir}/${gtf}\" > \"${output_dir}/${parent1}.gtf\""
echo "awk -v prefix=\"${parent2}_\" '{ \$1 = prefix \$1; print }' OFS='\t' \"${genome_dir}/${gtf}\" > \"${output_dir}/${parent2}.gtf\""
echo "cat ${output_dir}/$parent1.gtf ${output_dir}/$parent2.gtf > ${output_dir}/combined.gtf"
echo "mkdir -p ${output_dir}/star_index"
echo "cp ${output_dir}/combined_pseudogenome.fa ${output_dir}/star_index/combined_pseudogenome.fa"
echo "cp ${output_dir}/combined.gtf ${output_dir}/star_index/combined.gtf"
echo "sbatch ${output_dir}/star_index.sh"
} > "$script_file"

# Create STAR index script
star_script_file="${output_dir}/star_index.sh"
{
echo '#!/bin/bash'
echo "#SBATCH --job-name=star_index"
echo "#SBATCH --output=${output_dir}/star_index.out"
echo "#SBATCH --error=${output_dir}/star_index.err"
echo '#SBATCH --cpus-per-task=6'
echo '#SBATCH -p common,scavenger'
echo "#SBATCH --chdir=${output_dir}/star_index"
echo '#SBATCH --mem=30G'
echo ''
echo "source /hpc/group/willislab/eo22/miniforge3/etc/profile.d/conda.sh"
echo "conda activate imprinting-align"
echo ''
echo "STAR --runThreadN 6 --runMode genomeGenerate \\"
echo "     --genomeDir . \\"
echo "     --genomeFastaFiles combined_pseudogenome.fa \\"
echo "     --sjdbGTFfile combined.gtf \\"
echo "     --genomeSAindexNbases 13"
} > "$star_script_file"

echo "Submitting ${script_file}..."
#sbatch "$script_file"