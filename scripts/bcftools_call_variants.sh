#!/usr/bin/env bash
set -euo pipefail

# Call variants using bcftools as a SLURM array job.
# Usage: bash 04_call_variants.sh <listfile> <bamfiles> <output_dir> <genome_dir> <genome> <prefix> <parent1> <parent2> <gene_bed> <min_mq> <min_bq> <qual_thresh> <cluster_size> <cluster_window> <partition>

#ensure that the GATK genome dictionary is present in the folder with the genome
#if not, use gatk CreateSequenceDictionary -R genome.fa -O genome.dict to create it

if [ "$#" -ne 15 ]; then
    echo "Usage: $0 listfile bamfiles output_dir genome_dir genome prefix parent1 parent2 gene_bed min_mq min_bq qual_thresh cluster_size cluster_window partition"
    exit 1
fi

list="$1"
bamfiles="$2"
output_dir="$3"
genome_dir="$4"
genome="$5"
prefix="$6"
p1="$7"
p2="$8"
gene_bed="$9"
min_mq="${10}"
min_bq="${11}"
qual_thresh="${12}"
cluster_size="${13}"
cluster_window="${14}"
partition="${15}"

if [ ! -f "$list" ]; then
    echo "Error: File '$list' not found."
    exit 2
fi

if [ ! -f "$bamfiles" ]; then
    echo "Error: File '$bamfiles' not found."
    exit 2
fi

if [ ! -f "$gene_bed" ]; then
    echo "Error: File '$gene_bed' not found."
    exit 2
fi

# Validate list file format: should have exactly 2 columns, no headers
if ! awk 'NF != 2 {exit 1}' "$list"; then
    echo "Error: List file '$list' must have exactly 2 columns (region chrname) with no headers."
    exit 2
fi

# Count chromosomes for array size
nchr=$(grep -v '^#' "$list" | wc -l)

if [ "$nchr" -eq 0 ]; then
    echo "Error: No valid lines in list file '$list'."
    exit 2
fi

mkdir -p "$output_dir"

job_script="${output_dir}/bcftools_${prefix}.array.sh"

cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=bcftools_${prefix}
#SBATCH --output=${output_dir}/Chr_%a.${prefix}.out
#SBATCH --error=${output_dir}/Chr_%a.${prefix}.err
#SBATCH --cpus-per-task=6
#SBATCH -p $partition
#SBATCH --chdir=${output_dir}
#SBATCH --mem=30G
#SBATCH --array=1-${nchr}

source $(conda info --base)/etc/profile.d/conda.sh
conda activate imprinting-align

# Get the chromosome for this array task
region=\$(awk -v i=\$SLURM_ARRAY_TASK_ID 'NR==i {print \$1}' ${list})
chrname=\$(awk -v i=\$SLURM_ARRAY_TASK_ID 'NR==i {print \$2}' ${list})

bcftools mpileup --threads 6 --redo-BAQ --min-MQ ${min_mq} --min-BQ ${min_bq} \
    --per-sample-mF --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
    -f ${genome_dir}/${genome} -b ${bamfiles} -I -r \$region \
  | bcftools call --multiallelic-caller \
  | bcftools filter --threads 6 --SnpGap 3 \
    -e "QUAL<${qual_thresh} || INFO/RPBZ<-2 || INFO/RPBZ>2 || INFO/SCBZ<-2 || INFO/SCBZ>2" -Oz -o Chr_\${chrname}.${prefix}.filtered.vcf.gz

tabix Chr_\${chrname}.${prefix}.filtered.vcf.gz

bcftools view -s "$p1,$p2" -v snps -T "$gene_bed" Chr_\${chrname}.${prefix}.filtered.vcf.gz | bcftools filter -e 'GT=="mis" || GT=="het" || (GT[0]=="ref" && GT[1]=="ref")' -Oz -o Chr_\${chrname}.${prefix}.coding_snps.vcf.gz

tabix Chr_\${chrname}.${prefix}.coding_snps.vcf.gz

gatk VariantFiltration -R "$genome_dir/$genome" -V Chr_\${chrname}.${prefix}.coding_snps.vcf.gz -O Chr_\${chrname}.${prefix}.snpcluster.vcf --cluster-size "$cluster_size" --cluster-window-size "$cluster_window"

grep -v "SnpCluster" Chr_\${chrname}.${prefix}.snpcluster.vcf | bgzip > Chr_\${chrname}.${prefix}.final_snps.vcf.gz
tabix Chr_\${chrname}.${prefix}.final_snps.vcf.gz
EOF

chmod +x "$job_script"
echo "Created array job script: $job_script"
echo "Submit manually with: sbatch $job_script"