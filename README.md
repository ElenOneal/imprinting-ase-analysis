# imprinting-ase-analysis

Scripts for detecting imprinted genes in the *Mimulus guttatus* species complex using allele-specific expression (ASE) in hybrid endosperm. 

## Table of Contents
- [Project Overview](#project-overview)
- [Pipeline Overview](#pipeline-overview)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Data Access](#data-access)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Results](#results)

## Project Overview

Genomic imprinting is a form of epigenetic regulation in which gene expression is biased toward one parental allele. It has evolved independently in mammals and flowering plants, and is theorized to result from conflict between maternal and paternal genomes over resource allocation to offspring. In plants, imprinting is largely restricted to the endosperm — a triploid nutritive tissue formed by fusion of the diploid maternal central cell with a haploid paternal gamete. Under normal biallelic expression, endosperm genes conform to a 2:1 maternal:paternal ratio. Imprinted genes deviate from this expectation: maternally expressed genes (MEGs) are biased toward the maternal allele, and paternally expressed genes (PEGs) toward the paternal allele.

This pipeline identifies MEGs and PEGs from reciprocal hybrid crosses within the *M. guttatus* species complex, using parent-of-origin SNPs to quantify allele-specific expression in hybrid endosperm transcriptomes.

**Note:** This pipeline was developed using *M. guttatus* IM767 v2 reference genome.  A newer genome assembly is available and updating the pipeline to include multiple genomes planned. The reference annotation file (`Mguttatusvar_IM767_887_v2.1.gene_exons.gff3`) 
can be downloaded from Phytozome: https://phytozome-next.jgi.doe.gov

## Pipeline Overview

**Stage 1 — SNP discovery (parental DNA-seq)**
- Align parental reads to reference genome with BWA
- Call variants with bcftools to identify parent-of-origin SNPs

**Note:** `call_pairwise_snps_w_genes.py` currently assumes Phytozome GFF3  annotation format. Users with other annotation formats may need to adjust the attribute parsing in `make_gene_dict` and `make_cds_dict`.

**Stage 2 — ASE quantification (hybrid endosperm RNA-seq)**
- Align hybrid endosperm reads with STAR
- Quantify allele-specific expression at parental SNPs using custom Python scripts

**Stage 3 — Imprinting classification (R)**
- Test deviation from expected 2:1 maternal:paternal expression ratio
- Classify genes as MEGs, PEGs, or biallelic
- Generate summary tables and figures

## Status 
This repository is under active development. 

The following scripts are implemented and tested:
- `01_trim_reads.sh` — Trimmomatic quality trimming
- `02_align_parents.sh` — BWA alignment of parental DNA-seq
- `03_filter_parent_bams.sh` — Filter parental BAM files
- `04_call_variants.sh` — bcftools SNP calling
- `05_snp_catalog.sh` — Extract coding SNPs per gene
- `06_index_star.sh` — Build STAR index
- `07_align_star.sh` — STAR alignment of hybrid endosperm RNA-seq
- `08_reindex_star.sh` — Re-index BAM files
- `09_realign_star.sh` — Re-align to reference
- `10_filter_rna.sh` — Filter RNA BAM files for high-quality, uniquely mapped reads
- `11_snpcount_array.sh` — Count parental alleles per gene (SLURM array job)
- `call_pairwise_snps_w_genes.py` — Extract parent-of-origin SNPs within genes
- `classify_parental_alleles.py` — Quantify allele-specific read counts per gene
- `12_cat_allelecounts.sh` — Concatenate per-chromosome allele count files into per-sample genome-wide files

## Repository Structure

```
imprinting-ase-analysis/
├── scripts/
│   ├── 01_trim_reads.sh
│   ├── 02_align_parents.sh
│   ├── 03_filter_parent_bams.sh
│   ├── 04_call_variants.sh
│   ├── 05_snp_catalog.sh
│   ├── 06_index_star.sh
│   ├── 07_align_star.sh
│   ├── 08_reindex_star.sh
│   ├── 09_realign_star.sh
│   ├── 10_filter_rna.sh
│   ├── 11_snpcount_array.sh
│   ├── 12_cat_allelecounts.sh
│   ├── R/
│   │   ├── edger_preanalysis.R
│   │   ├── filter_snpcounts.R
│   │   ├── get_megs_pegs.R
│   │   ├── imprinting_analysis.R
│   │   ├── merge_snpcounts.R
│   │   ├── plot_imprinted_genes.R
│   │   ├── plot_mds.R
│   │   └── snpcounts_by_parent.R
│   ├── python/
│   │   ├── call_pairwise_snps_w_genes.py
│   │   └── classify_parental_alleles.py
├── config/
│   ├── chromosomes.list         # Chromosome names
│   ├── im767_v2.genes.bed       # Gene coordinates
│   ├── parental_samples.tsv     # Parental sample metadata
│   ├── samples.tsv              # Hybrid endosperm sample metadata
│   └── TruSeq.fa                # Trimmomatic adapter sequences
├── envs/
│   └── imprinting-align.yml     # Conda environment
├── LICENSE
└── README.md
```

##  Data Access

The full dataset and figures will be available upon publication.  For now, this repository provides code and illustrative examples only.

## Getting Started

### Prerequisites

- [conda](https://docs.conda.io/en/latest/)
- Git

### Installation
```bash
git clone git@github.com:ElenOneal/imprinting-ase-analysis.git
cd imprinting-ase-analysis
conda env create -f envs/imprinting-align.yml
conda activate imprinting-align
```
## Usage

Each script is standalone and processes one stage of the pipeline. Scripts are numbered and should be run in order. Run individual scripts with appropriate arguments:

```bash
# Stage 1: Parental SNP discovery (DNA-seq)
bash scripts/01_trim_reads.sh <input_dir> <output_dir>
bash scripts/02_align_parents.sh <samples.tsv> <output_dir>
bash scripts/03_filter_parent_bams.sh <samples.tsv> <bam_dir> <output_dir>
bash scripts/04_call_variants.sh <samples.tsv> <output_dir>
bash scripts/05_snp_catalog.sh <vcf_dir> <gff3> <output_dir>

# Stage 2: Hybrid endosperm RNA-seq alignment and ASE quantification
bash scripts/06_index_star.sh <genome> <annotation> <output_dir>
bash scripts/07_align_star.sh <samples.tsv> <index_dir> <output_dir>
bash scripts/08_reindex_star.sh <samples.tsv> <bam_dir>
bash scripts/09_realign_star.sh <samples.tsv> <index_dir> <output_dir>
bash scripts/10_filter_rna.sh <samples.tsv> <output_dir>
bash scripts/11_snpcount_array.sh <samples.tsv> <snp_catalog> <bam_dir> <output_dir> <scripts_dir> <partition> <conda_env>
bash scripts/12_cat_allelecounts.sh <samples.tsv> <output_dir>
```

**Outputs:** 
- Stage 1: parent-of-origin SNP catalog (.bed format)
- Stage 2: gene-level parental allele counts per sample (_genecounts.txt)
- Final: HTML summary report ([NudatusImprinting.html](NudatusImprinting.html))

## Dependencies

Managed via conda (`envs/imprinting-align.yml`). Key tools:

| Tool | Version | Purpose |
|------|---------|---------|
| BWA | 0.7.19 | Parental DNA-seq alignment |
| bcftools | 1.23.1 | SNP calling |
| STAR | 2.7.11b | Endosperm RNA-seq alignment |
| Python | 3.8+ | ASE quantification |
| R | 4.2.3 | Statistical analysis and figures |

## License

MIT License. See [LICENSE](LICENSE) for details.

## Results

A rendered HTML report of the imprinting analysis is available at
[NudatusImprinting.html](NudatusImprinting.html).

## Acknowledgements
Willis Lab (Duke) / NSF / Dr. John Willis
Franks Lab (NCSU) / NSF / Miguel Flores, Dr. Robert Franks
