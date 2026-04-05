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

## Project Overview

Genomic imprinting is a form of epigenetic regulation in which gene expression is biased toward one parental allele. It has evolved independently in mammals and flowering plants, and is theorized to result from conflict between maternal and paternal genomes over resource allocation to offspring. In plants, imprinting is largely restricted to the endosperm — a triploid nutritive tissue formed by fusion of the diploid maternal central cell with a haploid paternal gamete. Under normal biallelic expression, endosperm genes conform to a 2:1 maternal:paternal ratio. Imprinted genes deviate from this expectation: maternally expressed genes (MEGs) are biased toward the maternal allele, and paternally expressed genes (PEGs) toward the paternal allele.

This pipeline identifies MEGs and PEGs from reciprocal hybrid crosses within the *M. guttatus* species complex, using parent-of-origin SNPs to quantify allele-specific expression in hybrid endosperm transcriptomes.

**Note:** This pipeline was developed using *M. guttatus* IM767 v2 reference genome. 
A newer genome assembly is available and updating the pipeline to include multiple genomes planned.

## Pipeline Overview

**Stage 1 — SNP discovery (parental DNA-seq)**
- Align parental reads to reference genome with BWA
- Call variants with bcftools to identify parent-of-origin SNPs

**Note:** `call_pairwise_snps_w_genes.py` currently assumes Phytozome GFF3 
annotation format. Users with other annotation formats may need to adjust 
the attribute parsing in `make_gene_dict` and `make_cds_dict`.

**Stage 2 — ASE quantification (hybrid endosperm RNA-seq)**
- Align hybrid endosperm reads with STAR
- Quantify allele-specific expression at parental SNPs using custom Python scripts

**Stage 3 — Imprinting classification (R)**
- Test deviation from expected 2:1 maternal:paternal expression ratio
- Classify genes as MEGs, PEGs, or biallelic
- Generate summary tables and figures

## Repository Structure

imprinting-ase-analysis/
├── scripts/
│   ├── 01_align_parental.sh       # BWA alignment of parental DNA-seq
│   ├── 02_call_snps.sh            # bcftools SNP calling
│   ├── 03_align_endosperm.sh      # STAR alignment of endosperm RNA-seq
│   ├── 04_quantify_ase.py         # Allele-specific expression quantification
│   └── 05_classify_imprinting.R   # MEG/PEG classification and figures
├── config/
│   └── config.yml                 # All parameters and paths
├── data/
│   └── example/                   # Small example dataset
├── envs/
│   └── environment.yml            # Conda environment
├── run_pipeline.sh                # Single entry point
├── LICENSE
└── README.md


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
conda env create -f envs/environment.yml
conda activate imprinting
```
## Usage

All parameters are set in `config/config.yml`. To run the full pipeline on the example data:
```bash
bash run_pipeline.sh config/config.yml
```

To run individual stages:
```bash
bash scripts/01_align_parental.sh config/config.yml
bash scripts/02_call_snps.sh config/config.yml
bash scripts/03_align_endosperm.sh config/config.yml
python scripts/04_quantify_ase.py --config config/config.yml
Rscript scripts/05_classify_imprinting.R config/config.yml
```

**Outputs:** A tab-delimited table of MEGs and PEGs with expression ratios and statistical test results, and diagnostic figures in `results/figures/`.

## Dependencies

Managed via conda (`envs/environment.yml`). Key tools:

| Tool | Version | Purpose |
|------|---------|---------|
| BWA | 0.7.17 | Parental DNA-seq alignment |
| bcftools | 1.15 | SNP calling |
| STAR | 2.7.10 | Endosperm RNA-seq alignment |
| Python | 3.8+ | ASE quantification |
| R | 4.0.5 | Statistical analysis and figures |

## Data Access

Raw sequencing data and full results will be deposited upon publication. This repository provides code and a small example dataset for illustrative purposes.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgements
Willis Lab (Duke) / NSF / Dr. John Willis
Franks Lab (NCSU) / NSF / Miguel Flores, Dr. Robert Franks
