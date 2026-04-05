#!/usr/bin/env python

####### DEPENDENCIES ###########################################

import sys
import os
import numpy as np
from timeit import default_timer as timer
import collections
from collections import OrderedDict
import re
import time
import logging

##########FILES and ARGUMENTS NEEDED for PROGRAM to RUN############
# Usage: python call_pairwise_snps_w_genes.py <snpfile> <outprefix> <minqual> <gtffile>
# Depth thresholds are computed automatically from data (min=10, max=mean+2SD per parent)

snpfile   = sys.argv[1]  # text file made with bcftools query, 8 columns
outprefix = sys.argv[2]
minqual   = sys.argv[3]  # minimum mapping quality for site
gtffile   = sys.argv[4]  # e.g. Mguttatusvar_IM767_887_v2.1.gene_exons.gff3

######################CODE#########################################

# Input file format (from bcftools query):
# contig  pos  ref  alt  QUAL  MQ  genotype:allelicdepth:depth  genotype:depth
# Chr_01  14205  T  G  219.035  60  0/1;1  1/1;11

minqual = int(minqual)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

missing    = ["./.", ".", ".|."]
hets       = ["0/1", "0|1", "1/2", "1|2", "0/2", "0|2"]
homozygous = ["1/1", "1|1", "0/0", "0|0", "2/2", "2|2"]
homoref    = ["0/0"]
homoalt    = ["1/1", "2/2"]
badhets    = ["0/2", "2/3"]

out1      = open(outprefix + '_genic_snps.bed', 'w')
out2      = open(outprefix + '_coding_sequence_snps.bed', 'w')
genecounts = open(outprefix + '_genic_snp_counts.txt', 'w')
cdscounts  = open(outprefix + '_coding_sequence_snp_counts.txt', 'w')
depths = open(outprefix + '_depths.txt', 'w')


def compute_depth_thresholds(variantfile, min_abs=10):
    """
    Compute per-parent depth thresholds from the query file.
    mindepth = min_abs (absolute minimum)
    maxdepth = min(mean_p1 + 2*sd_p1, mean_p2 + 2*sd_p2)
    Uses the lower of the two parents' upper thresholds (conservative).
    """
    depths_p1 = []
    depths_p2 = []

    with open(variantfile, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            try:
                d1 = int(cols[6].split(':')[1])
                d2 = int(cols[7].split(':')[1])
                depths_p1.append(d1)
                depths_p2.append(d2)
            except (IndexError, ValueError):
                continue

    if not depths_p1 or not depths_p2:
        raise ValueError("Could not compute depth thresholds: no valid depth data found.")

    mean_p1, sd_p1 = np.mean(depths_p1), np.std(depths_p1)
    mean_p2, sd_p2 = np.mean(depths_p2), np.std(depths_p2)

    max_p1 = int(mean_p1 + 2 * sd_p1)
    max_p2 = int(mean_p2 + 2 * sd_p2)

    computed_min = min_abs
    computed_max = min(max_p1, max_p2)

    logger.info('Parent 1 depth: mean=%.1f, SD=%.1f, upper threshold=%d' % (mean_p1, sd_p1, max_p1))
    logger.info('Parent 2 depth: mean=%.1f, SD=%.1f, upper threshold=%d' % (mean_p2, sd_p2, max_p2))
    logger.info('Using mindepth=%d, maxdepth=%d' % (computed_min, computed_max))
    
    depths.write('Parent\tMeanDepth\tSD\tUpperThreshold\n')
    depths.write('Parent1\t%.1f\t%.1f\t%d\n' % (mean_p1, sd_p1, max_p1))
    depths.write('Parent2\t%.1f\t%.1f\t%d\n' % (mean_p2, sd_p2, max_p2))

    return computed_min, computed_max


def make_gene_dict(genefile):

    outdict = {}
    genes = []

    with open(genefile, "r") as INFILE:
        for line in INFILE:
            cols = line.split()
            linestring = "".join(line)
            if "gene" in linestring:
                gene_name = linestring.split("\t")[8].split(';')[1]
                gene_name = gene_name.strip('Name="')
                gene_name = gene_name.strip('\n')
                genes.append(gene_name)
                scaffold = cols[0]
                begin = int(cols[3])
                end = int(cols[4])
                gene_range = range(begin, end + 1)
                for site in gene_range:
                    siteID = "\t".join([scaffold, str(site)])
                    outdict[siteID] = gene_name
            else:
                continue
    return (outdict, genes)


def make_cds_dict(genefile):

    outdict = {}

    with open(genefile, "r") as INFILE:
        for line in INFILE:
            cols = line.split()
            linestring = "".join(line)
            if "CDS" in linestring or "UTR" in linestring:
                gene_name = linestring.split("\t")[8].split(';')[1]
                gene = gene_name.strip('Parent="')
                gene_string = gene.split('.')
                gene_string = gene_string[0:2]
                gene_name = '.'.join(gene_string)
                scaffold = cols[0]
                begin = int(cols[3])
                end = int(cols[4])
                gene_range = range(begin, end + 1)
                for site in gene_range:
                    siteID = "\t".join([scaffold, str(site)])
                    outdict[siteID] = gene_name
            else:
                continue
    return outdict


def check_all_greater(lst, number):
    return all(element >= number for element in lst)


def check_all_lesser(lst, number):
    return all(element <= number for element in lst)


def call_homo_snps(ref, alt, gT):

    alleles = []
    hetType = 'None'

    if len(alt) == 1:
        for g in gT:
            if g.split('/')[0] == '0':
                alleles.append(ref)
            else:
                alleles.append(alt)
    elif len(alt) == 3:
        alt1 = alt.split(',')[0]
        alt2 = alt.split(',')[1]
        for g in gT:
            if g.split('/')[0] == '1':
                alleles.append(alt1)
            elif g.split('/')[1] == '2':
                alleles.append(alt2)
            else:
                alleles.append('X')

    return (alleles, hetType)


def make_parental_snp_catalog(variantfile, mindepth, maxdepth):

    out1.write('contig\tbedPosition\tPosition\tMatAllele\tPatAllele\tHetType\tGene\n')
    out2.write('contig\tbedPosition\tPosition\tMatAllele\tPatAllele\tHetType\tGene\n')

    fh = logging.FileHandler(outprefix + '.log')
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    no_sites    = 0
    gene_sites  = 0
    coding_sites = 0
    low_depth   = 0
    high_depth  = 0

    start_time = time.time()

    gene_dict = make_gene_dict(gtffile)
    site_dict = gene_dict[0]
    all_genes = gene_dict[1]
    logger.info('--- %s seconds to assemble gene dictionary ---' % (time.time() - start_time))
    genelist_time = time.time()

    cds_dict = make_cds_dict(gtffile)
    logger.info('--- %s seconds to assemble coding sequence dictionary ---' % (time.time() - genelist_time))

    scaffold_list = ['Chr_01']
    full_genes_seen = []
    coding_seq_seen = []

    new_scaffold_time = time.time()

    with open(variantfile, 'r') as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            cols = line.replace('\n', '').split('\t')
            contig = cols[0]
            pos    = int(cols[1])
            bedpos = pos - 1
            ref    = cols[2]
            alt    = str(cols[3])
            QUAL   = cols[4]
            mq     = int(cols[5])

            if contig not in list(OrderedDict.fromkeys(scaffold_list)):
                scaffold_list.append(contig)
                logger.info('--- %s seconds to complete %s ---' % (time.time() - new_scaffold_time, scaffold_list[-2]))
                new_scaffold_time = time.time()

            if QUAL == '.' or float(QUAL) < minqual: continue
            if mq < minqual: continue
            if ref == 'N': continue
            if len(ref) > 1: continue
            if len(alt) > 3: continue
            if len(alt) == 2: continue
            if '*' in alt: continue

            no_sites += 1

            gData = [cols[6], cols[7]]

            genotypes = [i.split(':')[0] for i in gData if i.split(':')[0] not in missing]

            if len(genotypes) < 2: continue
            if len(set(genotypes)) < 2: continue

            genotypes = [i for i in genotypes if i not in badhets]

            if len(genotypes) < 2: continue

            depths = [int(i.split(':')[1]) for i in gData]

            if not check_all_greater(depths, mindepth):
                low_depth += 1
                continue

            if not check_all_lesser(depths, maxdepth):
                high_depth += 1
                continue

            if len(alt) == 1:
                if len(set.intersection(set(genotypes), set(hets))) == 0:
                    allele_info = call_homo_snps(ref, alt, genotypes)
                    #print(contig, pos, ref, alt, gData, allele_info)
                else:
                    continue
            elif len(alt) == 3:
                if len(set.intersection(set(genotypes), set(hets))) == 0:
                    allele_info = call_homo_snps(ref, alt, genotypes)
                else:
                    continue
            else:
                continue

            alleles = allele_info[0]
            hT      = allele_info[1]

            if 'X' in alleles:
                continue

            alleles = '\t'.join(alleles)
            site    = "\t".join([contig, str(pos)])

            if site in site_dict.keys():
                gene_sites += 1
                gene1 = site_dict[site]
                full_genes_seen.append(gene1)
                out1.write(contig + '\t' + str(bedpos) + '\t' + str(pos) + '\t' + alleles + '\t' + hT + '\t' + gene1 + '\n')

            if site in cds_dict.keys():
                coding_sites += 1
                gene2 = cds_dict[site]
                coding_seq_seen.append(gene2)
                out2.write(contig + '\t' + str(bedpos) + '\t' + str(pos) + '\t' + alleles + '\t' + hT + '\t' + gene2 + '\n')

    full_gene_set        = set(full_genes_seen)
    no_unique_full_genes = len(list(full_gene_set))
    snps_per_gene        = collections.Counter(full_genes_seen)

    genes_zero = set(all_genes) - full_gene_set
    snps_per_gene.update({i: 0 for i in genes_zero})
    snps_per_gene = collections.OrderedDict(sorted(snps_per_gene.items()))

    coding_seq_set   = set(coding_seq_seen)
    no_coding_genes  = len(list(coding_seq_set))
    snps_per_coding_gene = collections.Counter(coding_seq_seen)

    for key, value in snps_per_gene.items():
        genecounts.write(key + '\t' + str(value) + '\n')

    for key, value in snps_per_coding_gene.items():
        cdscounts.write(key + '\t' + str(value) + '\n')

    minute_time  = float(time.time() / 60)
    minute_start = float(start_time / 60)
    logger.info('There are %d genes with distinguishing parental SNPs in genic regions' % no_unique_full_genes)
    logger.info('There are %d genes with distinguishing parental SNPs in coding regions' % no_coding_genes)
    logger.info('There are %d sites with low parental depth' % low_depth)
    logger.info('There are %d sites with high parental depth' % high_depth)
    logger.info('--- %d minutes to complete task ---' % (minute_time - minute_start))


# ── Entry point ────────────────────────────────────────────────────────────────

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

minqual = int(minqual)

logger.info('Computing depth thresholds from data...')
mindepth, maxdepth = compute_depth_thresholds(snpfile)

make_parental_snp_catalog(snpfile, mindepth, maxdepth)