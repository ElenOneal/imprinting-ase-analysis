#!/usr/bin/env python

# Count reads assigned to parent1 and parent2 per gene using pysam pileup.
# Counts at most one SNP per read per gene to avoid double-counting.
# Genes with zero counts for both parents are not written to output.
#
# Based on Count_Variants.py by Stefan Wyder (2015-09-01).
# Modified by E. Oneal.
#
# Usage: python classify_parental_alleles.py <bamfile> <snpfile> <outprefix>
#
# snpfile columns: chrom  bedpos  pos  parent1_allele  parent2_allele  gene
# outprefix: prefix for output files (<outprefix>_genecounts.txt, <outprefix>.log)

from sys import argv
import pysam
from collections import defaultdict
import time
import logging

##########################################################
# Quality cutoffs
##########################################################

MIN_BASE_QUAL   = 20
MIN_MAPPING_QUAL = 60

##########################################################
# Arguments
##########################################################

BamFile   = argv[1]
SnpFile   = argv[2]
outprefix = argv[3]

##########################################################
# Logging
##########################################################

logging.basicConfig(level=logging.DEBUG, format='%(message)s')
logger = logging.getLogger()
fh = logging.FileHandler(outprefix + '.log')
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)

start_time = time.time()

##########################################################
# Data structures
##########################################################

# Use sets instead of lists for O(1) lookup and lower memory usage
ReadCounted  = defaultdict(set)
SNPpositions = []
SNPinfo      = {}
parent1Reads = {}
parent2Reads = {}

##########################################################
# Read SNP file
##########################################################

with open(SnpFile, 'r') as snpPos:
    for line in snpPos:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        Chrom, bedpos, Pos, Parent1Allele, Parent2Allele, HetType, OverlappingGene = line.split('\t')
        SNPpositions.append([Chrom, int(bedpos)])
        SNPinfo[(Chrom, int(bedpos))] = line.split('\t')

logger.info('--- %.1f seconds to read SNP file (%d positions) ---' % (
    time.time() - start_time, len(SNPpositions)))

##########################################################
# Pileup over BAM file
##########################################################

samfile = pysam.AlignmentFile(BamFile, 'rb')

for QueryChrom, QueryPos in SNPpositions:
    OverlappingGene = ''
    for pileupcolumn in samfile.pileup(
            QueryChrom, QueryPos, QueryPos + 1,
            truncate=True, max_depth=10000, stepper='samtools'):

        if pileupcolumn.pos != QueryPos:
            continue

        # Positions are 0-based in pysam
        fields = SNPinfo.get((QueryChrom, QueryPos), ['', '', '', '', '', '', ''])
        Chrom, bedpos, Pos, Parent1Allele, Parent2Allele, HetType, OverlappingGene = fields[:7] + [''] * (7 - len(fields[:7]))

        countparent1  = 0
        countparent2  = 0
        basecount = {}

        for pileupread in pileupcolumn.pileups:
            Alignment = pileupread.alignment
            ReadName  = Alignment.query_name

            if pileupread.query_position is None:
                continue

            BaseQual    = ord(Alignment.qual[pileupread.query_position]) - 33
            MappingQual = Alignment.mapping_quality

            if (not pileupread.is_del
                    and not pileupread.is_refskip
                    and not Alignment.is_secondary
                    and ReadName not in ReadCounted[OverlappingGene]
                    and BaseQual    >= MIN_BASE_QUAL
                    and MappingQual >= MIN_MAPPING_QUAL):

                Allele = Alignment.query_sequence[pileupread.query_position]

                if Allele == Parent2Allele:
                    countparent2 += 1
                    ReadCounted[OverlappingGene].add(ReadName)
                elif Allele == Parent1Allele:
                    countparent1 += 1
                    ReadCounted[OverlappingGene].add(ReadName)

                # Track base counts at this position
                basecount[Allele] = basecount.get(Allele, 0) + 1

        parent1Reads[OverlappingGene] = parent1Reads.get(OverlappingGene, 0) + countparent1
        parent2Reads[OverlappingGene] = parent2Reads.get(OverlappingGene, 0) + countparent2

    logger.info('--- %.1f seconds to count %s pos %s gene %s ---' % (
        time.time() - start_time, QueryChrom, QueryPos, OverlappingGene))

samfile.close()
logger.info('--- %.1f seconds to complete pileup ---' % (time.time() - start_time))

##########################################################
# Read full gene list from SNP file
##########################################################

all_genes = set()
with open(SnpFile, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        all_genes.add(fields[6])  # OverlappingGene is 7th column

##########################################################
# Write output — include 0s for genes with no counts
##########################################################

with open(outprefix + '_genecounts.txt', 'w') as outfile:
    outfile.write('\t'.join(['Gene', 'parent1_reads', 'parent2_reads']) + '\n')
    for gene in sorted(all_genes):
        parent1 = parent1Reads.get(gene, 0)
        parent2 = parent2Reads.get(gene, 0)
        outfile.write('\t'.join([gene, str(parent1), str(parent2)]) + '\n')
logger.info('--- %.1f seconds total ---' % (time.time() - start_time))
