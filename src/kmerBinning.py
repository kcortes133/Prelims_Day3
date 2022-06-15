# Author: Katherina Cortes
# Date: June 13, 2022
# Purpose: to filter out reads that have viral origins from the other
#   reads in a metagenomic sample by looking at kmers

from src import alignment
import time


# @param reads: list of all reads in sample
# @param kLen: int specifying kmer size
# @param kmers: dictionary of kmers of size kLen from viruses and corresponding virusIDS
# @returns top: dict of reads with over 6000 matching kmers from virus genomes and their count
# @returns totHits: dict of reads and viruses with shared kmers and counts
#
# get reads likely to be viral in origin
def subsetReads(reads, kLen, kmers):
    hits = {}
    count = 0
    top = {}
    totHits = {}
    interval = len(reads)/10
    s = time.time()
    for r in reads:
        hits[count] = 0
        totHits[count] = {}
        for i in range(kLen, len(r) - kLen):
            k = r[i: i + kLen]
            if k in kmers:
                hits[count] += len(kmers[k])
                if k in totHits:
                    for v in kmers[k]:
                        totHits[count][v]+=1
                else:
                    for v in kmers[k]:
                        totHits[count][v] = 1
        if hits[count] > 6000:
            top[count] = hits[count]
        count += 1
    return top, totHits


# @param top: dict of reads with over 6000 matching kmers from virus genomes and their count
# @param totHits: dict of reads and viruses with shared kmers and counts
# @param reads: list of all reads in sample
# @param vGDB: dictionary of virus IDs and names and theri genomic sequence
# @returns virusesFound: dictionary of virus IDs and names and counts
# @returns viralNames: dictionary of virus ID and corresponding name
def getVirusMatches(top, totHits, reads, vGDB):
    # alignment for reads with most hits for viralGenomes they matched to
    virusesFound = {}
    viralNames = {}

    for t in top:
        topVG = sorted(totHits[t].items(), key=lambda item: item[1], reverse=True)[:3]
        maxS = 0
        match = ''
        for a in range(len(topVG)):
            viralGenomeName = topVG[a][0]
            score = alignment.globalAlign(reads[t], vGDB[viralGenomeName])
            maxS = max(maxS, score)
            if maxS == score:
                match = viralGenomeName.split(' ')[0]
                match = match.split('.')[0]
                viralNames[match] = ''.join(viralGenomeName)
        if match in virusesFound:
            virusesFound[match] += 1
        else:
            virusesFound[match] = 1
    return virusesFound, viralNames

