import numpy as np
from collections import Counter
from src import alignment
import time


def subsetReads(reads, kLen, kmers):
    hits = {}
    count = 0
    top = {}
    totHits = {}
    interval = len(reads)/10
    s = time.time()
    for r in reads:
        if count % int(interval) == 0:
            print(time.time()-s)
        hits[count] = 0
        totHits[count] = {}
        for i in range(kLen, len(r) - kLen):
            k = r[i: i + kLen]
            if k in kmers:
                #hits[count] += 1
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
    return top, totHits, hits



def getVirusMatches(top, totHits, reads, vGDB):
    # alignment for reads with most hits for viralGenomes they matched to
    virusesFound = {}
    viralNames = {}

    for t in top:
        top5VG = sorted(totHits[t].items(), key=lambda item: item[1], reverse=True)[:3]
        maxS = 0
        match = ''
        for a in range(len(top5VG)):
            viralGenomeName = top5VG[a][0]
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

