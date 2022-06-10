import sys, time, random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src import dataExploration, database, alignment

from collections import Counter


annotationFile = 'Outputs/SRR12464727_Annotation.csv'
testF = 'Inputs/SRR12464727.fastq'

#reads = dataExploration.getReads(testF)
# reads = 39,335,609
# time to read file -> 221 seconds ~ 3.6 mins


# 14,858 viral genomes
viralDB1 = 'viralDB/viral.1.genomic.fna'
viralDB2 = 'viralDB/viral.2.genomic.fna'
viralDB3 = 'viralDB/viral.3.genomic.fna'
viralDB4 = 'viralDB/viral.4.genomic.fna'
vGDB = {}
vGDB.update(dataExploration.getViralDB(viralDB1))
vGDB.update(dataExploration.getViralDB(viralDB2))
vGDB.update(dataExploration.getViralDB(viralDB3))
vGDB.update(dataExploration.getViralDB(viralDB4))

kLen = 10
kmers = database.getkmers(vGDB, kLen)

kmerBinning = True
subset = 0
count = 0
s = time.time()
if kmerBinning:
    reads = dataExploration.getReads(testF)
    scores = []
    print(time.time()-s)
    s1 = time.time()
    # TODO: make dictionary {read: {viralGenome:#hits, viralGenome:hits}, ....}
    hits = {}
    totHits = {}
    for r in reads:
        possibleVGs = []
        for i in range(kLen, len(r)-kLen):
            k = r[i: i+kLen]
            if k in kmers:
                possibleVGs.extend(kmers[k])
        hits[count] = len(possibleVGs)
        totHits = Counter(possibleVGs)
        #if count%10000==0: print(count)
        count +=1
    print(time.time()-s)
    print(time.time()-s1)
    top = dict(sorted(hits.items(), key=lambda item: item[1], reverse=True)[:10])
    # TODO : alignment for reads with most hits for viralGenomes they matched to
    for t in top:
        print(top[t], len(totHits[t]))
        for align in range(3):
            print(reads[t], totHits[t][align])

    #x = hits.values()
    #print(len(list(filter(lambda x: x>6000, sorted(x)))))
    #x2 = sorted(x, reverse=True)[:3000]

    #plt.hist(x)
    #plt.show()

    #plt.hist(x2)
    #plt.show()




align = False
if align:
    reads = dataExploration.getReads(testF)
    annotationDF = pd.DataFrame(np.zeros((len(reads), 4)), columns=['LocalAlignment', 'LA_Score', 'GlobalAlignment', 'GA_Score' ])
    annotationDF.to_csv(annotationFile)
    annotationDF = pd.read_csv(annotationFile)
    for i in range(15):
        randI = random.randint(25,len(reads))
        maxS = -1
        match = ''
        for v in vGDB:
            #score = alignment.globalAlign(reads[i], vGDB[v])
            score = alignment.localAlign(reads[randI], vGDB[v])
            maxS = max(score, maxS)
            if maxS == score:
                match = v
        annotationDF.at[randI, 'LocalAlignment'] = match
        annotationDF.at[randI, 'LA_Score'] = maxS
        print(randI, match, maxS)


    annotationDF.to_csv(annotationFile, index=False)