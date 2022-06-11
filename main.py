import sys, time, random, argparse, csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src import dataExploration, database, alignment

from collections import Counter

parser = argparse.ArgumentParser(description='')
parser.add_argument('--unzip-viralDB', metavar='unzipVDB', type=bool, default=False, help='unzip viral genome database')

args = parser.parse_args()

annotationFile = 'Outputs/SRR12464727_Annotation.csv'
testF = 'Inputs/SRR12464727.fastq'

#reads = dataExploration.getReads(testF)
# reads = 39,335,609
# time to read file -> 221 seconds ~ 3.6 mins


# 14,858 viral genomes
# TODO: make commands to unzip the ones in compressed folder
# TODO: make small input files
viralDB1 = 'viralDB/genomes/viral.1.genomic.fna'
viralDB2 = 'viralDB/genomes/viral.2.genomic.fna'
viralDB3 = 'viralDB/genomes/viral.3.genomic.fna'
viralDB4 = 'viralDB/genomes/viral.4.genomic.fna'
vGDB = {}
vGDB.update(dataExploration.getViralDB(viralDB1))
vGDB.update(dataExploration.getViralDB(viralDB2))
vGDB.update(dataExploration.getViralDB(viralDB3))
vGDB.update(dataExploration.getViralDB(viralDB4))

kLen = 10
kmers = database.getkmers(vGDB, kLen)

kmerBinning = True
subset = 0
s = time.time()
if kmerBinning:
    count = 0
    virusOutFile = 'Outputs/virusCount.csv'
    reads = dataExploration.getReads(testF)
    print(len(reads))
    scores = []
    print('Time to get reads: ', time.time()-s)
    s1 = time.time()
    hits = {}
    totHits = {}
    top = {}
    for r in reads:
        possibleVGs = []
        for i in range(kLen, len(r)-kLen):
            k = r[i: i+kLen]
            if k in kmers:
                possibleVGs.extend(kmers[k])
        hits[count] = len(possibleVGs)
        if len(possibleVGs) > 4000:
            top[count] = len(possibleVGs)
        totHits[count] = Counter(possibleVGs)
        count +=1
    print('Time to subset reads: ',time.time()-s1)
    print('Total time:  ', time.time()-s)
    s2 = time.time()
    #top = dict(sorted(hits.items(), key=lambda item: item[1], reverse=True)[:3000])
    # alignment for reads with most hits for viralGenomes they matched to
    virusesFound = {}
    print(len(top))
    for t in top:
        top3VG = sorted(totHits[t].items(), key=lambda item:item[1], reverse=True)[:3]
        allVG = list(totHits[t].items())
        maxS = 0
        match = ''
        # TODO: align all not just top 3
        for a in range(len(allVG)):
            viralGenomeName = allVG[t][a][0]
            genomeHits = allVG[t][a][1]
            score = alignment.globalAlign(reads[t], vGDB[viralGenomeName])
            maxS = max(maxS, score)
            if maxS == score:
                match = viralGenomeName.split(' ')[0]
                match = match.split('.')[0]
        if match in virusesFound:
            virusesFound[match]+=1
        else:
            virusesFound[match] = 1

    print('Time to align top reads: ', time.time()-s2)
    print('Total time elapsed: ', time.time()-s)
    print(virusesFound)
    with open(virusOutFile, 'w') as of:
        for key in virusesFound.keys():
            of.write("%s,%s\n"%(key,virusesFound[key]))


    x = hits.values()
    print(len(list(filter(lambda x: x>4000, sorted(x)))))
    x2 = top.values()

    plt.hist(x)
    plt.show()

    plt.hist(x2)
    plt.show()
    # TODO: check covid in viral genomicd db, look at length of top viruses & genomes
    # TODO: make pie chart, nested pie chart????
    # TODO: host chart, virus chart, taxonomy information
    # TODO: given argument can input already created file




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