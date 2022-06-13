# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose:


import time, random, argparse, os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src import dataExploration, database, alignment, plotGenerator, kmerBinning, outFiles
from collections import Counter
import random

parser = argparse.ArgumentParser(description='')
parser.add_argument('--kBinning', metavar='kBinning', type=bool, default=False, help='find reads likely to be'
                                                                                           ' viral and align to viral genomes')
parser.add_argument('--kLen', metavar='kLen', type=int, default=10, help='length of k for binning')
parser.add_argument('--align', metavar='align', type=bool, help='align all reads to viral genomes,'
                                                                ' only use if number of reads is small')
parser.add_argument('--readsFile', metavar='readsF', type=str, default='Example/Inputs', help='fastq file of reads')
parser.add_argument('--outputFile', metavar='outputF', type=str, default='Example/Outputs', help='output directory')
parser.add_argument('--visualize', metavar='visualize', type=bool, default=False, help='')
parser.add_argument('--viralSampleFile', metavar='vsFile', type=str, default='', help='create figures for already made'
                                                                                      ' viral community characterization file')

args = parser.parse_args()

def main():
    annotationFile = 'Outputs/SRR12464727_Annotation.csv'
    testF = 'Inputs/SRR12464727.fastq'
    run = 'SRR12464727'
    annotationFile = 'Outputs/SRR12432009_Annotation.csv'
    testF = 'Inputs/SRR12432009.fastq'
    run = 'SRR12432009'

    # TODO: make small example input and output files
    # TODO: add input and output arguments

    # unzip viral genome files
    # will do automatically if genomes file doesnt exist and
    # 14,858 viral genomes
    # unzip the viral genomes in compressed folder if not unzipped
    if not os.path.exists('viralDB/genomes'):
        os.makedirs('viralDB/genomes')
        vDir = 'viralDB/compressed/viral.'
        for i in range(1,5):
            vfile = vDir + str(i) + '.1.genomic.fna.gz'
            ovFile = vDir + str(i) + '.genomic.fna'
            database.extract(vfile, ovFile)

    vGDB = database.makeVDB()
    s = time.time()
    if args.kBinning:
        kmers = database.getkmers(vGDB, args.kLen)
        virusOutFile = 'Outputs/virusCount' + run +'_all.csv'
        s1 = time.time()
        top = {}
        totHits = {}
        hits = {}
        directory = 'Inputs/' + run + 'reads/'

        #if not os.path.exists(directory):
        #    dataExploration.readstoFiles(testF, run, 20)

        virusesFound = {}
        viralNames = {}
        reads = dataExploration.getReads(testF)
        reads = random.sample(reads, int(len(reads)/4))
        print(time.time() - s)
        top, totHits, hits = kmerBinning.subsetReads(reads, args.kLen, kmers)
        print('Time to subset reads: ',time.time()-s1)
        print(time.time()- s)
        virusesFound, viralNames = kmerBinning.getVirusMatches(top, totHits, reads, vGDB)

        #for fileName in os.listdir(directory):
        #    f = os.path.join(directory + fileName)
        #    print(f)
        #    reads = dataExploration.getReadsinFiles(f)
        #    topTemp, totHitsTemp, hitsTemp = kmerBinning.subsetReads(reads, args.kLen, kmers)
        #    top.update(topTemp)
        #    totHits.update(totHitsTemp)
        #    hits.update(hitsTemp)
        #    virusesFoundTemp, viralNamesTemp = kmerBinning.getVirusMatches(top, totHits, reads, vGDB)

        #    virusesFound.update(virusesFoundTemp)

        #    viralNames.update(viralNamesTemp)
        print('Possible Viruses found: ', len(top))
        print('Time to subset reads: ',time.time()-s1)
        print('Total time:  ', time.time()-s)
        s2 = time.time()
        # alignment for reads with most hits for viralGenomes they matched to
        #virusesFound, viralNames = kmerBinning.getVirusMatches()

        print('Time to align top reads: ', time.time()-s2)
        print('Total time elapsed: ', time.time()-s)
        viralHostsDB = database.getViralHosts()

        outFiles.writeViruses(virusOutFile, virusesFound, viralHostsDB, viralNames)
        plotGenerator.pieChart(virusesFound)

        # TODO: check covid in viral genomicd db, look at length of top viruses & genomes
        # TODO: given argument can input already created file

    if args.align:
        reads = dataExploration.getReads(testF)
        annotationDF = pd.DataFrame(np.zeros((len(reads), 4)), columns=['GlobalAlignment', 'GA_Score'])
        annotationDF.to_csv(annotationFile)
        annotationDF = pd.read_csv(annotationFile)
        for i in range(15):
            randI = random.randint(25,len(reads))
            maxS = -1
            match = ''
            for v in vGDB:
                score = alignment.globalAlign(reads[i], vGDB[v])
                maxS = max(score, maxS)
                if maxS == score:
                    match = v
            annotationDF.at[randI, 'Global_Alignment'] = match
            annotationDF.at[randI, 'GA_Score'] = maxS

        annotationDF.to_csv(annotationFile, index=False)


main()