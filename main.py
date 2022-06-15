# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose: get information about the viral community in a metagenomics sample


import time, argparse, os

import numpy as np
import pandas as pd
from src import dataExploration, database, alignment, plotGenerator, kmerBinning, outFiles
import random

parser = argparse.ArgumentParser(description='')
parser.add_argument('--kBinning', metavar='kBinning', type=bool, default=False, help='find reads likely to be'
                                                                                           ' viral and align to viral genomes')
parser.add_argument('--kLen', metavar='kLen', type=int, default=10, help='length of k for binning')
parser.add_argument('--align', metavar='align', type=bool, default=False, help='align all reads to all viral genomes,'
                                                                'NOT RECOMMENDED, especially for large files'
                                                                ' only use if number of reads is small')
parser.add_argument('--readsFile', metavar='readsF', type=str, default='Example/Input/SRR12464727_example.fastq', help='fastq file of reads')
parser.add_argument('--run', metavar='run', type=str, default='SRR12464727_example', help='run ID')
parser.add_argument('--randSubset', metavar='rSub', type=int, default=4, help='fraction of random reads to use for analysis')
parser.add_argument('--visualize', metavar='visualize', type=bool, default=False, help='Visualize viral community from file')
parser.add_argument('--viralSampleFile', metavar='vsFile', type=str, default='Example/Output/virusCountSRR12464727_example_all.csv',
                    help='create figures for already made viral community characterization file')

args = parser.parse_args()

def main():
    # unzip viral genome files
    # will do automatically if genomes file doesnt exist and
    # 14,858 viral genomes
    # unzip the viral genomes in compressed folder if not unzipped
    # make file to write outputs to
    if not os.path.exists('Outputs/'):
        os.makedirs('Outputs')

    # unzip viral genome database
    if not os.path.exists('viralDB/genomes'):
        os.makedirs('viralDB/genomes')
        vDir = 'viralDB/compressed/viral.'
        uvDir = 'viralDB/genomes/viral.'
        for i in range(1,5):
            vfile = vDir + str(i) + '.1.genomic .fna.gz'
            ovFile = uvDir + str(i) + '.genomic.fna'
            database.extract(vfile, ovFile)

    # make viral genome data base into dictionary
    vGDB = database.makeVDB()
    s = time.time()

    # get viral community information using kmers
    if args.kBinning:
        # get kmers from viruses
        kmers = database.getkmers(vGDB, args.kLen)
        virusOutFile = 'Outputs/virusCount' + args.run +'_all.csv'

        # get reads from file
        reads = dataExploration.getReads(args.readsFile)
        # randomly subset reads
        reads = random.sample(reads, int(len(reads)/args.randSubset))
        print(len(reads))
        print('Tme to get reads from file: ',time.time() - s)
        # get reads likely to be from viral genomes

        top ={}
        totHits ={}
        stepS = int(len(reads)/10)
        for i in range(10):
            iterReads = reads[stepS*i:stepS*(1+i)-1]
            topT, totHitsT = kmerBinning.subsetReads(iterReads, args.kLen, kmers)
            print('Time to subset' + str(i+1) + '/10 of reads: ', time.time()-s)
            top.update(topT)
            totHits.update(totHitsT)
        print('Time to subset reads: ',time.time()-s)
        # get viruses that have high alignment score to reads
        virusesFound, viralNames = kmerBinning.getVirusMatches(top, totHits, reads, vGDB)

        print('Possible Viruses found: ', len(top))
        print('Total time:  ', time.time()-s)
        viralHostsDB = database.getViralHosts()

        # write viral information to file
        outFiles.writeViruses(virusOutFile, virusesFound, viralHostsDB, viralNames)
        plotGenerator.pieChart(virusesFound)

    # make pie charts about virus count and host count
    if args.visualize:
        vtitle = 'Viruses in ' + args.run
        htitle = 'Hosts in ' + args.run
        plotGenerator.pieFromFile(args.viralSampleFile, 10, vtitle)
        plotGenerator.pieFromFile_Hosts(args.viralSampleFile, 10, htitle)

    # get viruses by aligning every read to every viral genome
    # not recommended unless you have a small set of reads
    if args.align:
        reads = dataExploration.getReads(args.readsFile)
        annotationDF = pd.DataFrame(np.zeros((len(reads), 4)), columns=['GlobalAlignment', 'GA_Score'])
        annotationDF.to_csv(args.outputFile)
        annotationDF = pd.read_csv(args.outputFile)
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

        annotationDF.to_csv(args.outputFile, index=False)


main()