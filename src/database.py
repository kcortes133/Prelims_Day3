# Author: Katherina Cortes
# Date: June 10, 2022
# Purpose: get viral information except genome

import gzip, shutil
import pandas as pd
from src import dataExploration



def extract(zipF, unzipF):
    with gzip.open(zipF, 'rb') as f_in:
        with open(unzipF, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return


def getkmers(vDB, kLen):
    kmers = {}
    for v in vDB:
        vG = vDB[v]
        for i in range(kLen, len(vG) - kLen):
            k = vG[i:i+kLen]
            if k in kmers:
                kmers[k].append(v)
            else: kmers[k] = [v]
    return kmers


def getViralHosts():
    vFile = 'viralDB/taxid10239.nbr'
    hostDB = pd.read_table(vFile, comment='#', delimiter='\t', names=["Representative","Neighbor","Host","Selected lineage","Taxonomy name","Segment name"])
    hostDB = hostDB[['Representative', 'Host']]
    return dict(zip(hostDB['Representative'], hostDB['Host']))


def getViralInfoDB():
    vFile = 'viralDB/virushostdb.tsv'
    vDB = pd.read_table(vFile, delimiter='\t')
    return vDB


def makeVDB():
    viralDB1 = 'viralDB/genomes/viral.1.genomic.fna'
    viralDB2 = 'viralDB/genomes/viral.2.genomic.fna'
    viralDB3 = 'viralDB/genomes/viral.3.genomic.fna'
    viralDB4 = 'viralDB/genomes/viral.4.genomic.fna'
    vGDB = {}
    vGDB.update(dataExploration.getViralDB(viralDB1))
    vGDB.update(dataExploration.getViralDB(viralDB2))
    vGDB.update(dataExploration.getViralDB(viralDB3))
    vGDB.update(dataExploration.getViralDB(viralDB4))
    return vGDB


