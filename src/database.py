import gzip, shutil
import pandas as pd

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
            if len(k) != kLen:
                print('wtf', k)
            if k in kmers:
                kmers[k].append(v)
            else: kmers[k] = [v]
    print(len(kmers))
    return kmers


def getViralHosts():
    vFile = 'viralDB/taxid10239.nbr'
    hostDB = pd.read_table(vFile, comment='#', delimiter='\t', names=["Representative","Neighbor","Host","Selected lineage","Taxonomy name","Segment name"])
    hostDB = hostDB[['Representative', 'Host']]
    return hostDB

def getViralInfoDB():
    vFile = 'viralDB/virushostdb.tsv'
    vDB = pd.read_table(vFile, delimiter='\t')
    print(vDB)
    return vDB


getViralHosts()
getViralInfoDB()