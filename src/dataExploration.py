# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose: read in files

import matplotlib.pyplot as plt
import numpy as np


# get reads from file
# @param readF: file name containing reads
# @returns reads: list of reads from file
def getReads(readF):
    reads = []
    next = False

    with open(readF, 'r') as f:
        lines = f.readlines()
    print(len(lines))

    count = 0
    for line in lines:
        line = line.strip()
        if next:
            reads.append(line)
            next = False
            count +=1
        if line.startswith('@'):
            next = True
        #if count == 3000000:
        #    break
    return reads

# get reads from file
# @param readF: file name containing reads
# @returns reads: list of reads from file
def getViralDB(readF):
    reads = {}
    label = ''
    next = False

    with open(readF, 'r') as f:
        lines = f.readlines()


    for line in lines:
        line = line.strip()
        if next:
            reads[label] = line
            next = False
        if line.startswith('>'):
            next = True
            label = line.strip('>')

    return reads
