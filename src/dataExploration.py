# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose: read in reads and viral genomes from files

import os


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
        line = line.strip()[:150]
        if next:
            reads.append(line)
            next = False
            count +=1
        if line.startswith('@'):
            next = True
    print(len(reads))
    return reads


# get reads from file
# @param readF: file name containing reads
# @returns reads: list of reads from file
def getReadsSubset(readF, stepNum, step):
    reads = []
    next = False

    with open(readF, 'r') as f:
        lines = f.readlines()
    stepSize = int(len(lines)/stepNum)
    print('total number of reads', len(lines))

    count = 0
    for l in range(stepSize*(step-1), min(stepSize*step, len(lines))):
        line = lines[l].strip()
        if next:
            reads.append(line)
            next = False
            count +=1
        if line.startswith('@'):
            next = True
    print('Length of reads:', len(reads))
    return reads


# get reads from file
# @param readF: file name containing reads
# @returns reads: list of reads from file
def readstoFiles(readF, run, stepNum):
    next = False
    if not os.path.exists('Inputs/' + run + 'reads'):
        os.makedirs('Inputs/' + run + 'reads')

    rFile = 'Inputs/'+run+'reads/' + run

    with open(readF, 'r') as f:
        lines = f.readlines()
        print(len(lines)/4)
    stepSize = int(len(lines)/stepNum)
    print('total number of reads', len(lines))

    for step in range(stepNum):
        with open(rFile+ '_' +str(step) +'.txt', 'w') as f:
            for l in range(stepSize*(step-1), min(stepSize*step, len(lines))):
                line = lines[l].strip()
                if next:
                    f.write(line+'\n')
                    next = False
                if line.startswith('@'):
                    next = True
    return


def getReadsinFiles(file):
    with open(file, 'r') as f:
        lines = f.readlines()

    reads = []
    for l in lines:
        reads.append(l.strip())
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
