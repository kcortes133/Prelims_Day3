import numpy as np


# mismatch penalty of -1
# indel penalty of -4
# match addition of +1
# align one sequence to another
def globalAlign(vSeq, refSeq):
    # match, mismatch, indel
    scorM = [3,-2, -4]
    # rows, columns
    alignM = np.zeros((len(vSeq), len(refSeq)))

    alignM[0][0] = scorM[0]
    for j in range(1, len(refSeq)):
        if vSeq[0] == refSeq[j]:
            match = scorM[0]
        else:
            match = scorM[1]
        alignM[0][j] = max(0, match)

    for i in range(1, len(vSeq)):
        if vSeq[i] == refSeq[0]:
            match = scorM[0]
        else:
            match = scorM[1]
        alignM[i][0] = max(0, match)

    for i in range(1, len(vSeq)):
        for j in range(1, len(refSeq)):

            if vSeq[i] == refSeq[j]: match = scorM[0]
            else: match = scorM[1]

            alignM[i][j] = max(alignM[i-1][j]+scorM[2], alignM[i][j-1]+scorM[2], alignM[i-1][j-1] + match)
    return max(alignM[-1][:].max(), alignM[:][-1].max())

# mismatch penalty of -1
# indel penalty of -4
# match addition of +1
# ignores flanking areas
def localAlign(vSeq, refSeq):
    # match, mismatch, indel
    scorM = [1,-1, -4]
    alignM = np.zeros((len(vSeq), len(refSeq)))

    alignM[0][0] = scorM[0]
    for j in range(1, len(refSeq)):
        if vSeq[0] == refSeq[j]:
            match = scorM[0]
        else:
            match = scorM[1]
        alignM[0][j] = max(0, match)

    for i in range(1, len(vSeq)):
        if vSeq[i] == refSeq[0]:
            match = scorM[0]
        else:
            match = scorM[1]
        alignM[i][0] = max(0, match)

    for i in range(1, len(vSeq)):
        for j in range(1, len(refSeq)):

            if vSeq[i] == refSeq[j]: match = scorM[0]
            else: match = scorM[1]

            alignM[i][j] = max(alignM[i-1][j]+scorM[2], alignM[i][j-1]+scorM[2], alignM[i-1][j-1] + match, 0)

    return alignM.max()

s1 = 'cat'
s2 = 'ccatb'

#print(globalAlign(s1, s2))
print(localAlign(s1, s2))
