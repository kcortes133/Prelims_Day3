# Author: Katherina Cortes
# Date: June 11, 2022
# Purpose: visualizations of labeled viruses

import numpy as np
import matplotlib.pyplot as plt


def filter(labels):
    labelsFiltered = {}
    labels = dict(sorted(labels.items(), key=lambda item: item[1], reverse=True))
    count = 0
    for l in labels:
        print(labels[l])
        if count < 9:
            labelsFiltered[l] = int(labels[l])
        else:
            if 'Other' in labelsFiltered:
                labelsFiltered['Other'] += int(labels[l])
            else: labelsFiltered['Other'] = int(labels[l])
        count +=1
    print(labelsFiltered)

    return labelsFiltered


def pieChart(labels):
    # TODO: make low count ones same slice
    viruses = list(labels.keys())
    counts = list(labels.values())
    fig, ax = plt.subplots(figsize=(6,3), subplot_kw=dict(aspect='equal'))

    def func(pct, allvals):
        absolute = int(np.round(pct/100.*np.sum(allvals)))
        return "{:.1f}%\n({:d} g)".format(pct, absolute)

    #wedges, texts, autotexts = ax.pie(counts, autopct=lambda pct: func(pct, counts),
    #                                  textprops=dict(color="w"))
    wedges, texts = ax.pie(counts, textprops=dict(color="w"))

    ax.legend(wedges, viruses,
              title="Viruses",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    plt.setp(texts, size=8, weight="bold")

    ax.set_title("Viruses present in sample")

    plt.show()


def pieFromFile(file):
    #file = 'Outputs/virusCountSRR12464727_all.csv'
    with open(file, 'r') as f:
        labels = {}
        lines = f.readlines()
        for l in lines:
            ls = l.strip().split(',')
            labels[ls[0]] = int(ls[2])

    filteredL = filter(labels)
    pieChart(filteredL)