# Author: Katherina Cortes
# Date: June 11, 2022
# Purpose: visualizations of labeled viruses

import matplotlib.pyplot as plt


# @params labels: dict of virus and their counts
# @params top: int, number of viruses to limit pie chart to
# @returns labelsFiltered: top viruses and counts, all other viruses added to Other
def filter(labels, top):
    labelsFiltered = {}
    labels = dict(sorted(labels.items(), key=lambda item: item[1], reverse=True))
    count = 0
    for l in labels:
        if count < top-1:
            labelsFiltered[l] = int(labels[l])
        else:
            if 'Other' in labelsFiltered:
                labelsFiltered['Other'] += int(labels[l])
            else: labelsFiltered['Other'] = int(labels[l])
        count +=1

    return labelsFiltered


# @param labels: dict of viruses and counts
# @displays: pie chart of viruses
def pieChart(labels, title):
    viruses = list(labels.keys())
    counts = list(labels.values())

    fig, ax = plt.subplots(figsize=(6,3), subplot_kw=dict(aspect='equal'))
    wedges, texts = ax.pie(counts, textprops=dict(color="w"))

    ax.legend(wedges, viruses,
              title="Viruses",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    plt.setp(texts, size=8, weight="bold")
    ax.set_title(title)
    plt.show()


# @params file: file name where virus count information is stored
# @params run: name of run to add to title of chart
# @displays: pie chart of viruses
def pieFromFile(file, top, title):
    with open(file, 'r') as f:
        labels = {}
        lines = f.readlines()
        for l in lines:
            ls = l.strip().split(',')
            labels[ls[0]] = int(ls[2])

    filteredL = filter(labels, top)
    pieChart(filteredL, title)
    return


# @params file: file name of virus information
# @params run: name of run to add to figure title
# @params top: highest count hosts to display
# @displays: pie chart of virus hosts and counts
# makes a pie chart of viral hosts and relative counts
def pieFromFile_Hosts(file, top, title):
    # get hosts and counts from file
    with open(file, 'r') as f:
        labels = {}
        lines = f.readlines()
        for l in lines:
            ls = l.strip().split(',')
            if ls[3] in labels:
                labels[ls[3]] += int(ls[2])
            else:
                labels[ls[3]] = int(ls[2])

    # only get highest count hosts
    filteredL = filter(labels, top)
    print(filteredL)
    pieChart(filteredL, title)
    return


def count(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        print(len(lines))
        c = 0
        for l in lines:
            ls = l.strip().split(',')
            c += int(ls[2])
    print(c)
    return c

count('Outputs/virusCountSRR12464727_all.csv')
count('Outputs/virusCountSRR12432009_all.csv')
