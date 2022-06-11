import numpy as np
import matplotlib.pyplot as plt


def pieChart(labels):
    viruses = list(labels.keys())
    counts = list(labels.values())
    print(viruses)
    print(counts)

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

    #bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    #kw = dict(arrowprops=dict(arrowstyle="-"),
    #          bbox=bbox_props, zorder=0, va="center")

    #for i, p in enumerate(wedges):
    #    ang = (p.theta2 - p.theta1) / 2. + p.theta1
    #    y = np.sin(np.deg2rad(ang))
    #    x = np.cos(np.deg2rad(ang))
    #    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    #    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    #    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    #    ax.annotate(viruses[i], xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
    #                horizontalalignment=horizontalalignment, **kw)

    plt.show()


file = 'Outputs/virusCount1.csv'

with open(file, 'r') as f:
    labels = {}
    lines = f.readlines()
    for l in lines:
        ls = l.strip().split(',')
        print(ls[0], ls[1])
        labels[ls[0]] = ls[1]

pieChart(labels)