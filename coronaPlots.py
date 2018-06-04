# -*- coding: utf-8 -*-
"""
Created on Fri May 26 15:03:53 2017

Plotting functions

@author: vgor
"""
import pylab as plt
import numpy as np
from collections import Counter
from matplotlib.patches import Rectangle
    
def intersectionPlot(overlappH, overlapTemp, colors):
    ax = plt.subplot(111)
    x0 = 0
    x1 = 0
    for i in range(5):
        ax.add_patch(Rectangle((x0, 0.5), overlappH[i, 0], 1, color = colors[-i]))
        ax.add_patch(Rectangle((x1, 2.5), overlapTemp[i, 0], 1, color = colors[-i]))
        ax.text(x0 + overlappH[i, 0] / 2, 1, "{:.0f}\n\n{:.1f}%".format(*overlappH[i, :]), va = "center", ha = "center", fontsize = 9)
        ax.text(x1 + overlapTemp[i, 0] / 2, 3, "{:.0f}\n\n{:.1f}%".format(*overlapTemp[i, :]), va = "center", ha = "center", fontsize = 9)
        ax.text(x0 + overlappH[i, 0] / 2, 0.45, "{}".format(overlappH.shape[0] - i), va = "top", ha = "center", fontsize = 10, color = colors[-i])
        ax.text(x1 + overlapTemp[i, 0] / 2, 2.45, "{}".format(overlapTemp.shape[0] - i), va = "top", ha = "center", fontsize = 10, color = colors[-i])
        x0 += overlappH[i, 0]
        x1 += overlapTemp[i, 0]
    
    plt.ylim(0, 4)
    plt.xlim(0, max(x0, x1)*1.05)
    plt.yticks([1, 3], ["pH", "Temperature"], va = "center", ha = "right")
    plt.text(-0.005*max(x0, x1), 0.45, "# of conditions", va = "top", ha = "right", fontsize = 10)
    plt.text(-0.005*max(x0, x1), 2.45, "# of conditions", va = "top", ha = "right", fontsize = 10)
    plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on', labelbottom='on')
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14, pad = 8)

def directionPlot(omniT, omniP):
    countsT = Counter(omniT["Direction"].astype(str))
    countsP = Counter(omniP["Direction"].astype(str))
    
    colors = {"Down": "r", "nan": "y", "Up": "g"}
    ax = plt.subplot(111)
    x0 = 0
    x1 = 1
    for label in ["Down", "nan", "Up"]:
        ax.add_patch(Rectangle((0.5, x0), 1.5, countsT[label], color = colors[label]))
        ax.add_patch(Rectangle((3.5, x1), 1.5, countsP[label], color = colors[label]))
        ax.text(1.25, x0 + countsT[label] / 2, label, va = "center", ha = "center", fontsize = 12)
        ax.text(4.25, x1 + countsP[label] / 2, label, va = "center", ha = "center", fontsize = 12)
        x0 += countsT[label]
        x1 += countsP[label]
    
    plt.ylim(0, max(omniT.shape[0], omniP.shape[0])*1.05)
    plt.xlim(0, 5.5)
    plt.xticks([1.25, 4.25], ["Temperature", "pH"], va = "center", ha = "center")
    plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on', labelbottom='on')
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14, pad = 8)

def proteinAbundancePlot(abundanceData, phData, tempData):
    #add number column for plotting
    abundanceData["N"] = range(abundanceData.shape[0])

    #mark protein present
    abundanceData["inTemp"] = abundanceData.index.map(lambda x: x in tempData["GeneID"].values).values.astype(bool)
    abundanceData["inPH"] = abundanceData.index.map(lambda x: x in phData["GeneID"].values).values.astype(bool)
    
    #indices for existing entries
    indBoth = np.logical_and(abundanceData["inTemp"], abundanceData["inPH"])
    indT = np.logical_and(abundanceData["inTemp"], ~ abundanceData["inPH"])
    indP = np.logical_and(~ abundanceData["inTemp"], abundanceData["inPH"])
    
    plt.semilogy(abundanceData["N"], abundanceData[("conc", "cavg")], "-", ms = 2, color = "k", alpha = 0.5, label = "Reported concentration")
    plt.semilogy(abundanceData.loc[indT, "N"], abundanceData.loc[indT, ("conc", "cavg")], "o", ms = 4, label = "Temperature", alpha = 0.4)
    plt.semilogy(abundanceData.loc[indP, "N"], abundanceData.loc[indP, ("conc", "cavg")], "o", ms = 4, label = "pH", alpha = 0.4)
    plt.semilogy(abundanceData.loc[indBoth, "N"], abundanceData.loc[indBoth, ("conc", "cavg")], "o", ms = 4, label = "Both", alpha = 0.4)
    plt.xlabel("Gene", fontsize = 14)
    plt.ylabel("Concentration in plasma, mg/mL", fontsize = 14)
    sel = np.linspace(0, abundanceData.shape[0] - 1, 25).astype(int)
    plt.xticks(sel, abundanceData.index[sel], va = "top", ha = "center", rotation = "vertical")
    plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on', labelbottom='on')
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14, pad = 8)
    plt.legend(fontsize=12)
