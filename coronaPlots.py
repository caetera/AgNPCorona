# -*- coding: utf-8 -*-
"""
Created on Fri May 26 15:03:53 2017

Plot functions

@author: vgor
"""
import pylab as plt
import numpy as np
from pandas import DataFrame, Series
from collections import defaultdict
from sklearn.decomposition import PCA

def scaleLabels(x, y, xs = 0, ys = 1, num = 10):
    a = (x - y) / (xs - ys)
    b = x - a * xs
       
    return a * np.linspace(0, 1, num) + b

def cumulativeBarplot(hist, bins, colors, names, **kwargs):
    """
    Create barplot with data series stacked up and legend reversed and located outside the plot
    hist - dataframe with columns corresponding to bars and rows corresponding to series
    colors - list of color codes
    names - list of data series names
    """
    if len(hist) != len(colors) or len(hist) != len(names):
        raise ValueError("Number of rows, length of colors and length of names should be all equal")
        
    bottom = hist.apply(np.cumsum, axis = 0)
    
    plt.bar(range(hist.shape[1]), hist.iloc[0, :], color = colors[0], label = names[0], align = "center")
    for i in range(1, len(hist)):
        plt.bar(range(hist.shape[1]), hist.iloc[i, :], bottom = bottom.iloc[i -1, :], color = colors[i],\
                    label = names[i], align = "center")
    
    #drawing means
    if not bins is None:
        means = (hist.T * (bins[:-1] + bins[1:])/2).apply(np.sum, axis = 1).values
        if kwargs.has_key("means"):
            low = min(means.min(), kwargs["means"].min())
            high = max(means.max(), kwargs["means"].max())
            plt.plot(range(1, hist.shape[1], 3), (kwargs["means"] - (high + low)/2)/(high - low) * 0.7 + 0.5,\
                lw = 3, ms = 0, color = "#A02C2C", label = "Control")
        else:
            low = means.min()
            high = means.max()
        
        plt.plot((means - (high + low)/2)/(high - low) * 0.7 + 0.5,\
                lw = 3, ms = 0, color = "#002255", label = "Sample")
        
        if kwargs.has_key("scaleY"):
            ylabels = ["{:.2f}".format(z/kwargs["scaleY"]) for z in scaleLabels(low, high, 0.15, 0.85, 10)]
        else:
            ylabels = ["{:.2f}".format(z) for z in scaleLabels(low, high, 0.15, 0.85, 10)]
        plt.yticks(np.linspace(0,1,10), ylabels, va = "center")
        
    plt.ylim(0, 1)
    plt.xlim(-0.6, hist.shape[1] - 0.4)
    plt.xticks(range(hist.shape[1]), hist.columns, rotation = "vertical", ha = "center")
    plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on')
    plt.tick_params(axis = 'both', which = 'major', labelsize = 18, pad = 8)
    handles, labels = plt.gca().get_legend_handles_labels()
    newHandles = handles[:1:-1] + [plt.Line2D([0], [0], linewidth = 0)] + handles[1::-1]
    newLabels = labels[:1:-1] + ["Average values"] + labels[1::-1]
    plt.legend(newHandles, newLabels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 18)
    
    return means
    
def pIPlot(mdc, colnames, nBins = 10, colors = None, **kwargs):
    """
    Create plot for pI values, using relative abundances as weights
    mdc - connected data for protein abundances and protein characteristics
    colnames - columns to include
    nBins - number of bins to split the data into
    colors - color to be used (if defined should be not less than nBins)
    """
    hist = DataFrame()
    bins = np.linspace(mdc["pI"].min(), mdc["pI"].max(), nBins)
    for colname in colnames:
        col = mdc.loc[np.isfinite(mdc[colname]), :]
        hist[colname] = np.histogram(col["pI"], bins = bins, weights = col[colname]/col[colname].sum())[0]
    
    if colors is None:
        colors = [plt.get_cmap("hsv")(1.*i/len(hist)) for i in range(len(hist))]
    
    names = ["{:.2f} - {:.2f}".format(a, b) for (a,b) in zip(bins[:-1], bins[1:])]
    plt.title("Protein pI\n")
    return cumulativeBarplot(hist, bins, colors, names, **kwargs)

def MWPlot(mdc, colnames, nBins = 10, colors = None, **kwargs):
    """
    Create plot for MW values on log scale, using relative abundances as weights
    mdc - connected data for protein abundances and protein characteristics
    colnames - columns to include
    nBins - number of bins to split the data into
    colors - color to be used (if defined should be not less than nBins)
    """
    hist = DataFrame()
    bins = np.exp(np.linspace(np.log(mdc["MW"].min()), np.log(mdc["MW"].max()), nBins))
    for colname in colnames:
        col = mdc.loc[np.isfinite(mdc[colname]), :]
        hist[colname] = np.histogram(col["MW"], bins = bins, weights = col[colname]/col[colname].sum())[0]

    if colors is None:
        colors = [plt.get_cmap("hsv")(1.*i/len(hist)) for i in range(len(hist))]
    
    names = ["{:.0f} - {:.0f}".format(a, b) for (a,b) in zip(bins[:-1]/1000., bins[1:]/1000.)]
    plt.title("Protein MW (kDa)\n")
    return cumulativeBarplot(hist, bins, colors, names, scaleY = 1000, **kwargs)

def GravyPlot(mdc, colnames, nBins = 10, colors = None, **kwargs):
    """
    Create plot for GRAVY values , using relative abundances as weights
    mdc - connected data for protein abundances and protein characteristics
    colnames - columns to include
    nBins - number of bins to split the data into
    colors - color to be used (if defined should be not less than nBins)
    """
    hist = DataFrame()
    bins = np.linspace(mdc["GRAVY"].min(), mdc["GRAVY"].max(), nBins)
    for colname in colnames:
        col = mdc.loc[np.isfinite(mdc[colname]), :]
        hist[colname] = np.histogram(col["GRAVY"], bins = bins, weights = col[colname]/col[colname].sum())[0]

    if colors is None:
        colors = [plt.get_cmap("hsv")(1.*i/len(hist)) for i in range(len(hist))]
    
    names = ["{:.2f} - {:.2f}".format(a, b) for (a,b) in zip(bins[:-1], bins[1:])]
    plt.title("Protein GRAVY\n")
    return cumulativeBarplot(hist, bins, colors, names, **kwargs)
        
def AAPlot(mdc, colnames, nBins = 10, colors = None):
    """
    Create amino acid content plot, using relative abundances as weights
    mdc - connected data for protein abundances and protein characteristics
    colnames - columns to include
    nBins - number of bins to split the data into
    colors - color to be used (if defined should be not less than nBins)
    """
    df = DataFrame()
    for colname in colnames:
            col = mdc.loc[np.isfinite(mdc[colname]), :]
            weights = col[colname]/col[colname].sum()
            AA = defaultdict(float)
            for i in col.index:
                for aa, content in col.loc[i, "AAC"].items():
                    AA[aa] += weights[i] * content
            
            fs = Series(AA)
            fs = fs/fs.sum()
            df[colname] = fs.copy()
            
    if colors is None:
        colors = [plt.get_cmap("hsv")(1.*i/nBins) for i in range(nBins)]
        colors = (colors*3)[::3] + [(1.0, 1.0, 1.0)]
        
    df.fillna(0, inplace = True)
    df.sort(df.columns[0], ascending = False, inplace = True)
    df = df[:nBins]        
    df.loc["Other", :] = 1 - df.sum(axis = 0)
    cumulativeBarplot(df, None, colors[:nBins + 1], df.index)
    plt.title("Protein AA Content\n")
    plt.ylim(0, (1 - np.min(df.loc["Other", :])) * 1.05)

def GOPlot(mdc, colnames, goColumn, nBins = 10, colors = None):
    """
    Create plot for GO terms, using relative abundances as weights
    Only GO terms with relative total abundance greater that 0.01 are plotted
    mdc - connected data for protein abundances and protein characteristics
    colnames - columns to include
    goColum - GO column to use
    nBins - number of GO terms (the most abundant) to preserve
    colors - color to be used (if defined should be not less than the number of GO terms
            after the filtering + 1)
    """    
    df = DataFrame()
    for colname in colnames:
            col = mdc.loc[np.isfinite(mdc[colname]), :]
            weights = col[colname]/col[colname].sum()
            functions = defaultdict(float)
            for i in col.index:
                for term in col.loc[i, goColumn]:
                    functions[term] += weights[i]
            
            fs = Series(functions)
            fs = fs/fs.sum()
            df[colname] = fs.copy()
    
    if colors is None:
        colors = [plt.get_cmap("hsv")(1.*i/nBins) for i in range(nBins)]
        colors = (colors*3)[::3] + [(1.0, 1.0, 1.0)]
        
    df.fillna(0, inplace = True)
    df.sort(df.columns[0], ascending = False, inplace = True)
    df = df[:nBins]
    df.loc["Other", :] = 1 - df.sum(axis = 0)
    cumulativeBarplot(df, None, colors[:nBins + 1], df.index)
    plt.title("Protein {}\n".format(goColumn))
    plt.ylim(0, (1 - np.min(df.loc["Other", :])) * 1.05)

def PCAplot(data, labels):
    """
    Create PCA plot
    """
    pC = PCA(2)
    pca_data = pC.fit_transform(data[data.apply(lambda x: ~np.any(np.isnan(x)), axis = 1)].T)
    
    for i in range(pca_data.shape[0] / 3):
        plt.plot(pca_data[i*3:i*3+3, 0], pca_data[i*3:i*3+3, 1], marker = "o", lw = 0, label = labels[i], ms = 10)
    
    plt.legend(loc = "lower left", fontsize = 18, numpoints = 1)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 18, pad = 8)
    plt.xlabel("PC1 ({:.1f} %)".format(pC.explained_variance_ratio_[0] * 100), fontsize = 18)
    plt.ylabel("PC2 ({:.1f} %)".format(pC.explained_variance_ratio_[1] * 100), fontsize = 18)
    
