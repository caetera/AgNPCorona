
# coding: utf-8

# #### Setup the environment
import numpy as np
from pickle import load, dump
from multiprocessing import Pool
from pandas import DataFrame, read_csv
from scipy.optimize import curve_fit
from scipy.stats import norm

def gaussian(x, height=1, position=0, width=1):
    return height * np.exp(-1 * np.power(x - position, 2) / (2 * width**2))

def shrinkage(sample1, sample2):
    span1 = np.percentile(sample1, [10, 90])
    span2 = np.percentile(sample2, [10, 90])

    return(span1[1] - span1[0]) / (span2[1] - span2[0])

def permutTest(sample1, sample2, function = None, ntest = 100000):
    """
    Permutation test
    """
    if function is None:
        function = lambda a, b: np.median(a) - np.median(b)

    sampleSize = len(sample1)
    effectSize = function(sample1, sample2)

    population = np.concatenate((sample1, sample2))

    effects = []
    for n in range(ntest):
        np.random.shuffle(population)
        effects.append(function(population[:sampleSize], population[sampleSize:]))

    yt, xt = np.histogram(effects, bins=100)
    xt = (xt[:-1] + xt[1:])/2
    opt = curve_fit(gaussian, xt, yt, p0=(4000, np.mean(effects), np.std(effects)))[0]

    if effectSize > opt[1]:
        return 2*norm(*opt[1:]).sf(effectSize)
    else:
        return 2*norm(*opt[1:]).cdf(effectSize)

def runTest(args):
    pcData, mask, featureID = args
    return permutTest(pcData.loc[mask, featureID].values,
                      pcData.loc[~mask, featureID].values, shrinkage)


# #### Prepare for the test
workers = Pool(24)
pcData = load(open("./pcData.pickle", "rb"))
groupIndex = read_csv("./omniIndex.csv")


# #### Run the test
result = DataFrame(index = pcData.columns)

for yy in ["pH", "Temp"]:
    mask = groupIndex[yy].values

    result[yy] = workers.map(runTest, [(pcData, mask, ID) for ID in pcData.columns])


# #### Cleaning
workers.close()
dump(result, open("./permResultShrink.pickle", "wb"))
