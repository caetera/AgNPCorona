# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:42:16 2017

Data processing and figures

@author: vgor
"""
import re
import numpy as np
import pylab as plt
from coronaHelpers import readPercOut, cleanIDs, getProteinData, parseIndex
from coronaPlots import proteinAbundancePlot, intersectionPlot, directionPlot, lesseningDirectionPlot
from pandas import DataFrame, Series, read_csv, concat
from collections import Counter
from pyteomics.pylab_aux import scatter_trend
from scipy.optimize import curve_fit
from pyteomics import fasta
from pickle import dump, load
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess


#Working directories
workDir = "E:\\RawData\\20161120_NP\\"
resDir = workDir + "res\\"
fastaDir = "E:\\Fasta\\"

#Mapping UNIMOD ID to trivial name
unimodMap = {"1": "Acetyl",
             "4": "Carbamidomethyl",
             "5": "Carbamyl",
             "35": "Oxidation"}

def parseSamples(name):
    """
    Parse sample information  from  OpenMS peptide file
    (Corresponndence between abundance and sample)
    """
    with open("{}{}.peptides.csv".format(workDir, name)) as fin:
        lineNr = 0
        while lineNr < 2:
            fin.readline()
            lineNr += 1
        desc = fin.readline()
        pairs = re.findall(r"(\d+):.+?\\(\w+)_\d+\.", desc)
        return pairs

def quantifyProteins(name):
    """
    TopN Quantition using peptides from OpenMS
    name - the name of the file to work with
    """
    def quantifyProtein(peptidelist, workcolumns, aggFunc = np.nanmedian, N = 3):
        """
        Perform the quantification for individual protein
        aggFunc - function used for aggreagtion (ex. mean)
        N as in TopN
        """
        data = iqpeptides.reindex(index=peptidelist)
    
        result = Series()
    
        for c in workcolumns:
            topN = data.sort_values(c, ascending = False)[:N]
            if np.all(np.isfinite(topN[c])):
                result[c] = aggFunc(topN[c])
    
        return result
    
    def countPeptides(peptidelist, workcolumns):
        """
        Count peptides with quantification for individual protein
        """
        data = iqpeptides.reindex(index=peptidelist)
    
        result = Series()
    
        for c in workcolumns:
            result[c] = np.sum(np.isfinite(data[c]))
    
        return result


    #load percolator results
    psms, peptides, proteins = readPercOut("{}{}.pout.xml".format(workDir, name))
    
    #load OpenMS peptides results
    qpeptides = read_csv("{}{}.peptides.csv".format(workDir, name), sep = "\t", skiprows = 3)
    
    #load OpenMS column descrptions
    pairs = {"abundance_" + k: v for k,v in parseSamples(name)}
    
    #rename the columns
    qpeptides.columns = [pairs[c] if c in pairs.keys() else c for c in qpeptides.columns]
    
    #convert sequences to OpenMS notation
    peptides["Sequence"] = peptides["ID"].apply(lambda s: re.sub(r"\[UNIMOD:(\d+)\]",\
                             lambda i: "({})".format(unimodMap.get(i.groups()[0])), s))
    
    #filter by q-value and map the peptide qunatifications
    iqpeptides = peptides[peptides["q-value"] < 0.01].merge(qpeptides, left_on = "Sequence",\
                 right_on = "peptide", how = "inner")
    
    #remove unnecessary columns
    iqpeptides.drop(["Sequence", "peptide", "n_proteins", "charge", "PSMs", "CalcMass", "ExpMass", "protein"],\
                     axis = 1, inplace = True)
    
    #reindex using peptide ID
    iqpeptides.set_index("ID", drop = True, inplace = True, verify_integrity = True)
    
    #identified proteins (filter by q-value)
    iproteins = proteins[proteins["q-value"] < 0.01].copy()
    
    #quantify proteins
    iproteins = concat([iproteins, iproteins["Peptides"].apply(quantifyProtein, args = (sorted(pairs.values()),))], axis = 1) 
    
    #count quantified peptides per protein
    proteinsCount = iproteins["Peptides"].apply(countPeptides, args = (sorted(pairs.values()),)) 
    
    #find non quantified proteins (all Abundances are NaN)
    qind = np.any(np.isfinite(iproteins[list(pairs.values())].values), axis = 1)
    iproteins.drop(["Peptides", "q-value"], axis = 1, inplace = True)
    iproteins = iproteins[qind]
    
    #report information about identifications and quantifications
    print("""
          {}
          Identified peptides: {}
          Quantified peptides: {}
          Identified proteins: {}
          Quantified proteins: {}
          Median quantified peptides per protein: {}
          """.format(name, sum(peptides["q-value"] < 0.01), len(iqpeptides), sum(proteins["q-value"] < 0.01),
                      len(iproteins), np.median(proteinsCount.values)))
    
    return iproteins

def processExperiment(name):
    """
    Do the qunatification on a single experiments
    """
    #collect data from all replicates and contol
    d = quantifyProteins("{}".format(name))
    
    #names of columns with numerical data
    datacolumns = d.columns[1:]
    
    #log10 transformation, cleaning, saving the result into plain csv
    #protein ID (as reported by search engine vs quantification results)
    d[datacolumns] = d[datacolumns].apply(np.log10, axis = 1)
    d.replace(-np.inf, 0, inplace = True)
    d.set_index("ID", verify_integrity = True, inplace = True)
    d.to_csv("{}{}_raw.csv".format(workDir, name))      

def calcCoverage(proteinSequence, peptides):
    """
    Calculate protein coverage in percent and number of peptides
    
    NOTE: peptides are provided as IDs from percolator, i.e. contain modification
    information - [UNIMOD:(digits)]
    """
    #found aminoacids
    found = np.repeat(0, len(proteinSequence))
    
    numpeptides = 0
    
    for peptide in peptides:
        numpeptides += 1
        for match in re.finditer(re.sub(r"\[UNIMOD:(\d+)\]", "", peptide), proteinSequence):
            found[match.start():match.end()] = 1
    
    return found.mean() * 100, numpeptides

def exportProteinData(resData, name):
    """
    Export data for Protein table (supplementary material for paper)
    """
    #calculate protein data
    res = resData.copy(deep=True)
    res["CleanID"] = [cleanIDs(x) for x in res.index]
    mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
    
    #calculate coverage and number of peptides
    #load percolator results
    peptides, proteins = readPercOut("{}{}.pout.xml".format(workDir, name))[1:]
    
    #convert percolator IDs to UniprotIDs
    proteins["CleanID"] = [cleanIDs(x) for x in proteins["ID"]]
    #merging peptide data
    mdc = mdc.merge(proteins, how="inner", on="CleanID")
    
    #identified peptides
    ipeptidesIDs = set(peptides.loc[peptides["q-value"] < 0.01, "ID"])
    
    #add coverage and number of peptides
    mdc[["Coverage(%)", "#Peptides"]] = DataFrame(mdc.apply(lambda r: calcCoverage(r["Sequence"],\
       filter(lambda pID: pID in ipeptidesIDs, r["Peptides"])), axis=1).values.tolist())
    
    #prepare for export
    mdc["GOComponent"] = mdc["GOComponent"].str.join(", ")
    mdc["GOFunction"] = mdc["GOFunction"].str.join(", ")
    mdc["GOProcess"] = mdc["GOProcess"].str.join(", ")
    
    datacols = list(filter(lambda s: s.startswith(name[:4]), mdc.columns))
    newcols = ["UniprotID", "GeneID"] + datacols +\
    ["Name", "GOComponent", "GOFunction", "GOProcess", "Coverage(%)", "#Peptides", "MW", "pI", "GRAVY", "q-value"]
    
    mdc[newcols].to_excel("{}{}_allProteins.xlsx".format(resDir, name), index=False)

def getTopN(resT, labels, N):
    """
    Extract N most abundant proteins for each conditon
    resT - DataFrame as returned by the quantification routine
    labels - column names to process
    N - number of top proteins to extract
    """
    #collect indices of rows in the dataframe
    dd = set()
    for label in labels:
        dd = dd.union(set(resT.sort_values(label, ascending = False).index[:N]))
    
    return resT.loc[list(dd), :].fillna(0)

def sigmoid(x, x0, k):
    """
    Sigmoid function definition
    """
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return y
     
def get_sigmoid_fit(row, x):
    """
    Perform fittinig with the sigmoid function on one row (protein amount in 5 conditions)
    x - is defined externally and should be scaled to [0, 1]
    Returns float number, with absolute value corresponding to the intercept (in scaled coordiantes)
    and sign corresponding to the direction of change
    """
    y = row[:5].values.astype(float)
    y = (y - min(y)) / (max(y) - min(y)) #scale y to [0, 1]
    try:
        opt, cov = curve_fit(sigmoid, x, y, [0.5, 0.5 if y[0] < y[-1] else -0.5])
        #if the fit is good  -> return the fit
        if 1 - np.power((sigmoid(x, *opt) - y), 2).sum() / np.power((y - y.mean()), 2).sum() > 0.5:
            return opt[0] * np.sign(opt[1])
        #bad fit -> return nothing
        else:
            return np.nan
    except:
        return np.nan

def plot_sigmoid_fit(x, y, ystd, xlab = "", ylab = ""):
    """
    Create the plot of sigmoid fit
    x, y - coordinates
    ystd - standard deviation of y (for error bars)
    xlab, ylab - labels for x and y 
    """
    #scale x and y to [0, 1]
    x_scaled = (x - x.min()) / (x.max() - x.min())
    y_scaled = (y - y.min()) / (y.max() - y.min())
    
    plt.errorbar(x_scaled, y_scaled, yerr = ystd / (y.max() - y.min()), fmt = "bo")
    #perform fitting
    try:
        opt, cov = curve_fit(sigmoid, x_scaled, y_scaled, [0.5, 0.5 if y[0] < y[-1] else -0.5])
        #Total and unexplained variance
        varTot = np.power((y_scaled - y_scaled.mean()), 2).sum()
        varErr = np.power((sigmoid(x_scaled, *opt) - y_scaled), 2).sum()
        #plotting
        plt.plot(np.linspace(0, 1, 50), sigmoid(np.linspace(0, 1, 50), *opt), "g-")
        plt.axvline(opt[0], ls = "--", color = "r")
        plt.text(opt[0]+0.02, 1.0, "{:.1f}".format(x.min() + opt[0] * (x.max() - x.min())), color = "r")
        plt.text(0.9, 0.5, "$R^2 = {:.2f}$".format(1 - varErr/varTot))
    except Exception:
        plt.text(0.5, 0.5, "Fit not found", color = "k", ha = "center")
    plt.xticks(x_scaled, x)
    plt.yticks([0, 0.5, 1], ["{:.2f}".format(z) for z in [y.min(), y.mean(), y.max()]])
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.1, 1.1)
    plt.xlabel(xlab)
    plt.ylabel(ylab)


#%%Analysis
#Perform quantification
processExperiment("Temperature")
    
processExperiment("pH")

#%%Perform the analysis with control subtraction
#pH experiments
#collect the quan data and merging
md = read_csv("{}pH_raw.csv".format(workDir))
md.set_index("ID", drop = True, inplace = True)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
md.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, md.index)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#join technical replicates
res = DataFrame(index = md.index)
for n in ["pH{}_{}".format(c,r) for c in ["3", "5", "6_6", "7", "9"] for r in ["B1", "B2", "B3", "C"]]:
        res[n] = \
            md.loc[:, ["{}T{}".format(n, i) for i in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#export protein data
exportProteinData(res, "pH")

#subtract protein levels in control samples
res.fillna(0, inplace = True)
for c in ["3", "5", "6_6", "7", "9"]:
    for r in ["B1", "B2", "B3"]:
        res["pH{}_{}".format(c, r)] = res["pH{}_{}".format(c, r)] - res["pH{}_C".format(c)]
    
    res.drop("pH{}_C".format(c), axis = 1, inplace = True)

#remove negative values (after subtraction)
res[res <= 0] = np.nan
res = res[res.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#retrieve the properties and map them to proteins
res["CleanID"] = [cleanIDs(x) for x in res.index]
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#Summarize the data at the condition level (i.e. join biological replicates)
#Join biological replicates
res2 = DataFrame(index = res.index)
for c in ["3", "5", "6_6", "7", "9"]:
    res2["pH{}".format(c)] = \
            res.loc[:, ["pH{}_{}".format(c, r) for r in ["B1", "B2", "B3"]]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#extract UniprotID
res2["CleanID"] = [cleanIDs(x) for x in res2.index]

#Quantification overlap
overlap = Counter(np.sum(np.isfinite(res2[res2["CleanID"] != ""].drop("CleanID", axis = 1)), axis = 1))
overlappH = np.array([[overlap[k], float(overlap[k])/(res2["CleanID"] != "").sum() * 100] for k in sorted(overlap.keys(), reverse = True)])

#retrieve the properties and map them to proteins
mdc2 = res2[res2["CleanID"] != ""].merge(getProteinData(res2.loc[res2["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc2.drop(["CleanID", "Sequence"], axis = 1, inplace = True)
#save the results for later use
mdc2.to_csv("{}pH_conditionDetail.csv".format(resDir), index = False)


#Temperature experiments
#collect the quan data
md = read_csv("{}Temperature_raw.csv".format(workDir))
md.set_index("ID", drop = True, inplace = True)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
md.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, md.index)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#join technical replicates
res = DataFrame(index = md.index)
for n in ["Temp{}{}".format(c,r) for c in "ABCDE" for r in "123C"]:
        res[n] = \
            md.loc[:, ["{}_{}".format(n, i) for i in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#export protein data
exportProteinData(res, "Temperature")

#subtract protein levels in control samples
res.fillna(0, inplace = True)
for c in "ABCDE":
    for r in "123":
        res["Temp{}{}".format(c, r)] = res["Temp{}{}".format(c, r)] - res["Temp{}C".format(c)]
    
    res.drop("Temp{}C".format(c), axis = 1, inplace = True)

#remove negative values (after subtraction)
res[res <= 0] = np.nan
res = res[res.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#retrieve the properties and map them to proteins
res["CleanID"] = [cleanIDs(x) for x in res.index]
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#Join biological replicates
res2 = DataFrame(index = res.index)
for c in "ABCDE":
    res2["Temp{}".format(c)] = \
            res.loc[:, ["Temp{}{}".format(c, r) for r in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#extract UniprotID
res2["CleanID"] = [cleanIDs(x) for x in res2.index]

#Quantification overlap
overlap = Counter(np.sum(np.isfinite(res2[res2["CleanID"] != ""].drop("CleanID", axis = 1)), axis = 1))
overlapTemp = np.array([[overlap[k], float(overlap[k])/(res2["CleanID"] != "").sum() * 100] for k in sorted(overlap.keys(), reverse = True)])

#retrieve the properties and map them to proteins
mdc2 = res2[res2["CleanID"] != ""].merge(getProteinData(res2.loc[res2["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc2.drop(["CleanID", "Sequence"], axis = 1, inplace = True)
#save data for later
mdc2.to_csv("{}\\Temp_conditionDetail.csv".format(resDir), index = False)

#Intersection of quantified proteins
colors = ["#cf7b43", "#7d6acb", "#75b140", "#c459b6", "#54a676", "#cc4b42", "#6797d0", "#a5903e", "#c45d83"]

plt.figure(figsize = (9,6))
intersectionPlot(overlappH, overlapTemp, colors)    
plt.savefig("{}QProtIntersection.svg".format(resDir), bbox_inches = "tight")

#%%Identify differentially abundant proteins
#Performing fit with sigmoid curve
def calcSigmoids(mdc2, xr, label):
    #scale x before perfoming fit
    x = (xr - xr.min())/ (xr.max() - xr.min())

    #fit itself
    mdc2["SigmaFit"] = mdc2[mdc2.iloc[:, :5].apply(lambda z: np.all(np.isfinite(z)), axis = 1)].apply(get_sigmoid_fit, args = (x,),  axis = 1)
    mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "Direction"] = mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "SigmaFit"].apply(lambda t: "Up" if t > 0 else "Down")

    mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "Crit"+label] = mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "SigmaFit"].apply(lambda t: xr[0] + abs(t) * (xr[-1] - xr[0]))

#pH experiment
#Real pH values
x = np.array([4.9, 6.1, 6.8, 7.7, 8.9])

#calculating fits
calcSigmoids(mdc2, x, "pH")
mdc2[np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)].to_csv("{}omnipHS.csv".format(resDir), index = False)
mdc2[np.isfinite(mdc2["SigmaFit"])].to_csv("{}sigmoidpHS.csv".format(resDir), index = False)

#plotting sigmoid fits
cInd = dict(zip(mdc2.columns[:5], [[t.startswith(s) for t in mdc.columns] for s in mdc2.columns[:5]]))
for i, row in mdc2.loc[np.isfinite(mdc2["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc2.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "pH", "log10(Abundance)")
    plt.savefig("{}SigmaFitspHS\\{}.svg".format(resDir, row["GeneID"]))
    plt.close()

#Plotting flat fits
#select proteins quantified in all conditions
mdc3 = mdc2[np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)]

cInd = dict(zip(mdc3.columns[:5], [[t.startswith(s) for t in mdc.columns] for s in mdc3.columns[:5]]))
for i, row in mdc3.loc[~np.isfinite(mdc3["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc3.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "pH", "log10(Abundance)")
    plt.savefig("{}SigmaFitspHS-\\{}.svg".format(resDir, row["GeneID"]))
    plt.close()


#Temeperature experiment
#Real temperature values
x = np.array([4., 17., 30., 41., 47.])

#calculating fits
calcSigmoids(mdc2, x, "Temp")
mdc2[np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)].to_csv("{}omniTempS.csv".format(resDir), index = False)
mdc2[np.isfinite(mdc2["SigmaFit"])].to_csv("{}sigmoidTempS.csv".format(resDir), index = False)

#plottining sigmoid fits
cInd = dict(zip(mdc2.columns[:5], [[t.startswith(s) for t in mdc.columns] for s in mdc2.columns[:5]]))
for i, row in mdc2.loc[np.isfinite(mdc2["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc2.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "Temperature", "log10(Abundance)")
    plt.savefig("{}SigmaFitsTempS\\{}.svg".format(resDir, row["GeneID"]))
    plt.close()

#plotting flat fits
#select proteins quantified in all conditions
mdc3 = mdc2[np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)]

cInd = dict(zip(mdc3.columns[:5], [[t.startswith(s) for t in mdc.columns] for s in mdc3.columns[:5]]))
for i, row in mdc3.loc[~np.isfinite(mdc3["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc3.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "Temperature", "log10(Abundance)")
    plt.savefig("{}SigmaFitsTempS-\\{}.svg".format(resDir, row["GeneID"]))
    plt.close()

#direction in persistent fraction
omniT = read_csv("{}omniTempS.csv".format(resDir))
omniP = read_csv("{}omnipHS.csv".format(resDir))

plt.figure(figsize = (9,6))
directionPlot(omniT, omniP)    
plt.savefig("{}DirectionPersistent.svg".format(resDir), bbox_inches = "tight")


#%%Top20
mdc2 = read_csv("{}pH_conditionDetail.csv".format(resDir))
getTopN(mdc2, ["pH3", "pH5", "pH6_6", "pH7", "pH9"], 20).to_csv("{}Top20pH.csv".format(resDir), index  = False)

mdc2 = read_csv("{}Temp_conditionDetail.csv".format(resDir))
getTopN(mdc2, ["TempA", "TempB", "TempC", "TempD", "TempE"], 20).to_csv("{}Top20Temp.csv".format(resDir), index  = False)


#%%TPP and pI correlation
#Load ThermoProteomeProfiling data
TS4 = read_csv("{}TableS4.csv".format(resDir))[["Protein NAME", "meltP_Vehicle_Expt2"]]
#load critcal temperature
sigmT = read_csv("{}sigmoidTempS.csv".format(resDir))

#idenitidy overlapping proteins
overlap = TS4[np.isfinite(TS4["meltP_Vehicle_Expt2"])].merge(sigmT[np.isfinite(sigmT["SigmaFit"])], left_on = "Protein NAME", right_on = "GeneID")

#plot correlation
plt.figure(figsize = (8,5))
scatter_trend(overlap["CritTemp"], overlap["meltP_Vehicle_Expt2"],)
plt.xlabel("Critical Temperature")
plt.ylabel("Melting Temperature")
plt.savefig("{}meltPTemp.svg".format(resDir))

#pI to critical pH correlation
mdc2 = read_csv("{}sigmoidpHS.csv".format(resDir))
plt.figure(figsize = (8,5))
selector = np.isfinite(mdc2["SigmaFit"])
scatter_trend(mdc2.loc[selector, "CritpH"], mdc2.loc[selector, "pI"])
plt.xlabel("Critical pH")
plt.ylabel("pI")
plt.savefig("{}pIpH.svg".format(resDir))

#%%Plasma Proteome plot
ppConcentration = read_csv("{}plasma_protein_concentrations.csv".format(resDir))

massAbr = { "g": 1,
           "mg": 1e-3,
           "ug": 1e-6,
           "µg": 1e-6,
           "ng": 1e-9,
           "pg": 1e-12,
           "fg": 1e-15}

volAbr = { "l": 1,
          "dl": 1e-1,
          "ml": 1e-3}


def parseConcentration(s):
    """
    Convert different formats of concentration to common one
    """
    if s.find("±") != -1:
        med, dev = [float(x) for x in s.split("±")]
        return np.array([med - dev, med + dev])
    elif s.find("-") != -1:
        low, high =[ float(x) for x in s.split("-")]
        return np.array([low, high])
    elif s.find("+") != -1:
        med, dev = [float(x) for x in s.split("+")]
        return np.array([med - dev, med + dev])
    else:
        return np.array([float(s)])
    
#parse units into coefficients
units = ppConcentration["Unit"].str.split("/", expand = True).apply\
           (lambda u: massAbr.get(u[0], np.nan) / volAbr.get(u[1].lower(), np.nan), axis = 1)

#remove lines with unparsable unit (it is only one, BTW)
ppConcentration = ppConcentration[np.isfinite(units)]

#make all concentation comparable
ppConcentration["conc"] = ppConcentration["concentration"].apply(parseConcentration) * units

#summarize different concentration
#pandas 0.20+ syntax 
def cmin(x):
    return np.concatenate(x.values).min()

def cmax(x):
    return np.concatenate(x.values).max()

def cavg(x):
    return np.concatenate(x.values).mean()

#group by Gene
abundanceData = ppConcentration.groupby("Gene Symbol").agg({"conc": [cmin, cmax, cavg]})
abundanceData.sort_values(("conc", "cavg"), ascending = False, inplace = True)

#read data from experiments
tempData = read_csv("{}Temp_conditionDetail.csv".format(resDir))
phData = read_csv("{}pH_conditionDetail.csv".format(resDir))

#all proteins
plt.figure(figsize = (9, 6))
proteinAbundancePlot(abundanceData, phData, tempData)
plt.savefig("{}plasmaConteration_all.svg".format(resDir))

#persistent proteins
plt.figure(figsize = (9, 6))
proteinAbundancePlot(abundanceData, \
                     phData[np.all(np.isfinite(phData.loc[:, "pH3":"pH9"]), axis = 1)], \
                     tempData[np.all(np.isfinite(tempData.loc[:, "TempA":"TempE"]), axis = 1)])
plt.savefig("{}plasmaConteration_omni.svg".format(resDir))


#%%Plasma Proteome Database building
accRe = re.compile(r"\w+\|(\w+)\|.+")

proteinIDs = [accRe.match(protein[0]).groups()[0] for protein in fasta.read("{}PPD.fasta".format(fastaDir))]

#update calmodulin
proteinIDs.remove("P62158")
proteinIDs.append("P0DP23")

protData = getProteinData(proteinIDs)

dump(protData, open("{}PPD_data.pickle".format(resDir), "wb"))


#%%PhysChem
featureTable, featureDescriptions = parseIndex("{}AAindex.txt".format(resDir))
#featureTable has only AA indices (544 x 20 matrix)
#featureDescriptions has description and reference information (544 x 2)
#both indixed by feature code

#get protein data
protData = load(open("{}PPD_data.pickle".format(resDir), "rb"))

#persistent proteins
phData = read_csv("{}omnipHS.csv".format(resDir))
tData = read_csv("{}omniTempS.csv".format(resDir))

protData["pH"] = protData["UniprotID"].apply(lambda s: s in phData["UniprotID"].values)
protData["Temp"] = protData["UniprotID"].apply(lambda s: s in tData["UniprotID"].values)

#create amino acid composition table for each protein, ommiting U
#since it is not present in Kyoto index
#Nproteins x 20 matrix
df = protData[["UniprotID", "AAC"]]
aacs = DataFrame(df["AAC"].tolist()).drop("U", axis = 1).fillna(0)[featureTable.columns]

#calculate 544 features for each protein
#Nproteins x 544 matrix
Xdata = DataFrame(np.dot(aacs.values, featureTable.values.T) / aacs.values.sum(axis = 1).reshape(-1, 1),\
          index = df["UniprotID"], columns = featureTable.index)

#Permutation test is running as a separate multithreaded process
#the code is in shrinkTest.py
#it results in a p-value estimation for each of 544 protein features in both experiments
#Preparing the input for it
dump(Xdata, open("{}pcData.pickle".format(resDir), "wb"))#proteins properties

protData[["pH", "Temp"]].to_csv("{}omniIndex.csv".format(resDir), index = False)#persistent proteins index

#load permutation test results
shrink = load(open("{}permResultShrink.pickle".format(resDir), "rb"))

for case in ["pH", "Temp"]:
    Ydata = protData[case].values
    #shrinkage mesure
    spanCorona = np.percentile(Xdata.loc[Ydata], [10,90], axis = 0)
    spanAll = np.percentile(Xdata, [10,90], axis = 0)
    medCorona = np.median(Xdata.loc[Ydata], axis = 0)
    medAll = np.median(Xdata, axis = 0)
    shrink["Lessening({})".format(case)] = (spanCorona[1,] - spanCorona[0,])/(spanAll[1,] - spanAll[0,])
    shrink["Difference({})".format(case)] = (medCorona - medAll)/(spanAll[1,] - spanAll[0,])
    #Benjamini-Hochberg correction
    shrink["FDR({})".format(case)] = multipletests(shrink[case], method="fdr_bh")[1]
    #special column used by Cytoscape to color only significant nodes
    shrink["Color({})".format(case)] = 1.0
    shrink.loc[shrink["FDR({})".format(case)] < 0.005, "Color({})".format(case)] = \
        shrink.loc[shrink["FDR({})".format(case)] < 0.005, "Lessening({})".format(case)]
    
    shrink.rename({case: "p-value({})".format(case)}, axis = "columns", inplace = True)

#reordering columns
shrink = shrink[["{}({})".format(a, b) for b in ["pH", "Temp"]
                                        for a in ["Difference", "Lessening", "p-value", "FDR", "Color"]]]

#Prepare the table with node information for Cytoscape
classNames = {"A": "Alpha and turn propensities",
              "B": "Beta propensity",
              "C": "Composition",
              "H": "Hydrophobicity",
              "O": "Other properties",
              "P": "Physicochemical properties",
              "nan": "Undefined"}

#reading protein properties from Tomii et al.
featureClasses = read_csv("{}AAindex_classes.txt".format(resDir), sep = "  ", header = None, engine = "python")
featureClasses.columns = ["Class", "Name", "Description"]
featureClasses["Class"] = featureClasses["Class"].str.slice(0, 1)
featureClasses.drop("Description", axis=1, inplace = True)
featureClasses.index = featureClasses["Name"]

#mapping the properties to results
shrink["Class"] = featureClasses.reindex(shrink.index)["Class"].apply(lambda s: classNames[str(s)])

concat([shrink, featureDescriptions], axis = 1).to_excel("{}shrinkTest.xlsx".format(resDir))

#Data for minimum spanning tree
#pairwise distances between all 544 protein properties
dd = 1 - np.abs(np.corrcoef(featureTable.values))
cyImport = DataFrame([[featureTable.index[ind1], featureTable.index[ind2], dd[ind1, ind2]]
             for ind1 in range(dd.shape[0]) for ind2 in range(dd.shape[1]) if ind1 < ind2],
                columns = ["Source", "Target", "Distance"])

cyImport["Correlation"] = 1 - cyImport["Distance"]

cyImport.to_csv("{}AAcorr.csv".format(resDir), index = False)


#refined assesment of significantly changing properties
#shrinkdata contains column with refined classes for significantly changing properties
#and column M that is 1 or -1
#  1 - direct changes - i.e. bigger value corresponds to increase in property;
# -1 - reverse changes - i.e. bigger value corresponds to decrease in property
shrinkdata = read_csv("{}multiplier.csv".format(resDir))
shrinkdata.index = shrinkdata["ID"]
shrink["Class"] = shrinkdata["Class"]
#M includes only alpha, turn, and beta propensities that are significantly changing
shrink["M"] = shrinkdata["M"]

#for hydrophobicity M can be calculated knowing the difference between index of I (hydrophobic)
#and D (hydrophilic)
ind = np.logical_and(shrink["Class"] == "Hydrophobicity", shrink["FDR(Temp)"] < 0.005)
shrink.loc[ind, "M"] = (((featureTable.loc[ind, "I"] - featureTable.loc[ind, "D"]) > 0).astype(int) - 0.5)*2

plt.figure(figsize=(9,6))
lesseningDirectionPlot(shrink)
plt.savefig("{}Direction_Lessening.svg".format(resDir))

#plots for span comparison
for experiment in ["Temp", "pH"]:
    fdrCol = "FDR({})".format(experiment)
    lessCol = "Lessening({})".format(experiment)
    
    topFeatures = shrink[shrink[fdrCol] < 0.005]
    Ydata = protData[experiment].values
    
    #boxplots
    for feature, row in topFeatures.iterrows():
        plt.figure(figsize = (8, 6))
        plt.boxplot([Xdata.loc[:, feature], Xdata.loc[Ydata, feature]], whis = [10, 90])
        plt.title(featureDescriptions.loc[feature, "Description"], ha = "center", va = "bottom", 
                  fontsize = 12)
        plt.figtext(0.5, 0.5, "Lessening = {:.2f}\n\np-val(BH) = {:.2e}".format(row[lessCol],\
                                                                                row[fdrCol]),\
                                                        va="center", ha="center", fontsize = 12)
        plt.xticks([1, 2], ["Background", "Persistent"])
        plt.tick_params(axis = 'both', which = 'major', labelsize = 12, pad = 8)
        plt.tight_layout()
        if row[lessCol] > 1:
            plt.savefig("{}lessening{}/broad_boxplot_{}.png".format(resDir, experiment, feature))
        else:
            plt.savefig("{}lessening{}/narrow_boxplot_{}.png".format(resDir, experiment, feature))
        
        plt.close()
    
    #histograms with smoothing
    for feature, row in topFeatures.iterrows():
        plt.figure(figsize = (8, 6))
        dAll = Xdata.loc[:, feature].values.reshape(-1, 1)
        dBound = Xdata.loc[Ydata, feature].values.reshape(-1, 1)
    
        values, edges = np.histogram(dAll, 50, density=True)
        centers = (edges[1:] + edges[:-1])/2
        lowessFit = lowess(values, centers, is_sorted=True, frac=0.2, it=0)
        plt.plot(lowessFit[:, 0], lowessFit[:, 1], label="Background")
        values, edges = np.histogram(dBound, edges,density=True)
        lowessFit = lowess(values, centers, is_sorted=True, frac=0.2, it=0)
        plt.plot(lowessFit[:, 0], lowessFit[:, 1], label="Persistent")
        plt.title(featureDescriptions.loc[feature, "Description"], ha = "center", va = "bottom")
        plt.figtext(0.1, 0.5, "Lessening = {:.2f}\n\np-val(BH) = {:.2e}".format(row[lessCol],\
                                                                                row[fdrCol]),\
                                                            va="center", ha="left", fontsize = 12)
        
        plt.title(featureDescriptions.loc[feature, "Description"], ha = "center", va = "bottom", fontsize = 12)
        plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off',
                        top='off', labelleft='off', labelbottom='on')
        plt.tick_params(axis = 'both', which = 'major', labelsize = 12, pad = 8)
        plt.tight_layout()
        plt.legend(fontsize = 12)
        
        if row[lessCol] > 1:
            plt.savefig("{}lessening{}/broad_hist_{}.png".format(resDir, experiment, feature))
        else:
            plt.savefig("{}lessening{}/narrow_hist_{}.png".format(resDir, experiment, feature))
        
        plt.close()
