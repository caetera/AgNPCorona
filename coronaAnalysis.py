# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:42:16 2017

Data processing and figure generation

@author: vgor
"""
from __future__ import print_function
import re
import numpy as np
import pylab as plt
import requests
from coronaHelpers import readPercOut
from coronaPlots import *
from lxml import etree
from pandas import DataFrame, Series, read_csv, concat
from collections import defaultdict
from pyteomics.pylab_aux import scatter_trend
from pyteomics.mass import fast_mass
from pyteomics.electrochem import pI
from gravy import gravy, AAContent
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind as ttest

#Mapping UNIMOD ID to trivial name
unimodMap = {"1": "Acetyl",
             "4": "Carbamidomethyl",
             "5": "Carbamyl",
             "35": "Oxidation"}


def quantifyProteins(name):
    """
    TopN Quantition using peptides from OpenMS
    name - the name of the file to work with
    """
    def quantifyProtein(peptidelist, aggFunc = np.nanmedian, N = 3):
        """
        Perform the quantification for individual protein
        aggFunc - function used for aggreagtion (ex. mean)
        N as in TopN
        """
        data = iqpeptides.ix[peptidelist, :]
    
        result = Series()
    
        for i in "123":
            topN = data.sort("abundance_" + i, ascending = False)[:N]
            if np.all(np.isfinite(topN["abundance_" + i ])):
                result["abundance" + i] = aggFunc(topN["abundance_" + i ])
    
        return result

    #load percolator results
    psms, peptides, proteins = readPercOut("E:\\RawData\\20161120_NP\\{}.pout.xml".format(name))
    
    #load OpenMS peptides results
    qpeptides = read_csv("E:\\RawData\\20161120_NP\\{}.peptides.csv".format(name), sep = "\t", skiprows = 3)
    
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
    iproteins[["{}_{}".format(name, i) for i in "123"]] = iproteins["Peptides"].apply(quantifyProtein)
    
    #find non quantified proteins (all Abundances are NaN)
    qind = np.any(np.isfinite(iproteins[["{}_{}".format(name, i) for i in "123"]].values), axis = 1)
    iproteins.drop(["Peptides", "q-value"], axis = 1, inplace = True)
    iproteins = iproteins[qind]
    
    #report information about identifications and quantifications
    #identified peptides, quantified peptides, identified proteins, quantified proteins
    print(name, sum(peptides["q-value"] < 0.01), len(iqpeptides), sum(proteins["q-value"] < 0.01), len(iproteins))
    
    #make interreplicate correlation plots
    plt.figure(figsize = (7, 9))
    plt.subplot(3, 1, 1)
    ind = np.logical_and(iproteins["{}_1".format(name)] > 0, iproteins["{}_2".format(name)]> 0)
    scatter_trend(np.log10(iproteins.loc[ind, "{}_1".format(name)]), np.log10(iproteins.loc[ind, "{}_2".format(name)]))
    plt.xlabel("log10({}_1 abundance)".format(name))
    plt.ylabel("log10({}_2 abundance)".format(name))
    plt.subplot(3, 1, 2)
    ind = np.logical_and(iproteins["{}_1".format(name)] > 0, iproteins["{}_3".format(name)]> 0)
    scatter_trend(np.log10(iproteins.loc[ind, "{}_1".format(name)]), np.log10(iproteins.loc[ind, "{}_3".format(name)]))
    plt.xlabel("log10({}_1 abundance)".format(name))
    plt.ylabel("log10({}_3 abundance)".format(name))
    plt.subplot(3, 1, 3)
    ind = np.logical_and(iproteins["{}_2".format(name)] > 0, iproteins["{}_3".format(name)]> 0)
    scatter_trend(np.log10(iproteins.loc[ind, "{}_2".format(name)]), np.log10(iproteins.loc[ind, "{}_3".format(name)]))
    plt.xlabel("log10({}_2 abundance)".format(name))
    plt.ylabel("log10({}_3 abundance)".format(name))    
    
    plt.savefig("E:\\RawData\\20161120_NP\\{}_reps.png".format(name))
    plt.close()
    
    return iproteins

def processExperimentT(name):
    """
    Do the qunatification on temperature experiments
    """
    #collect data from all replicates and contol
    d = quantifyProteins("{}1".format(name))
    d = d.merge(quantifyProteins("{}2".format(name)), on = "ID", how = "outer")
    d = d.merge(quantifyProteins("{}3".format(name)), on = "ID", how = "outer")
    d = d.merge(quantifyProteins("{}C".format(name)), on = "ID", how = "outer")
    
    #names of columns with numerical data
    datacolumns = d.columns[1:]
    
    #log10 transformation, cleaning, saving the result into plain csv
    #protein ID (as reported by search engine vs quantification results)
    d[datacolumns] = d[datacolumns].apply(np.log10, axis = 1)
    d.replace(-np.inf, 0, inplace = True)
    d.set_index("ID", verify_integrity = True, inplace = True)
    d.to_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))      
    
def processExperimentP(name):
    """
    Do the qunatification on pH experiments
    """
    #collect data from all replicates and control
    d = quantifyProteins("{}_B1".format(name))
    d = d.merge(quantifyProteins("{}_B2".format(name)), on = "ID", how = "outer")
    d = d.merge(quantifyProteins("{}_B3".format(name)), on = "ID", how = "outer")
    d = d.merge(quantifyProteins("{}_C".format(name)), on = "ID", how = "outer")
    
    #names of the columns containing numerical data
    datacolumns = d.columns[1:]
    
    #log10 transformation, cleaning, saving the result into plain csv
    #protein ID (as reported by search engine vs quantification results)
    d[datacolumns] = d[datacolumns].apply(np.log10, axis = 1)
    d.replace(-np.inf, 0, inplace = True)
    d.set_index("ID", verify_integrity = True, inplace = True)
    d.to_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))
    
def splitID(s):
    """
    Extract Uniprot ID(s) from protein name
    """
    ids = [x.split("|")[1] for x in s.split(",")]
    
    return " ".join(ids)

def queryUniprot(accessionList, resFormat = "xml"):
    """
    Retrieve list of accessions from Uniprot in specific format
    resFormat - the resulting format as in Uniprot API
    """
    params = {
        'from': 'ACC',
        'to': 'ACC',
        'format': resFormat,
        'query':' '.join(accessionList)
        }

    return requests.post('http://www.uniprot.org/uploadlists/', params)

def cleanIDs(idString):
    """
    Convert search engine protein IDs to UniprotIDs
    Exclude contaminants
    """
    #the format of ID is like ??|UNIPROTID|???, ??|UNIPROTID|???, ....
    elements = [e.split("|") for e in idString.split(",")]
    elements = filter(lambda z: len(z) == 3, elements)
    if len(elements) > 1: #more than one protein ID (happens when proteins are indistinguishable)
        accNrs = map(lambda s: s[1].replace("CONT_", ""), elements)
        if all([accNrs[0] == a for a in accNrs]): #all IDs are the same
            if accNrs[0] == "P02768": #map BSA to HSA
                return accNrs[0]
            else:
                #Internal check
                print(idString, accNrs)
                return ""
        else:
            if accNrs[0] == "P60712": #Asign as actin
                return accNrs[1]
            else:
                #Internal check
                print(idString, accNrs)
                return ""
                    
    else:
        if elements[0][1].startswith("CONT_"): #contaminants
            return ""
        else:
            return elements[0][1] #regular protein
            
def getDisulfidePositions(entry, ns):
    """
    Parse positions of cysteines bound by disulfide bonds from Uniprot XML
    """
    positions = map(lambda ee: ee.get("position"), entry.findall("./n:feature[@type=\"disulfide bond\"]/n:location/*", namespaces = ns))
    
    positions = filter(lambda ee: not (ee is None), positions)
    
    return set([int(p) for p in positions])
    
def getProteinData(ids):
    """
    Retrieve protein data from Uniprot and parse it in to ready to use format
    """
    #create namespace dictionary for further use
    ns = {"n": "http://uniprot.org/uniprot"}
    #retrieve the data
    print("Retriving {} entries from Uniprot".format(len(ids)))
    xmlData = etree.fromstring(queryUniprot(ids).content)
    print("Uniprot done. Parsing")
    #read XML and extract basic data
    data = []
    for entry in xmlData.findall("n:entry", namespaces = ns):
        accNr = entry.find("./n:accession", namespaces = ns).text
        geneID = entry.find("./n:gene/n:name[@type=\"primary\"]", namespaces = ns).text
        fullname = entry.find(".//n:fullName", namespaces = ns).text
        sequence = entry.find("./n:sequence", namespaces = ns).text
        sequence = sequence.replace("\n", "").strip()
        nDisulfide = len(entry.findall("./n:feature[@type=\"disulfide bond\"]", namespaces = ns))
        bC = getDisulfidePositions(entry, ns)
        goTerms = entry.xpath("./n:dbReference[@type=\"GO\"]/n:property[@type=\"term\"]/@value", namespaces = ns)
        goData = defaultdict(list)
        for term in goTerms:
            goData[term[0]].append(term[2:])
        
        data.append((accNr, geneID, fullname, goData["C"], goData["F"], goData["P"], sequence, nDisulfide, bC))
    
    data = DataFrame(data, columns = ["UniprotID", "GeneID", "Name", "GOComponent", "GOFunction", "GOProcess", "Sequence", "Disulfides", "BoundC"])
    
    #calculate some additional properties
    data["MW"] = data["Sequence"].apply(fast_mass)
    data["pI"] = data["Sequence"].apply(pI)
    data["GRAVY"] = data["Sequence"].apply(gravy)
    data["AAC"] = data["Sequence"].apply(AAContent, frequencies = False)
    
    #check if all IDs are assigned
    idCheck = set(data["UniprotID"]) == set(ids)
    print("ID check: {}".format(idCheck))
    if not idCheck:
        print("IDs missing in Uniprot: {}".format(set(data["UniprotID"]).difference(set(ids))))
        print("IDs missing in query: {}".format(set(ids).difference(set(data["UniprotID"]))))
    
    return data

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
        dd = dd.union(set(resT.sort(label, ascending = False).index[:N]))
    
    return resT.loc[list(dd), :].fillna(0)

def sigmoid(x, x0, k):
    """
    Sigmoid function definition
    """
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return y
     
def get_sigmoid_fit(row):
    """
    Perform fittinig with the sigmoid function on one row (protein amount in 5 conditions)
    x - is defined externally and should be scaled to [0, 1]
    Returns float number, with absolute value corresponding to the intercept (in scaled coordiantes)
    and sign corresponding to the direction of change
    """
    y = row[:5].values.astype(float)
    y = (y - min(y)) / (max(y) - min(y)) #scale y to [0, 1]
    try:
        opt, cov = curve_fit(sigmoid, x, y, [0.5, 0.5])
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
        opt, cov = curve_fit(sigmoid, x_scaled, y_scaled, [0.5, 0.5])
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
for name in ["Temp" + x for x in "ABCDE"]:
    processExperimentT(name)
    
for name in ["pH" + x for x in ["3", "5", "6_6", "7", "9"]]:
    processExperimentP(name)

#%%Compute average values for control samples
#pH experiment
#read quantification results and merge them together
data = []
for name in ["pH" + x for x in ["3", "5", "6_6", "7", "9"]]:
    d = read_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))
    #keep only columns that correspond to control
    d = d.iloc[np.all(d.iloc[:, -4:-1] != 0, axis = 1), [0, -3, -2, -1]]
    d.set_index("ID", drop = True, inplace = True)
    data.append(d)

md = concat(data, join = "outer", axis = 1)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#Summarize the data at the biologiacal repicate level (i.e. join technical repeats)
res = DataFrame(index = md.index)
for c in ["3", "5", "6_6", "7", "9"]:
    res["pH{}".format(c)] = \
        md.loc[:, ["pH{}_C_{}".format(c, i) for i in "123"]]\
        .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res.index)

#retrieve the properties and map them to proteins
res["CleanID"] = map(cleanIDs, res.index)
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#calculate average values for control samples
#the plot routine returns average values
colnames = ["pH{}".format(c) for c in ["3", "5", "6_6", "7", "9"]]
pImeansP = pIPlot(mdc, colnames, colors = colors)
MWmeansP = MWPlot(mdc, colnames, colors = colors)
GravyMeansP = GravyPlot(mdc, colnames, colors = colors)

#Temperature experiments
#read quantification and merge them together
data = []
for name in ["Temp" + x for x in "ABCDE"]:
    d = read_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))
    #keep the columns corresponding to control
    d = d.iloc[np.all(d.iloc[:, -4:-1] != 0, axis = 1), [0, -3, -2, -1]]
    d.set_index("ID", drop = True, inplace = True)
    data.append(d)

md = concat(data, join = "outer", axis = 1)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#Summarize the data at the biologiacal repicate level (i.e. join technical repeats)
res = DataFrame(index = md.index)
for n in ["Temp{}".format(c) for c in "ABCDE"]:
        res[n] = \
            md.loc[:, ["{}C_{}".format(n, i) for i in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res.index)

#retrieve the properties and map them to proteins
res["CleanID"] = map(cleanIDs, res.index)
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#calculate average values for control samples
#the plot routine returns average values
colnames = ["Temp{}".format(c) for c in "ABCDE"]
pImeansT = pIPlot(mdc, colnames, colors = colors)
MWmeansT = MWPlot(mdc, colnames, colors = colors)
GravyMeansT = GravyPlot(mdc, colnames, colors = colors)

#%%Perform the analysis with control subtraction
#pH experiments
#collect the quan data and merging
data = []
for name in ["pH" + x for x in ["3", "5", "6_6", "7", "9"]]:
    d = read_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))
    #keep the columns present in sample or control
    d = d.iloc[np.any(d.iloc[:, 1:10] != 0, axis = 1), :]
    d.set_index("ID", drop = True, inplace = True)
    data.append(d)

md = concat(data, join = "outer", axis = 1)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#join technical replicates
res = DataFrame(index = md.index)
for c in ["3", "5", "6_6", "7", "9"]:
    for r in ["B1", "B2", "B3", "C"]:
        res["pH{}_{}".format(c, r)] = \
            md.loc[:, ["pH{}_{}_{}".format(c, r, i) for i in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#subtract protein levels in control samples
res.fillna(0, inplace = True)
for c in ["3", "5", "6_6", "7", "9"]:
    for r in ["B1", "B2", "B3"]:
        res["pH{}_{}".format(c, r)] = res["pH{}_{}".format(c, r)] - res["pH{}_C".format(c)]
    
    res.drop("pH{}_C".format(c), axis = 1, inplace = True)

#remove negative values (after subtraction)
res[res <= 0] = np.nan
res = res[res.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#extract the data for heatmap (built by external tool)
zz = res.copy()
zz.columns = ["{} - {}".format(a,b) for a in ["pH4.9", "pH6.1", "pH6.8", "pH7.7", "pH8.9"] for b in "123"]
zz.to_csv("E:\\RawData\\20161120_NP\\res\\data.csv")

#make PCA plot 
plt.figure(figsize = (9,6))
PCAplot(res, ["4.9", "6.1", "6.8", "7.7", "8.9"])
plt.savefig("E:\\RawData\\20161120_NP\\res\\PCA_pH.svg", bbox_inches = "tight")

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res.index)

#retrieve the properties and map them to proteins
res["CleanID"] = map(cleanIDs, res.index)
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#create distribution plots
colnames = ["pH{}_B{}".format(c, r) for c in ["3", "5", "6_6", "7", "9"] for r in "123"]
colors = ["#cf7b43", "#7d6acb", "#75b140", "#c459b6", "#54a676", "#cc4b42", "#6797d0", "#a5903e", "#c45d83"]

plt.figure(figsize = (9,6))
pIPlot(mdc, colnames, colors = colors, means = pImeansP)
plt.savefig("E:\\RawData\\20161120_NP\\res\\pI_pH.svg", bbox_inches = "tight")
plt.figure(figsize = (9,6))
MWPlot(mdc, colnames, colors = colors, means = MWmeansP)
plt.savefig("E:\\RawData\\20161120_NP\\res\\MW_pH.svg", bbox_inches = "tight")
plt.figure(figsize = (9,6))
GravyPlot(mdc, colnames, colors = colors, means = GravyMeansP)
plt.savefig("E:\\RawData\\20161120_NP\\res\\GRAVY_pH.svg", bbox_inches = "tight")

#Summarize the data at the condition level (i.e. join biological replicates)
res2 = DataFrame(index = res.index)
for c in ["3", "5", "6_6", "7", "9"]:
    res2["pH{}".format(c)] = \
            res.loc[:, ["pH{}_B{}".format(c, r) for r in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res2.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res2.index)

#retrieve the properties and map them to proteins
res2["CleanID"] = map(cleanIDs, res2.index)
mdc2 = res2[res2["CleanID"] != ""].merge(getProteinData(res2.loc[res2["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc2.drop(["CleanID", "Sequence"], axis = 1, inplace = True)
#save the results for later use
mdc2.to_csv("E:\\RawData\\20161120_NP\\res\\pH_conditionDetail.csv", index = False)

#Temperature experiments
#collect the quan data and merging
data = []
for name in ["Temp" + x for x in "ABCDE"]:
    d = read_csv("E:\\RawData\\20161120_NP\\{}_raw.csv".format(name))
    #keep the columns present in sample or control
    d = d.iloc[np.any(d.iloc[:, 1:10] != 0, axis = 1), :]
    d.set_index("ID", drop = True, inplace = True)
    data.append(d)

md = concat(data, join = "outer", axis = 1)

#remove empty lines
md = md[md.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#join technical replicates
res = DataFrame(index = md.index)
for n in ["Temp{}{}".format(c,r) for c in "ABCDE" for r in "123C"]:
        res[n] = \
            md.loc[:, ["{}_{}".format(n, i) for i in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#subtract protein levels in control samples
res.fillna(0, inplace = True)
for c in "ABCDE":
    for r in "123":
        res["Temp{}{}".format(c, r)] = res["Temp{}{}".format(c, r)] - res["Temp{}C".format(c)]
    
    res.drop("Temp{}C".format(c), axis = 1, inplace = True)

res[res <= 0] = np.nan
res = res[res.apply(lambda x: np.any(np.isfinite(x)), axis = 1)]

#save the data for heatmap construction (external tool)
zz = res.copy()
zz.columns = ["{} - {}".format(a,b) for a in ["4C", "17C", "30C", "41C", "47C"] for b in "123"]
zz.to_csv("E:\\RawData\\20161120_NP\\res\\data.csv")

#PCA plot
plt.figure(figsize = (9,6))
PCAplot(res, ["4", "17", "30", "41", "47"])
plt.savefig("E:\\RawData\\20161120_NP\\res\\PCA_Temp.svg", bbox_inches = "tight")

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res.index)

#retrieve the properties and map them to proteins
res["CleanID"] = map(cleanIDs, res.index)
mdc = res[res["CleanID"] != ""].merge(getProteinData(res.loc[res["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc.drop(["CleanID", "Sequence"], axis = 1, inplace = True)

#build distribution plots
colnames = ["Temp{}{}".format(c, r) for c in "ABCDE" for r in "123"]
colors = ["#cf7b43", "#7d6acb", "#75b140", "#c459b6", "#54a676", "#cc4b42", "#6797d0", "#a5903e", "#c45d83"]

plt.figure(figsize = (9,6))
pIPlot(mdc, colnames, colors = colors, means = pImeansT)
plt.savefig("E:\\RawData\\20161120_NP\\res\\pI_Temp.svg", bbox_inches = "tight")
plt.figure(figsize = (9,6))
MWPlot(mdc, colnames, colors = colors, means = MWmeansT)
plt.savefig("E:\\RawData\\20161120_NP\\res\\MW_Temp.svg", bbox_inches = "tight")
plt.figure(figsize = (9,6))
GravyPlot(mdc, colnames, colors = colors, means = GravyMeansT)
plt.savefig("E:\\RawData\\20161120_NP\\res\\GRAVY_Temp.svg", bbox_inches = "tight")

#Join biological replicates
res2 = DataFrame(index = res.index)
for c in "ABCDE":
    res2["Temp{}".format(c)] = \
            res.loc[:, ["Temp{}{}".format(c, r) for r in "123"]]\
            .apply(lambda x: np.nansum(x)/sum(np.logical_and(np.isfinite(x), x != 0)), axis = 1)

#Correct the accesion number of Calmodulin (was changed in Uniprot after the search)
res2.index = map(lambda s: 'sp|P0DP23|CALM_HUMAN' if s == 'sp|P62158|CALM_HUMAN'  else s, res2.index)

#retrieve the properties and map them to proteins
res2["CleanID"] = map(cleanIDs, res2.index)
mdc2 = res2[res2["CleanID"] != ""].merge(getProteinData(res2.loc[res2["CleanID"] != "", "CleanID"]),\
                                         left_on = "CleanID", right_on = "UniprotID", how = "inner")
                                         
mdc2.drop(["CleanID", "Sequence"], axis = 1, inplace = True)
#save data for later
mdc2.to_csv("E:\\RawData\\20161120_NP\\res\\Temp_conditionDetail.csv", index = False)


#%%Identify differentially abundant proteins
#Performing fit with sigmoid curve
#Real temperature values
x = np.array([4., 17., 30., 41., 47.])
#Real pH values
x = np.array([4.9, 6.1, 6.8, 7.7, 8.9])
#scale x before perfoming fit
x = (x - x.min())/ (x.max() - x.min())

#fit itself
mdc2["SigmaFit"] = mdc2[mdc2.iloc[:, :5].apply(lambda x: np.all(np.isfinite(x)), axis = 1)].apply(get_sigmoid_fit, axis = 1)
mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "Direction"] = mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "SigmaFit"].apply(lambda t: "Up" if t > 0 else "Down")
#pH
mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "CritpH"] = mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "SigmaFit"].apply(lambda t: 4.9 + abs(t) * 4)
mdc2[np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)].to_csv("E:\\RawData\\20161120_NP\\res\\omnipHS.csv", index = False)
mdc2[np.isfinite(mdc2["SigmaFit"])].to_csv("E:\\RawData\\20161120_NP\\res\\sigmoidpHS.csv", index = False)
#temperature
mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "CritTemp"] = mdc2.loc[np.isfinite(mdc2["SigmaFit"]), "SigmaFit"].apply(lambda t: 4.0 + abs(t) * 43)
mdc2[np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)].to_csv("E:\\RawData\\20161120_NP\\res\\omniTempS.csv", index = False)
mdc2[np.isfinite(mdc2["SigmaFit"])].to_csv("E:\\RawData\\20161120_NP\\res\\sigmoidTempS.csv", index = False)

#plotting sigmoid fits
#pH
x = np.array([4.9, 6.1, 6.8, 7.7, 8.9])
cInd = dict(zip(mdc2.columns[:5], [map(lambda t: t.startswith(s), mdc.columns) for s in mdc2.columns[:5]]))
for i, row in mdc2.loc[np.isfinite(mdc2["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc2.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "pH", "log10(Abundance)")
    plt.savefig("E:\\RawData\\20161120_NP\\res\\SigmaFitspHS\\{}.svg".format(row["GeneID"]))
    plt.close()

#Flat fits
#select proteins quantified in all conditions
mdc3 = mdc2[np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)]

x = np.array([4.9, 6.1, 6.8, 7.7, 8.9])
cInd = dict(zip(mdc3.columns[:5], [map(lambda t: t.startswith(s), mdc.columns) for s in mdc3.columns[:5]]))
for i, row in mdc3.loc[~np.isfinite(mdc3["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc3.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "pH", "log10(Abundance)")
    plt.savefig("E:\\RawData\\20161120_NP\\res\\SigmaFitspHS-\\{}.svg".format(row["GeneID"]))
    plt.close()


#Temeperature
x = np.array([4., 17., 30., 41., 47.])
cInd = dict(zip(mdc2.columns[:5], [map(lambda t: t.startswith(s), mdc.columns) for s in mdc2.columns[:5]]))
for i, row in mdc2.loc[np.isfinite(mdc2["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc2.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "Temperature", "log10(Abundance)")
    plt.savefig("E:\\RawData\\20161120_NP\\res\\SigmaFitsTempS\\{}.svg".format(row["GeneID"]))
    plt.close()

#Flat fits
#select proteins quantified in all conditions
mdc3 = mdc2[np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)]

x = np.array([4., 17., 30., 41., 47.])
cInd = dict(zip(mdc3.columns[:5], [map(lambda t: t.startswith(s), mdc.columns) for s in mdc3.columns[:5]]))
for i, row in mdc3.loc[~np.isfinite(mdc3["SigmaFit"]), :].iterrows():
    err = [np.nanstd(mdc.loc[mdc["UniprotID"] == row["UniprotID"], cInd[t]]) for t in mdc3.columns[:5]]
    plt.figure()
    plt.title("{} ({})\n".format(row["Name"], row["GeneID"]))
    plot_sigmoid_fit(x, row[:5].values.astype(float), err, "Temperature", "log10(Abundance)")
    plt.savefig("E:\\RawData\\20161120_NP\\res\\SigmaFitsTempS-\\{}.svg".format(row["GeneID"]))
    plt.close()

#%%Top20
mdc2 = read_csv("E:\\RawData\\20161120_NP\\res\\pH_conditionDetail.csv")
getTopN(mdc2, ["pH3", "pH5", "pH6_6", "pH7", "pH9"], 20).to_clipboard(index  = False)

mdc2 = read_csv("E:\\RawData\\20161120_NP\\res\\Temp_conditionDetail.csv")
getTopN(mdc2, ["TempA", "TempB", "TempC", "TempD", "TempE"], 20).to_clipboard(index  = False)

#%%TPP and pI correlation
#Load ThermoProteomeProfiling data
TS4 = read_csv("E:\\RawData\\20161120_NP\\res\\TableS4.csv")[["Protein NAME", "meltP_Vehicle_Expt2"]]
#load critcal temperature
sigmT = read_csv("E:\\RawData\\20161120_NP\\res\\sigmoidTempS.csv")

#idenitidy overlapping proteins
overlap = TS4[np.isfinite(TS4["meltP_Vehicle_Expt2"])].merge(sigmT[np.isfinite(sigmT["SigmaFit"])], left_on = "Protein NAME", right_on = "GeneID")

#plot correlation
plt.figure(figsize = (8,5))
scatter_trend(overlap["CritTemp"], overlap["meltP_Vehicle_Expt2"],)
plt.xlabel("Critical Temperature")
plt.ylabel("Melting Temperature")
plt.savefig("E:\\RawData\\20161120_NP\\res\\meltPTemp.svg")


#pI to critical pH correlation
mdc2 = read_csv("E:\\RawData\\20161120_NP\\res\\sigmoidpHS.csv")
plt.figure(figsize = (8,5))
selector = np.isfinite(mdc2["SigmaFit"])
scatter_trend(mdc2.loc[selector, "CritpH"], mdc2.loc[selector, "pI"])
plt.xlabel("Critical pH")
plt.ylabel("pI")
plt.savefig("E:\\RawData\\20161120_NP\\res\\pIpH.svg")

#%%Other 
#Identification overlaps
def overlapBarplot(labels, protFormat):
    points = []
    errors = []
    names = []
    for a in range(len(labels) - 1):
        labelA = labels[a]
        labelB = labels[a+1]
        differenceA = [len(proteins[protFormat.format(labelA, i)].difference(proteins[protFormat.format(labelA, j)]))/\
                        float(len(proteins[protFormat.format(labelA, i)]))\
                        for (i,j) in [(1,2), (2,3), (1,3)]]        
        differenceB = [len(proteins[protFormat.format(labelA, i)].difference(proteins[protFormat.format(labelB, j)]))/\
                        float(len(proteins[protFormat.format(labelA, i)]))\
                        for i in (1,2,3) for j in (1,2,3)]
        
        names.extend(["{0} vs {0}".format(labelA), "{} vs {}".format(labelA, labelB)])
        points.extend([np.mean(differenceA), np.mean(differenceB)])
        errors.extend([np.std(differenceA), np.std(differenceB)])
        
        print(ttest(differenceA, differenceB, equal_var = False))
    
    plt.figure(figsize = (9, 7))
    plt.bar(range(len(names)), points, align = "center", alpha = 0.7)
    plt.errorbar(range(len(names)), points, yerr = errors, ecolor = "k", elinewidth = 2, fmt = None)
    plt.xticks(range(len(names)), names, rotation = "vertical")
    plt.ylabel("Difference fraction")
    plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on')

labels = "ABCDE"
proteins = {}
for e in labels:
    for r in "123":
        name = "Temp{}{}".format(e,r)
        proteins[name] = set(res[np.isfinite(res[name])].index)
        
overlapBarplot(labels, "Temp{}{}")
plt.savefig("E:\\RawData\\20161120_NP\\res\\IDoverlapTemp.svg")


labels = ["3", "5", "6_6", "7", "9"]
proteins = {}
for e in labels:
    for r in ["B1", "B2", "B3"]:
        name = "pH{}_{}".format(e, r)
        proteins[name] = set(res[np.isfinite(res[name])].index)

overlapBarplot(labels, "pH{}_B{}")
plt.savefig("E:\\RawData\\20161120_NP\\res\\IDoverlappH.svg")

#disulfides
testSet = mdc2[np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)]
testSet["TotalC"] = testSet["AAC"].apply(lambda x: x["C"])
testSet["UnboundC"] = testSet["TotalC"] - testSet["BoundC"].apply(len)
testSet["BoundC"] = testSet["BoundC"].apply(len)

testSet["Direction"] = testSet["Direction"].astype(str)
plt.boxplot([z[1].values for z in testSet.groupby("Direction")["UnboundC"]])
print([z[0] for z in testSet.groupby("Direction")])


#%%SVM test
from sklearn.svm import SVC
clf = SVC(kernel = 'linear')
X = DataFrame(list(mdc2.loc[(np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)), "AAC"]))
X.drop("U", inplace = True, axis = 1)
X.fillna(0, inplace = True)
Y = mdc2.loc[(np.all(np.isfinite(mdc2.loc[:, "TempA":"TempE"]), axis = 1)), "SigmaFit"].apply(np.isnan).values.astype(int)

clf.fit(X.values, Y)
print(clf.score(X.values, Y))
print(DataFrame(clf.coef_ * 100, columns = X.columns))

X = DataFrame(list(mdc2.loc[(np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)), "AAC"]))
X.drop("U", inplace = True, axis = 1)
X.fillna(0, inplace = True)
Y = mdc2.loc[(np.all(np.isfinite(mdc2.loc[:, "pH3":"pH9"]), axis = 1)), "SigmaFit"].apply(np.isnan).values.astype(int)

clf.fit(X.values, Y)
print(clf.score(X.values, Y))
print(DataFrame(clf.coef_ * 100, columns = X.columns))

from SupervisedPCA import SupervisedPCAClassifier

spca = SupervisedPCAClassifier()
spca.fit(X.values, Y)
pca_data = spca.get_transformed_data(X)
plt.scatter(pca_data[Y == 1, 0], pca_data[Y == 1, 1], color = 'b')
plt.scatter(pca_data[Y == 0, 0], pca_data[Y == 0, 1], color = 'r')

from sklearn.decomposition import PCA
pC = PCA(2)
pca_data = pC.fit_transform(X)
plt.scatter(pca_data[Y == 1, 0], pca_data[Y == 1, 1], color = 'b')
plt.scatter(pca_data[Y == 0, 0], pca_data[Y == 0, 1], color = 'r')

#%%Poster plots
#Identified and quantified proteins
z = read_csv("E:\\RawData\\20161120_NP\\res\\data.csv", sep = "\t")

plt.figure(figsize = (8,5))
zz = z[z["Exp"] == "P"]
i = 0
label = []
for c, g in zz.groupby("Condition"):
    plt.bar(range(i, i + len(g)), g["IDs"], align = "center", color = "teal", alpha = 0.6, label = "Identified")
    plt.bar(range(i, i + len(g)), g["Quant"], align = "center", color = "olive", alpha = 0.6, label = "Quantified")
    i += len(g) + 1
    label += ["1", "2", "3", "C", ""]

plt.tick_params(axis='both', which='both', left='off', right='off' , bottom='off', top='off', labelleft='on')
plt.xlim(-1, i - 1)
plt.legend()
plt.ylabel("Numeber of proteins")
plt.xticks(range(0,i), label)
plt.savefig("E:\\RawData\\20161120_NP\\res\\pHOverview.svg")
