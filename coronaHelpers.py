# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 20:46:25 2016

Helper scripts for corona analysis

@author: vgor
"""

import re
import requests
import numpy as np
from pandas import DataFrame
from lxml import etree
from pyteomics.openms import featurexml
from pyteomics import mzid
from pyteomics.mass import fast_mass
from pyteomics.electrochem import pI
from gravy import gravy, AAContent
from collections import defaultdict

#proton
proton = 1.007276

#regular expressions
psmidRE = re.compile(r".+(?:\\)(\w+)_(SII_\d+_\d+)")

#functions
def createPSM(element):
    """
    Create dictionary representation of PSM element from Percolator XML output
    """
    defns = element.nsmap[None] #default namespace
    
    psm = {}
    psm["ID"] = element.get("{{{}}}psm_id".format(defns))
    psm["SVMScore"] = float(element.find("{{{}}}svm_score".format(defns)).text)
    psm["q-value"] = float(element.find("{{{}}}q_value".format(defns)).text)
    psm["PEP"] = float(element.find("{{{}}}pep".format(defns)).text)
    psm["ExpMass"] = float(element.find("{{{}}}exp_mass".format(defns)).text)
    psm["CalcMass"] = float(element.find("{{{}}}calc_mass".format(defns)).text)
    psm["Protein"] = element.find("{{{}}}protein_id".format(defns)).text
    psm["p-value"] = float(element.find("{{{}}}p_value".format(defns)).text)
    pepSeq = element.find("{{{}}}peptide_seq".format(defns))
    psm["Sequence"] = "{}.{}.{}".format(pepSeq.get("n"), pepSeq.get("seq"), pepSeq.get("c"))
    
    return psm

def createPeptide(element):
    """
    Create dictionary representation of Peptide element from Percolator XML output
    """
    defns = element.nsmap[None] #default namespace
    
    peptide = {}
    peptide["ID"] = element.get("{{{}}}peptide_id".format(defns))
    peptide["SVMScore"] = float(element.find("{{{}}}svm_score".format(defns)).text)
    peptide["q-value"] = float(element.find("{{{}}}q_value".format(defns)).text)
    peptide["PEP"] = float(element.find("{{{}}}pep".format(defns)).text)
    peptide["ExpMass"] = float(element.find("{{{}}}exp_mass".format(defns)).text)
    peptide["CalcMass"] = float(element.find("{{{}}}calc_mass".format(defns)).text)
    peptide["Protein"] = element.find("{{{}}}protein_id".format(defns)).text
    peptide["p-value"] = float(element.find("{{{}}}p_value".format(defns)).text)
    peptide["PSMs"] = element.xpath(".//x:psm_id/text()", namespaces = {"x": defns}, smart_strings=False)
    
    return peptide

def createProtein(element):
    """
    Create dictionary representation of Protein element from Percolator XML output
    """
    defns = element.nsmap[None] #default namespace
    
    protein = {}
    protein["ID"] = element.get("{{{}}}protein_id".format(defns))
    protein["q-value"] = float(element.find("{{{}}}q_value".format(defns)).text)
    protein["Peptides"] = element.xpath("x:peptide_seq/@seq", namespaces = {"x": defns}, smart_strings=False)
    
    return protein
    
def readPercOut(filename):
    """
    Read psms, peptides and proteins from percolator output XML
    """
    pout = etree.parse(filename)
    
    psms = DataFrame(list(map(createPSM, pout.iterfind("//{*}psm"))))
    peptides = DataFrame(list(map(createPeptide, pout.iterfind("//{*}peptide"))))
    proteins = DataFrame(list(map(createProtein, pout.iterfind("//{*}protein"))))
        
    return psms, peptides, proteins
    
def createFeature(element):
    """
    Trim dictionary representation of feature from featureXML output
    """
    feature = {}
    keys = ["FWHM", "charge", "intensity", "overallquality", "label",\
                "spectrum_index", "spectrum_native_id"]
    newkeys = ["FWHM", "Charge", "Intensity", "Overallquality", "Label",\
                "spectrum_index", "spectrum_native_id"]
    
    for key, newkey in zip(keys, newkeys):
        feature[newkey] = element[key]
        
    feature["RT"] = element["position"][0]["position"]
    feature["mz"] = element["position"][1]["position"]
    hullX = [point["x"] for point in element["convexhull"][0]["pt"]]
    feature["RTmin"] = min(hullX)
    feature["RTmax"] = max(hullX)
    
    return feature
    
    
def readFeatures(filename):
    """
    Read features from featureXML format
    """
    features = DataFrame(map(createFeature, featurexml.read(filename)))
    
    return features

def readMzId(filename):
    """
    Collect Spectrum Identification Information from mzIdentML
    """
    siis = []
    
    for sir in mzid.MzIdentML(filename):
       for sii in sir["SpectrumIdentificationItem"]:
           result = {}
           result["SII"] = sii["id"]
           result["Charge"] = sii["chargeState"]
           result["CalcMZ"] = sii["calculatedMassToCharge"]
           result["ExpMZ"] = sii["experimentalMassToCharge"]
           result["PepRef"] = sii["peptide_ref"]
           result["SIRID"] = sir["id"]
           result["Scan"] = sir["scan number(s)"]
           result["RT"] = sir["scan start time"]
           result["SpectrumID"] = sir["spectrumID"]
           result["q-value"] = sii.get("MS-GF:QValue", np.nan)
           siis.append(result)
               
    df = DataFrame(siis)    
    df.index = df["SII"]
    return df

def mapQValues(idFilename, psmTable, outFilename):
    """
    Map percolator q-values to mzid file
    """
    psmTable.index = psmTable["SII"]
    
    print("Reding mzIdentML from {}".format(idFilename))
    idFile = etree.parse(idFilename)

    ns = idFile.getroot().nsmap[None]

    notfound = 0
    found = 0
    for sii in idFile.xpath("//a:DataCollection//a:SpectrumIdentificationItem", namespaces = {"a": ns}):
        siiID = sii.get("id")
        try:
            qvalue = str(psmTable.loc[siiID, "q-value"])
            found += 1
        except KeyError:
            qvalue = "1.0"
            notfound += 1
        
        #remove previous q-values if any
        for elem in sii.xpath("a:cvParam[@accession='MS:1002354']", namespaces = {"a": ns}):
            sii.remove(elem)
            
        sii.append(etree.Element("cvParam", 
                                 {"cvRef": "PSI-MS",
                                  "accession": "MS:1002354",
                                  "name": "PSM-level q-value",
                                  "value": qvalue}
                                  ))
    
    print("PSMs in input: {}\nFound: {}\nNotfound: {}".format(len(psmTable), found, notfound))
    print("Writing mzIdentML to {}".format(outFilename))
    idFile.write(outFilename, xml_declaration = True, pretty_print = True, encoding = "UTF-8")

def splitIDs(psms):
    """
    Split filename and SII from ID in percolator format
    """
    psms[["File", "SII"]] = psms["ID"].str.extract(psmidRE, expand = True)


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

    return requests.post('https://www.uniprot.org/uploadlists/', params)

def cleanIDs(idString):
    """
    Convert search engine protein IDs to UniprotIDs
    Exclude contaminants
    """
    #the format of ID is like ??|UNIPROTID|???, ??|UNIPROTID|???, ....
    elements = [e.split("|") for e in idString.split(",")]
    elements = [z for z in elements if len(z) == 3]
    if len(elements) > 1: #more than one protein ID (happens when proteins are indistinguishable)
        accNrs = [s[1].replace("CONT_", "") for s in elements]
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
            if accNrs[0] == "P62807" or accNrs[1] == "P62807": #Manual fix for two histone isoforms in the database
                return "P62807"
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
        try:
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
            
        except Exception as e:
            print("{}, while processing entry {} on line {}".format(repr(e), entry.tag, entry.sourceline))
    
    data = DataFrame(data, columns = ["UniprotID", "GeneID", "Name", "GOComponent", "GOFunction", "GOProcess", "Sequence", "Disulfides", "BoundC"])
    
    #remove non-standard AA from sequences
    data["Sequence"] = data["Sequence"].str.replace(r"[XBJ]", "")
    
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

def parseIndex(datafile):
    """
    Parsing of AA index file from Kyoto database
    """
    parseRE = re.compile(r"^H\s+(\w+)[\w\W]+?^D\s+([\w\W]+?)$[\W\w]+?^(R\s+[\W\w]+?)^C[\W\w]+?^I\s+([\w\W]+?)\/\/", flags = re.M)
    #Will match one feature entry in Kyoto database and return match object with 4 groups:
    #1 - Feature Code (1 word)
    #2 - Feature description (1 line)
    #3 - references (multiline, as PMID, authors, title, journal while some of lines could be ommited or empty)
    #4 - amino acid indices (3 lines)
    with open(datafile, "r") as data:
        vectors = []
        names = []
        descriptions = []
        references = []
        for feature in re.finditer(parseRE, data.read()):
            header, line1, line2 = list(map(str.split, feature.group(4).splitlines()))
            vector = {}
            for i in range(len(header)):#header line is AA1/AA2 referencing line 1 and line 2
                aa1, aa2 = header[i].split("/")
                try:
                    vector[aa1] = float(line1[i])
                except Exception:#missing values
                    vector[aa1] = 0.0
            
                try:
                    vector[aa2] = float(line2[i])
                except Exception:
                    vector[aa2] = 0.0
            
            vectors.append(vector)
            names.append(feature.group(1))
            descriptions.append(feature.group(2))
            references.append(" ".join([s[2:] for s in feature.group(3).splitlines()]))
        
    return DataFrame(vectors, index = names), DataFrame({"Description": descriptions, "Reference": references}, index = names)
