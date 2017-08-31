# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 20:46:25 2016

Helper scripts for corona analysis

@author: vgor
"""

import re
import numpy as np
from pandas import DataFrame
from lxml import etree
from pyteomics.openms import featurexml
from pyteomics import mzid

#proton
proton = 1.007276

#regular expressions
psmidRE = re.compile(r".+(?:\\)(\w+)_(SII_\d+_\d+)")

#classes
class RTtransformer(object):
    """
    Retention time transformation class
    """
    def __init__(self, trafoxml):
        self.transformation = "Unknown"
        self.transformationTable = None
        self.length = 0
        
        try:#load trafoXML
            trafo = etree.parse(trafoxml)
        except Exception as e:
            print "Error openining {}\n{}".format(trafoxml, e.message)
            return
        
        self.transformation = trafo.xpath("string(//Transformation/@name)")
        
        if self.transformation != "identity":#identity does not coantain transformation table
            fromPoints = np.array(trafo.xpath("//Transformation/Pairs/Pair/@from"), dtype = float)
            toPoints = np.array(trafo.xpath("//Transformation/Pairs/Pair/@to"), dtype = float)
            self.transformationTable = DataFrame({"from": fromPoints, "to": toPoints})
            self.transformationTable.sort("from", inplace = True)
            self.transformationTable.reset_index(drop = True, inplace = True)
            
            self._a, self._b = np.polyfit(fromPoints, toPoints, 1) #linear fit
            
            self.length = len(fromPoints)
    
    def __call__(self, point):
        """
        Transform one point
        """
        if self.transformation == "identity": #identity transformation
            return point
        
        #use linear fit for values outside transformation table
        if point < self.transformationTable.loc[0, "from"] or \
            point > self.transformationTable.loc[self.length - 1, "from"]:
            return self._a * point + self._b
        #use linear approximation of two closest points in the table otherwise
        else:
            upper = np.argwhere(self.transformationTable["from"] >= point)[0, 0]
            if upper > 0:
                x1, x2 = self.transformationTable.loc[upper-1:upper, "from"]
                y1, y2 = self.transformationTable.loc[upper-1:upper, "to"]
            else:#special case if the point == transformationTable[0, "from"]
                return self.transformationTable.loc[upper, "to"]
            
            return y1 + (y2 - y1) / (x2 - x1) * (point - x1)
    
    def __repr__(self):
        return "Transformation\nType: {}\nTable:\n{}\n".format(self.transformation, self.transformationTable)

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
    peptide["PSMs"] = element.xpath(".//x:psm_id/text()", namespaces = {"x": defns})
    
    return peptide

def createProtein(element):
    """
    Create dictionary representation of Protein element from Percolator XML output
    """
    defns = element.nsmap[None] #default namespace
    
    protein = {}
    protein["ID"] = element.get("{{{}}}protein_id".format(defns))
    protein["q-value"] = float(element.find("{{{}}}q_value".format(defns)).text)
    protein["Peptides"] = element.xpath("x:peptide_seq/@seq", namespaces = {"x": defns})
    
    return protein
    
def readPercOut(filename):
    """
    Read psms, peptides and proteins from percolator output XML
    """
    pout = etree.parse(filename)
    
    psms = DataFrame(map(createPSM, pout.iterfind("//{*}psm")))
    peptides = DataFrame(map(createPeptide, pout.iterfind("//{*}peptide")))
    proteins = DataFrame(map(createProtein, pout.iterfind("//{*}protein")))
    
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
    
    print "Reding mzIdentML from {}".format(idFilename)
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
    print "Writing mzIdentML to {}".format(outFilename)
    idFile.write(outFilename, xml_declaration = True, pretty_print = True, encoding = "UTF-8")

def splitIDs(psms):
    """
    Split filename and SII from ID in percolator format
    """
    psms["File"] = psms["ID"].apply(lambda s: psmidRE.match(s).groups()[0])
    psms["SII"] = psms["ID"].apply(lambda s: psmidRE.match(s).groups()[1])
    
