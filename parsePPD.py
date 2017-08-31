# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 13:28:30 2016

Compile Plasma Proteome Database into FASTA format

@author: vgor
"""

from os import walk, stat
from lxml import etree
from pandas import DataFrame
from io import StringIO
from pyteomics import fasta
import re
import requests
import numpy as np

xmlPath = "E:\\Fasta\\PlasmaProteomeDB\\Experimental_evidence\\"
fastaPath = "E:\\Fasta\\"

#ill-formed XML tag
invalidXMLTag = re.compile('(<\/*)(\w+(?:\s\w+(?!\s*=\s*".*?"))+)((?:\s\w+\s*=\s*".*?")*)(\/*>)')

def curePPDSyntax(filename):
    """
    Remove spaces in XML tags that come from PPD
    i.e. <Experimental evidence> become <ExperimentalEvidence>
    """
    with open(filename, "r") as fileIn:
        text = fileIn.read()
        curedText = ""
        pos = 0
        for match in invalidXMLTag.finditer(text):
            parts = match.groups()
            #groups are
            #0 - opening brace 
            #1 - element name with spaces
            #2 - attributes
            #3 - closing brace
            curedTag = parts[0] + "".join(map(str.capitalize, parts[1].split())) + parts[2] + parts[3]
            curedText += text[pos:match.start()] + curedTag
            pos = match.end()
            
        curedText += text[pos:]
        
    return curedText

def queryUniprot(accessionList):
    """
    Retrieve list of accessions from Uniprot in FASTA format
    """
    params = {
        'from':'ACC',
        'to':'ACC',
        'format':'fasta',
        'query':' '.join(accessionList)
        }

    return requests.post('http://www.uniprot.org/uploadlists/', params).text
    

data = []

for name in walk(xmlPath, False).next()[2]:
    if stat(xmlPath + name).st_size > 0: #ignore empty files
        d = etree.fromstring(curePPDSyntax(xmlPath + name))
        
        #collect protein information
        protein = {}
        protein["filename"] = name
        protein["ID"] = d.xpath("//protein")[0].get("ID")
        protein["Description"] = d.xpath("//protein/title")[0].text.strip()
        protein["GeneID"] = d.xpath("//protein/gene_symbol")[0].text.strip()
        protein["UniprotID"] = d.xpath("//protein/swissprot")[0].text.strip()
        protein["NrReference"] = len(d.xpath("//Experimental_evidence[@Reference]"))
        protein["NrPlasma"] = len(d.xpath("//Experimental_evidence[@Sample=\"Plasma\"]"))
    
        data.append(protein)

data = DataFrame(data)#convert to DataFrame for easier manin=pulation

#keep only entries with more than 1 experimental evidence in total and at least one in plasma
ids = data.loc[np.logical_and(data["NrPlasma"] > 0, data["NrReference"] > 1), "UniprotID"].values
print "Number of entries from PPD:", len(ids)

#remove entries whitout SwissProt identifier
ids = ids[np.logical_and(ids != "-", ids != "None")]
print "Removed empty SwissProt identifiers:", len(ids)

#Has only a single SwissProt id
singleIDs = filter(lambda x: len(x) == 6, ids)
print "Have single ID:", len(singleIDs)

#some ids are duplicates, thus we will use set
singleIDfasta = queryUniprot(set(singleIDs))

#collect data from fasta
proteins = [{"Description": p.description, "Sequence": p.sequence} 
                for p in fasta.read(StringIO(singleIDfasta))]


multipleIDs = filter(lambda x: len(x) > 6, ids)
print "Have multiple IDs:", len(multipleIDs)

multipleIDfasta = ""

for ID in multipleIDs:
    fastaText = queryUniprot(ID.split(","))
    multipleIDfasta += fastaText

#collect further data
proteins.extend([{"Description": p.description, "Sequence": p.sequence} 
                for p in fasta.read(StringIO(multipleIDfasta))])


proteins = DataFrame(proteins)
proteins["ID"] = proteins["Description"].apply(lambda s: s.split("|")[1])#extract UniprotID
proteins["db"] = proteins["Description"].apply(lambda s: s.split("|")[0])#extract database type (SwissProt or TrEMBL)
print "Number of proteins in fasta:", len(proteins)
proteins.drop_duplicates("ID", inplace = True)
print "Remove duplicates:", len(proteins)

print "SwissProt: {}; TrEMBL: {}".format(sum(proteins["db"] == "sp"), sum(proteins["db"] == "tr"))

#writing result
fasta.write(zip(proteins["Description"], proteins["Sequence"]), open(fastaPath + "PPD.fasta", "w"))
fasta.write_decoy_db(fastaPath + "PPD.fasta", fastaPath + "PPD+reverse.fasta", file_mode = "w")
fasta.write_decoy_db(fastaPath + "PPD.fasta", fastaPath + "DECOY_PPD.fasta", file_mode = "w", decoy_only = True)

#add contaminants
ppdProteins = [p for p in fasta.read(open(fastaPath + "PPD.fasta", "r"))]
contProteins = [p for p in fasta.read(open(fastaPath + "MQCont.fasta", "r"))]
fasta.write(ppdProteins + contProteins, open(fastaPath + "PPD_MQCont.fasta", "w"))
fasta.write_decoy_db(fastaPath + "PPD_MQCont.fasta", fastaPath + "PPD_MQCont+Rev.fasta", file_mode = "w")

