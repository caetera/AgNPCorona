# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 20:55:00 2017

@author: vgor
"""
from subprocess import check_call as call
from coronaHelpers import readPercOut, splitIDs, mapQValues
from pandas import DataFrame

parentfolder = "E:\\RawData\\20161120_NP"
fastaLocation = "E:\\Fasta\\PPD_MQCont+Rev.fasta"
msconvertMask = r'"E:\Pwiz\msconvert.exe" --32 --mzML -z --filter "peakPicking true 1-" --filter "threshold absolute 0 most-intense" -o "{pf}" --outfile "{pf}\{name}.mzML" "{pf}\{name}.raw"'
msgfMask = r'"C:\Program Files\Java\jre1.8.0_31\bin\java" -Xms128M -Xmx8096M -jar "E:\MSSoft\MSGFPlus\MSGFPlus.jar" -s "{pf}\{name}.mzML" -d "{db}" -o "{pf}\{name}.mzid" -t 10.0ppm -tda 0 -mod "{pf}\msgf.mods" -minCharge 2 -maxCharge 6 -inst 3 -thread 8 -m 3 -e 1 -ntt 2 -protocol 0 -minLength 6 -maxLength 40 -n 1 -addFeatures 1 -ti "0,1"'
pconvMask = r'"E:\MSSoft\percolator-converters-v3-01\bin\msgf2pin.exe" -o "{pf}\{name}.pin" -P DECOY_ -F "{db}" -c 2 -m 1 "{pf}\{name}.files"'
percMask = r'"E:\MSSoft\percolator-v3-01\bin\percolator.exe" -Y -c -g -P DECOY_ -f "{db}" -z trypsin -X "{pf}\{name}.pout.xml" "{pf}\{name}.pin" > "{pf}\{name}.pout"'

knimeCols = ["mzid" + a for a in "012"] + ["mzML" + a for a in "012"] + ["featureOut", "peptideOut", "proteinOut"]
knimeData = []

##Temperature experiments
for point in "ABCDE":
    for bRep in "123C":
        sname = "Temp{}{}".format(point, bRep)
        with open(parentfolder + "\\{}.files".format(sname), "w") as fout:
            for tRep in "123":
                name = "{}_{}".format(sname, tRep)
                fout.write("{pf}\\{name}.mzid\n".format(pf = parentfolder, name = name))
                execline = msconvertMask.format(pf = parentfolder, name = name)
                print execline
                call(execline, shell = True)
                execline = msgfMask.format(pf = parentfolder, name = name, db = fastaLocation)
                print execline
                call(execline, shell = True)
        
        execline = pconvMask.format(pf = parentfolder, name = sname, db = fastaLocation)
        print execline
        call(execline, shell = True)
        execline = percMask.format(pf = parentfolder, name = sname, db = fastaLocation)
        print execline
        call(execline, shell = True)
        
        psms = readPercOut("{}\\{}.pout.xml".format(parentfolder, sname))[0]
        splitIDs(psms)
        print "{} PSMs read from {}\\{}.pout.xml".format(len(psms), parentfolder, sname)
        for tRep in "123":
            name = "{}_{}".format(sname, tRep)
            mzidName = "{}\\{}.mzid".format(parentfolder, name)
            mapQValues(mzidName, psms[psms["File"] == name], mzidName)
        
        knimeData.append(["{}\\{}_{}.mzid".format(parentfolder, sname, n) for n in "123"] \
                    + ["{}\\{}_{}.mzML".format(parentfolder, sname, n) for n in "123"] \
                    + ["{}\\{}.featureXML".format(parentfolder, sname), \
                       "{}\\{}.peptides.csv".format(parentfolder, sname), \
                       "{}\\{}.proteins.csv".format(parentfolder, sname)])


#pH experiments
for point in ["3", "5", "6_6", "7", "9"]:
    for bRep in ["B1", "B2", "B3", "C"]:
        sname = "pH{}_{}".format(point, bRep)
        with open(parentfolder + "\\{}.files".format(sname), "w") as fout:
            for tRep in "123":
                name = "{}T{}".format(sname, tRep)
                print name
                fout.write("{pf}\\{name}.mzid\n".format(pf = parentfolder, name = name))
                execline = msconvertMask.format(pf = parentfolder, name = name)
                print execline
                call(execline, shell = True)
                execline = msgfMask.format(pf = parentfolder, name = name, db = fastaLocation)
                print execline
                call(execline, shell = True)
        
        execline = pconvMask.format(pf = parentfolder, name = sname, db = fastaLocation)
        print execline
        call(execline, shell = True)
        execline = percMask.format(pf = parentfolder, name = sname, db = fastaLocation)
        print execline
        call(execline, shell = True)
        
        psms = readPercOut("{}\\{}.pout.xml".format(parentfolder, sname))[0]
        splitIDs(psms)
        print "{} PSMs read from {}\\{}.pout.xml".format(len(psms), parentfolder, sname)
        for tRep in "123":
            name = "{}T{}".format(sname, tRep)
            mzidName = "{}\\{}.mzid".format(parentfolder, name)
            mapQValues(mzidName, psms[psms["File"] == name], mzidName)
        
        knimeData.append(["{}\\{}T{}.mzid".format(parentfolder, sname, n) for n in "123"] \
                    + ["{}\\{}T{}.mzML".format(parentfolder, sname, n) for n in "123"] \
                    + ["{}\\{}.featureXML".format(parentfolder, sname), \
                       "{}\\{}.peptides.csv".format(parentfolder, sname), \
                       "{}\\{}.proteins.csv".format(parentfolder, sname)])

DataFrame(knimeData, columns = knimeCols).to_csv("{}\\knimeInput.csv".format(parentfolder), index = False)