# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 20:55:00 2017

Batch execution of LC-MS analysis workflow

Version 2.0

@author: vgor
"""
from subprocess import check_call as call
from coronaHelpers import readPercOut, splitIDs, mapQValues

parentfolder = "E:\\RawData\\20161120_NP"
fastaLocation = "E:\\Fasta\\PPD_MQCont+Rev.fasta"
msconvertMask = r'"E:\Pwiz\msconvert.exe" --32 --mzML -z --filter "peakPicking true 1-" --filter "threshold absolute 0 most-intense" -o "{pf}" --outfile "{pf}\{name}.mzML" "{pf}\{name}.raw"'
msgfMask = r'"C:\Program Files\Java\jre1.8.0_152\bin\java" -Xms128M -Xmx8096M -jar "E:\MSSoft\MSGFPlus\MSGFPlus.jar" -s "{pf}\{name}.mzML" -d "{db}" -o "{pf}\{name}.mzid" -t 10.0ppm -tda 0 -mod "{pf}\msgf.mods" -minCharge 2 -maxCharge 6 -inst 3 -thread 8 -m 3 -e 1 -ntt 2 -protocol 0 -minLength 6 -maxLength 40 -n 1 -addFeatures 1 -ti "0,1"'
pconvMask = r'"E:\MSSoft\percolator-converters-v3-01\bin\msgf2pin.exe" -o "{pf}\{name}.pin" -P DECOY_ -F "{db}" -c 2 -m 1 "{pf}\{name}.files"'
percMask = r'"E:\MSSoft\percolator-v3-01\bin\percolator.exe" -Y -c -g -P DECOY_ -f "{db}" -z trypsin -X "{pf}\{name}.pout.xml" "{pf}\{name}.pin" > "{pf}\{name}.pout"'

##Temperature experiments
sampleNames = ["Temp{}{}_{}".format(point, bRep, tRep) for point in "ABCDE" for bRep in "123C" for tRep in "123"]
with open(parentfolder + "\\Temperature.files", "w") as fout:
    for name in sampleNames:
        fout.write("{pf}\\{name}.mzid\n".format(pf = parentfolder, name = name))
        execline = msconvertMask.format(pf = parentfolder, name = name)
        print(execline)
        call(execline, shell = True)
        execline = msgfMask.format(pf = parentfolder, name = name, db = fastaLocation)
        print(execline)
        call(execline, shell = True)
        
execline = pconvMask.format(pf = parentfolder, name = "Temperature", db = fastaLocation)
print(execline)
call(execline, shell = True)
execline = percMask.format(pf = parentfolder, name = "Temperature", db = fastaLocation)
print(execline)
call(execline, shell = True)

psms = readPercOut("{}\\Temperature.pout.xml".format(parentfolder))[0]
splitIDs(psms)
print("{} PSMs read from {}\\Temperature.pout.xml".format(len(psms), parentfolder))
for name in sampleNames:
    mzidName = "{}\\{}.mzid".format(parentfolder, name)
    mapQValues(mzidName, psms[psms["File"] == name], mzidName)

#pH experiments
sampleNames = ["pH{}_{}T{}".format(point, bRep, tRep) for point in ["3", "5", "6_6", "7", "9"]
                                                     for bRep in ["B1", "B2", "B3", "C"]
                                                     for tRep in "123"]

with open(parentfolder + "\\pH.files", "w") as fout:
    for name in sampleNames:
        fout.write("{pf}\\{name}.mzid\n".format(pf = parentfolder, name = name))
        execline = msconvertMask.format(pf = parentfolder, name = name)
        print(execline)
#        call(execline, shell = True)
        execline = msgfMask.format(pf = parentfolder, name = name, db = fastaLocation)
        print(execline)
#        call(execline, shell = True)
        
execline = pconvMask.format(pf = parentfolder, name = "pH", db = fastaLocation)
print(execline)
call(execline, shell = True)
execline = percMask.format(pf = parentfolder, name = "pH", db = fastaLocation)
print(execline)
call(execline, shell = True)

psms = readPercOut("{}\\pH.pout.xml".format(parentfolder))[0]
splitIDs(psms)
print("{} PSMs read from {}\\pH.pout.xml".format(len(psms), parentfolder))
for name in sampleNames:
    mzidName = "{}\\{}.mzid".format(parentfolder, name)
    mapQValues(mzidName, psms[psms["File"] == name], mzidName) 
            