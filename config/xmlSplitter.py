#!/usr/bin/env python  

"""
python script to make multiple xml config out of the main xml config which contains all the samples.
This allow to run only on a subset of the samples without having to recalculate the Eq lumi and makes sure that everything is in synch.
We want to make one xml per channel (bbEl, bbMu, ElEl, MuMu).
"""


Version = "V4"
inputFile = "FullSamples"+Version+".xml"



chanString= "testFaco"
outFile=chanString+Version+".xml"
dataSetTitleList=["Data ","WJets","Diboson","Singletop","TTJets_Dilept"]


with open(inputFile,"r") as input:
     with open(outFile,"wb") as output:
          # loop over the lines
         for line in input:
              # write all empty line and comments
              if "name" not in line:
                   output.write(line)
              if "name" in line:
                   # loop over the dataset
                   for dataSet in dataSetTitleList:
                        # write only the line matching with one of the dataset in the dataset list
                        if dataSet in line :
                             output.write(line)
                             continue
