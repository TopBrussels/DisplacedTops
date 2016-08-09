#!/usr/bin/env python  

"""
python script to make multiple xml config out of the main xml config which contains all the samples.
This allow to run only on a subset of the samples without having to recalculate the Eq lumi and makes sure that everything is in synch.
We want to make one xml per channel (bbEl, bbMu, ElEl, MuMu).
"""


Version = "V4"
inputFile = "FullSamples"+Version+".xml"



channels= ["ElEl","MuMu","bbEl","bbMu"]
dataSetTitleList=[
     ["DoubleEG","500_CTau-10","WJets","Diboson","SingleTop","QCDEMEnriched","QCDbcToE","TTJets_Dilept","DrellYann"],
     ["DoubleMuon","500_CTau-10","WJets","Diboson","SingleTop","QCDMuEnriched","TTJets_Dilept","DrellYann"],
     ["SingleElectron","500_CTau-10","WJets","Diboson","SingleTop","QCDEMEnriched","QCDbcToE","TTJets_Dilept","DrellYann"],
     ["SingleMuon","500_CTau-10","WJets","Diboson","SingleTop","QCDMuEnriched","TTJets_Dilept","DrellYann"]
     ]



i_chan=0
# loop over the channels
for chan in channels:
     
     # setting outfile name
     outFileYield="Yield_"+chan+Version+".xml"
     outFileTreeProc="TreeProc_"+chan+Version+".xml"

     # open input file and outputfile
     with open(inputFile,"r") as input:
          with open(outFileYield,"wb") as output:
               # loop over the lines of the input file
               for line in input:
                    # write all empty line and comments
                    if "name" not in line:
                         output.write(line)
                    if "name" in line:
                         # loop over the dataset
                         for dataSet in dataSetTitleList[i_chan]:
                              # write only the line matching with one of the dataset in the dataset list
                              if dataSet in line :
                                   output.write(line)
                                   continue
                              
     # removes the <data> from the file to have the treeproc file
     with open(outFileYield,"r") as input:
          with open(outFileTreeProc,"wb") as output:
               for line in input:
                    if "data>" not in line:
                         output.write(line)



     i_chan+=1
# eo loop over the channels 
