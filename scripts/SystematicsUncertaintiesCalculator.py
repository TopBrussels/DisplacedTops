"""
script to merge multiple json file into a single dictionnary writen in an other json file.
The script works as follows:
1) make a list of all the .json files in the current directory
2) create a list of dict where each element corresponds to the content of one json file
3) merge everything in a single dict and save it to a new json file
4) use the merged dictionary to produce the systematics uncertainty table with the total error.

qpython 30.08.2016
"""

import xml.etree.cElementTree as ET
import json
from math import sqrt
import os  
from tabulate import tabulate

# my import
import facoLib as fl


# dict with the cst systemtic uncertainties
dict_cstSystType={"lumi":0.05, "trigger_sf_electron":0.20, "trigger_sf_muon":0.10, "trk_sf_electron":0.24, "trk_sf_muon":0.12 }


# lepton dependent syst type
muonSystTypes=["XSWeight", "evt_puSF","sf_iso_muon","sf_id_muon","lumi","trigger_sf_muon", "trk_sf_muon"]
electronSystTypes=["XSWeight", "evt_puSF","sf_reco_electron","sf_id_electron","trigger_sf_electron", "trk_sf_electron"]


# loop over the channels
for chan in ["ElEl", "MuMu"]:

    # dictionary for the cross section which depends on the samples and on the systShift 
    Yield_dict = fl.getDictFromJson(chan, "", True) 
    print "the merged dict is \n"
    print Yield_dict


    if "ElEl" in chan:
        isElEl=True
	tree = ET.ElementTree(file='../config/ElElV4.xml')
	FinalState="At least two electrons"
        print FinalState
        systTypes=list(electronSystTypes)

    elif "MuMu" in chan:
        isMuMu=True
        tree = ET.ElementTree(file='../config/MuMuV4.xml')
        FinalState="At least two muons"
        systTypes=list(muonSystTypes)
        print FinalState

    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()
	    

    # double array for table writting 
    headers=["Dataset"]
    doubleArray = []
    
    # loop over the samples 
    for sampleName in sampleNames:
    
	    singleArray = [sampleName]
    
            SumUnc=0
            CurrentUnc=0
            # loop over the systematics types 
            for systType in systTypes :
                skimedSystType=systType.replace("_electron", "")
                skimedSystType=skimedSystType.replace("_muon", "")
                headers.append(skimedSystType)
    
    
                if systType in dict_cstSystType:
                     print "faco!!!", systType
                     uncertaintyMax = dict_cstSystType[systType]
		     if "trk" in systType and "NP" not in sampleName :
		         print "removing trk sf for non signal sample"
			 uncertaintyMax = 0
                else :
                    # one variable per yield
                    YieldDown= Yield_dict[sampleName+systType+"down"]
                    YieldUp= Yield_dict[sampleName+systType+"up"]
                    YieldCentral =  Yield_dict[sampleName+"Centralup"]
                    
                    # two diff wrt to central
                    diffDown=abs((YieldCentral-YieldDown)/YieldCentral)
                    diffUp=abs((YieldCentral-YieldUp)/YieldCentral)
        
                    # max of the two effs
                    uncertaintyMax=  max(diffDown,diffUp)
    
                uncertaintyMaxPercent = str(round(100*uncertaintyMax, 2))+" %"
		if uncertaintyMax == 0:
			uncertaintyMaxPercent = "-"
		
        
                # Sum the square of the uncertainties
                CurrentUnc=uncertaintyMax
                SumUnc=SumUnc+CurrentUnc*CurrentUnc
    
                # 
                singleArray.append(uncertaintyMaxPercent)
    
    
            # eo systType loop 
            
            # add the tot uncertainty for each samples
            headers.append("total")
            singleArray.append(str(round(100*sqrt(SumUnc), 2))+" %")
            
    
            doubleArray.append(singleArray)
    
    # eo loop over the samples 
    
    fl.makeTable("SystematicsTable_"+chan, doubleArray, headers, True, "List of the systemtic uncertainties for each sample in the " + chan + " channel.", False)
    
    
    

