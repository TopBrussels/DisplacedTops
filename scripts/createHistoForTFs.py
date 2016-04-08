############
# Pyroot macro to create d0 histograms from bb+El and bb+Mu Control region. These histograms will be used in the ClosureTestQCDEstimate.py script in order to caclculate Transfer Factors (TFs) to make the closure test.
# April 2016 by qpython@cern.ch 

import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt

#root file of QCD only tree                                                                                                
#date
date="24_3_2016"


# output root file
outfile = rt.TFile("rootFiles/d0forTFs.root",'RECREATE')


#leptons
leptons=["el","mu"]

# bool
debug = False


# define do histograms
electrond0=rt.TH1D("electrond0","electrond0",10,0.01,0.019999999)
muond0=rt.TH1D("muond0","muond0",10,0.01,0.019999999999)



# to be fixed if you want to run on MC instead of Data
treeTF = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
treeName="preCutTree"


for lepton in leptons:
    isElEl =False
    isMuMu =False


    if lepton == "el" :
        isElEl = True
    if lepton == "mu" :
        isMuMu = True
    
    
    rootTF =  treeTF.getroot()
    datasetsTF = rootTF.find('datasets')
    print "getting dataset to get the Transfer Factors"
    print "found  "  + str(len(datasetsTF)) + " datasets"
    datasetNamesTF = []
    
        
    for d in datasetsTF:
        # only for Data
        if d.attrib['add'] == '1' and "Data" in d.attrib['name'] :
            datasetNamesTF.append(str(d.attrib['name']))
            print str(d.attrib['name'])
            chTF = rt.TChain(treeName,treeName)
            sampleName=d.attrib['name']
    
    
            # fix the type of dataset (bgMC, signal or data)                                                                                           
            isBgMC = False
            isQCD = False
            isData = False
            
            if "Data" not in sampleName and "NP" not in sampleName:
                isBgMC = True
            if "QCD" in sampleName:
                isQCD = True
            if "Data" in sampleName:
                isData = True
    
    
            if (isMuMu):
                chTF.Add("/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"+date+"/_bbMu/DisplacedTop_Run2_TopTree_Study_Data_bbMu.root")
            if (isElEl):
                chTF.Add("/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"+date+"/_bbEl/DisplacedTop_Run2_TopTree_Study_Data_bbEl.root")
    
    
    
                #                    chTF.Add(rootFilebbMu)                                                                                         
            nevents=chTF.GetEntries()
    
    
    
                # calculate weight 
        
            if isData :
                lumivalue=float(d.attrib['EqLumi'])
                    
                weight= lumivalue / float(d.attrib['EqLumi'])
            if (1):
                print "lumivalue is " ,lumivalue
                print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                print "weight is " , weight
    
            ii=0
            # start of loop over events                                                                                                                
            for iev in chTF:
                    
#                if isMuMu:
#                    PileUpWeight=iev.evt_puSF_pc
#                if isElEl:
#                    PileUpWeight=iev.evt_puSF_pc
    
                # temporarily set to one
                PileUpWeight=1
                        
                LeptonWeight=1.0
                
                # bo the logic for the muon                                                                                                        
                if isMuMu:
                    LeptonWeight *= iev.sf_muon_pc[0]
                if isElEl:
                    LeptonWeight *= iev.sf_electron_pc[0]
    
                # temporarily set to one    
                LeptonWeight = 1
    
    
                # define shorter variable                                                                                                              
                if isMuMu:
                    d0mu=abs(iev.d0BeamSpot_muon_pc[0])
                    muond0.Fill(d0mu,1)
                    if (0.01 < d0mu) and (d0mu < 0.02) :
                        if (debug):
                            print "d0 is " , d0mu
                            print "weight is ", weight
                            print "PileUpWeight is ", PileUpWeight
                            print "LeptonWeight is ", LeptonWeight
                if isElEl :
                    d0el=abs(iev.d0BeamSpot_electron_pc[0])
                    electrond0.Fill(d0el,1)
    

    

# end of loop


# write the histo in the output file
outfile.cd()
electrond0.Write()
muond0.Write()
outfile.Close()
