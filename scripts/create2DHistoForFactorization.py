############
# Pyroot macro to create d0 histograms from bb+El and bb+Mu Control region. These histograms will be used in the ClosureTestQCDEstimate.py script in order to caclculate Transfer Factors (TFs) to make the closure test.
# April 2016 by qpython@cern.ch 

import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt


#channel                                                                                                                       
#channels=["_MuMu"]                                                                                                            
channels=["_ElEl"]                                                                                                            
#channels=["_ElEl","_MuMu"]


#root file  
#radix of the path
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"
#date
date="15_4_2016"


# output root file
outfile = rt.TFile("rootFiles/2D.root",'RECREATE')

# bool
debug = False

# lumi
lumivalue = 1
weight = 1

# define do histograms
electrond0VsElectronsd0=rt.TH2D("electrond0VsElectronsd0","electrond0VsElectronsd0", 50, 0.0, 0.05, 50, 0.0, 0.05)
muond0VsMuond0=rt.TH2D("muond0VsMuond0","muond0VsMuond0",50, 0.0, 0.05, 50, 0.0, 0.05)
muond0VsElectrond0=rt.TH2D("muond0VsElectrond0","muond0VsElectrond0",50, 0.0, 0.05, 50, 0.0, 0.05)



# to be fixed if you want to run on MC instead of Data

#tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
#treeName="doubleMuTree"

tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')
treeName="doubleElTree"

treeName="tree"



for chan in channels :
    isElEl =False
    isMuMu =False


    if chan == "_ElEl" :
        isElEl = True
    if chan == "_MuMu" :
        isMuMu = True
    
    
    root =  tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    
    # loop over the composite datasets (title)
    for i in range(0,1): # to be removed!!
        

        # loop over the dataset inside the composite dataset (name)
        for d in datasets:

            # only for Data
            if d.attrib['add'] == '1': 
#            if d.attrib['add'] == '1' and compositeDataset == str(d.attrib['title']): 
                datasetNames.append(str(d.attrib['name']))
                print str(d.attrib['name'])
                ch = rt.TChain(treeName,treeName)
                sampleName=d.attrib['name']

#                outfile = rt.TFile("rootFiles/"+sampleName+"2D.root",'RECREATE')

                # fix the type of dataset (bgMC, signal or data)                           
                isBgMC = False
                isSignal = False
                isData = False

                if "Data" not in sampleName and "NP" not in sampleName:
                    isBgMC = True
                if "NP" in sampleName:
                    isSignal = True
                if "Data" in sampleName:
                    isData = True

        
                ch.Add(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")
    
        
                nevents=ch.GetEntries()
        
        
                if isData :
                    lumivalue=float(d.attrib['EqLumi'])
                        
                    weight= lumivalue / float(d.attrib['EqLumi'])
                if (1):
                    print "lumivalue is " ,lumivalue
                    print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                    print "weight is " , weight
        
                ii=0
                # start of loop over events                                                                                                                
                for iev in ch:
    
                    if isMuMu:
                        PileUpWeight=iev.evt_puSF
                    if isElEl:
                        PileUpWeight=iev.evt_puSF
        
                    
                    LeptonWeight=1.0
                    
                    for ilept in range (0,2):
                        if isMuMu:
                            LeptonWeight *= iev.sf_muon[0]
                        if isElEl:
                            LeptonWeight *= iev.sf_electron[0]
                    
        
        
                    # define shorter variable                                                                                                              
                    if isMuMu:
                        d01mu=abs(iev.d0BeamSpot_muon[0])
                        d02mu=abs(iev.d0BeamSpot_muon[1])
                        muond0VsMuond0.Fill(d01mu,d02mu,weight*PileUpWeight*LeptonWeight)
                        if (debug):
                            print "d0 is " , d0mu
                            print "weight is ", weight
                            print "PileUpWeight is ", PileUpWeight
                            print "LeptonWeight is ", LeptonWeight
                    if isElEl :
                        d01el=abs(iev.d0BeamSpot_electron[0])
                        d02el=abs(iev.d0BeamSpot_electron[1])
#                        electrond0VsElectronsd0.Fill(d01el,d02el,weight*PileUpWeight*LeptonWeight)
                        electrond0VsElectronsd0.Fill(d01el,d02el)
        

    

# end of loop


# write the histo in the output file
outfile.cd()
muond0VsMuond0.Write()
electrond0VsElectronsd0.Write()
outfile.Close()
