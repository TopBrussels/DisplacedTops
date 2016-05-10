#!/usr/bin/env python
"""
This script will make a copy of some specific tree in a root file and fill them only if a certain condition is met.
This can be viewed as a skimer.
uses: python tree_trimmer.py 

This script was largely inspired by an other script that you can find there:
https://svnweb.cern.ch/trac/penn/browser/reece/rel/tree_trimmer/tags/tree_trimmer-00-01-00?order=name

"""

#------------------------------------------------------------------------------



## std

import optparse
import fnmatch
import os, sys
import xml.etree.cElementTree as ET


## ROOT

#import ROOT, rootlogo
import ROOT 
ROOT.gROOT.SetBatch(True)



#_____________________________________________________________________________



# path to root trees                                                                          
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"
#date="15_4_2016"
date="NoDisplacedTrigger"

# debug                                                                                                                                                                            
debug=False




# regions
#regions=["PCR", "DCR", "SR1","SR2", "SR3"]
regions=["PCR"]
# corresponding bool 
#bools=[True, True, True, True, True]
# This defines 5 inclusive regions and exclusIve region can be defined by requiring to to pass cut x and failling cut x+1
# The exclusive regions the following: Promt Control Region (PCR), Displaced CR (DCR), Singal Region (SR) 1, 2 and 3.



isElEl=False
isMuMu=True

# loading the xml
datasetNames = []

channels=["_ElEl","_MuMu"]
#channels=["_MuMu"]

# loop over the channel (lepton in final statue)                                       
for chan in channels:


    isElEl=False
    isMuMu=False

    treeName="tree"
    # getting the correct xml file depending on the final state
    if "MuMu" in chan:
        isMuMu=True
        tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')
        FinalState="At least two electrons"
        print FinalState
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/Yield_FullSamplesElMuV0.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    idataset=0


    for d in datasets:
        if d.attrib['add'] == '1' :
            datasetNames.append(str(d.attrib['name']))
            print str(d.attrib['name'])
            sampleName=d.attrib['name']
    
            # build the chain
            ch_in = ROOT.TChain(treeName,treeName)
            input_files = pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root"
            ch_in.Add(input_files)
    
            # max_events
            max_events = ch_in.GetEntries()
    
    
    
            # start new file for each sample and each channel
            new_file = ROOT.TFile(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+"Skimmed.root", 'RECREATE')
    
            # create one tree per regions in the current file
            
            trees=[]       
            directories=[]
            for i_reg in range(0,len(regions)):
    #            var=ROOT.TTree
                trees.append(ch_in.CloneTree(0))
#                directory=new_file.mkdir(ch_in.GetName()+regions[i_reg])
#                directories.append(directory)
#                trees[i_reg].SetName(ch_in.GetName()+regions[i_reg])
#                trees[i_reg].SetDirectory(directories[i_reg])
                
    #            print trees[i_reg]
    
    
            # bo loop over the event
            for i_event in range(max_events):
        
                i_entry = ch_in.LoadTree(i_event)
                ch_in.GetEntry(i_event)
                
    
                
                # Define one bool per cut on d0
    
                bools=[True, True]
                
                # make printout every 10 000 events
                if i_event % 10000 == 0:
                    print 'Processing event %i of %i' % (i_event, max_events)
    
                #  skip 99% of the events just to run faster
#                if not i_event % 100 == 0:
#                    continue
                    
                if isElEl:
                    nLept=ch_in.nElectrons
                    nLeptPair=ch_in.nElectronPairs
                if isMuMu:
                    nLept=ch_in.nMuons
                    nLeptPair=ch_in.nMuonPairs

                # loop over all the leptons to check the d0
                for ilept in range (0,nLept):
                    
                    # make the logic for the muon 
                    if isMuMu:
    
                        # if one of the leptons  is smaller than bound, the event fails                          
                        if abs(ch_in.d0BeamSpot_muon[ilept])> 0.01:
#                            bools[i_reg]=False
                            bools[0]=False
                    # eo the logic for the muon 
    
    
    
    
                    # make the logic for the electron 
                    if isElEl :
                        if abs(ch_in.d0BeamSpot_electron[ilept]) > 0.01:
#                            print "d0BeamSpot_electron[ilept] is " , ch_in.d0BeamSpot_electron[ilept]
#                            bools[i_reg]=False
                            bools[0]=False
                    # eo the logic for the electron


                

                # loop over all the lepton pairs and check the invmass
                for ileptPair in range (0,nLeptPair):
                    # reject if outside the Z peak

                    if isMuMu :
                        if  ch_in.invMass_mumu[ileptPair] <= 81.2 or 101.2 <= ch_in.invMass_mumu[ileptPair] :
                            bools[0]=False
    
                    if isElEl :
                        if  ch_in.invMass_elel[ileptPair] <= 81.2 or 101.2 <= ch_in.invMass_elel[ileptPair] :
                            bools[0]=False
    
                                    
                # fill the tree if the condition is passed
                if bools[0] ==True:
#                if (True):
                     trees[i_reg].Fill()
    
            # eo loop over the event 
    
            # to be FIXED only write SR for data.
            # to be FIXED logic for CR...
            # write the tree
            for i_reg in range(0,len(regions)):
#                directories[i_reg].cd()
                trees[i_reg].Write()
#                trees[i_reg].GetCurrentFile().Close()

            # closing files
            new_file.Close()
            ch_in.GetCurrentFile().Close()
    
    # eo loop over dataset
                                                       
        
# eo loop over final state
