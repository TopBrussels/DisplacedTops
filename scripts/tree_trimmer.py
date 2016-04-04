#!/usr/bin/env python
"""
This script will make a copy of some specific tree in a root file and fill them only if a certain condition is met.
This can be view as a skimer.
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
date="2_3_2016"

# debug                                                                                                                                                                            
debug=False




# regions
regions=["PCR", "DCR", "SR1","SR2", "SR3"]
#regions=["PCR", "DCR"]
# corresponding bounds
bounds=[0., 0.01, 0.02 , 0.5, 0.1]
# corresponding bool 
bools=[True, True, True, True, True]
# This defines 5 inclusive regions and exclusIve region can be defined by requiring to to pass cut x and failling cut x+1
# The exclusive regions the following: Promt Control Region (PCR), Displaced CR (DCR), Singal Region (SR) 1, 2 and 3.



isElEl=False
isMuMu=True

# loading the xml
datasetNames = []

channels=["_ElEl","_MuMu"]

# loop over the channel (lepton in final statue)                                       
for chan in channels:


    isElEl=False
    isMuMu=False

    # getting the correct xml file depending on the final state
    if "MuMu" in chan:
        isMuMu=True
        tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
        treeName="doubleMuTree"
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')
        treeName="doubleElTree"
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
            new_file = ROOT.TFile(pathTrunc+date+"/"+chan+"/"+sampleName+chan+"Skimmed.root", 'RECREATE')
    
            # create one tree per regions in the current file
            
            trees=[]       
            directories=[]
            for i_reg in range(0,len(regions)):
    #            var=ROOT.TTree
                trees.append(ch_in.CloneTree(0))
                directory=new_file.mkdir(ch_in.GetName()+regions[i_reg])
                directories.append(directory)
#                trees[i_reg].SetName(ch_in.GetName()+regions[i_reg])
                trees[i_reg].SetDirectory(directories[i_reg])
                
    #            print trees[i_reg]
    
    
            # bo loop over the event
            for i_event in xrange(max_events):
        
                i_entry = ch_in.LoadTree(i_event)
                ch_in.GetEntry(i_event)
                
    
                
                # Define one bool per cut on d0
    
                bools=[True, True, True, True, True]
                
                # make printout every 10 000 events
                if i_event % 10000 == 0:
                    print 'Processing event %i of %i' % (i_event, max_events)
    
                #  skip 99% of the events just to run faster
#                if not i_event % 100 == 0:
#                    continue
    
                 # loop over the 2 highest pt letpon
                for ilept in range (0,2):
                    
                    # make the logic for the muon 
                    if isMuMu:
    #                    LeptonWeight *= iev.sf_muon_mumu[ilept]
    
                        # if one of the leptons  is smaller than bound, the event fails                          
    
                        # looping over all the regions
                        for i_reg in range(0,len(regions)):
                            if (debug):
                                print "bound is ", bounds[i_reg]
                                print "bool is ", bools[i_reg]
                                print "region is ", regions[i_reg]
                            
                            if abs(ch_in.d0BeamSpot_muon_mumu[ilept]) <= bounds[i_reg]:
                                bools[i_reg]=False
                            else :
                                if (debug):
                                    print "one muon passed the cut to enter ", regions[i_reg]  
                                    print "d0 muon is " , ch_in.d0BeamSpot_muon_mumu[ilept]
                    # eo the logic for the muon 
    
    
    
    
                    # make the logic for the electron 
                    if isElEl :
    #                    LeptonWeight *= ch_in.sf_electron_elel[ilept]
    
                        # looping over all the regions
                        for i_reg in range(0,len(regions)):
                            if (debug):
                                print "bound is ", bounds[i_reg]
                                print "bool is ", bools[i_reg]
                                print "region is ", regions[i_reg]
    
                            if abs(ch_in.d0BeamSpot_electron_elel[ilept]) <= bounds[i_reg]:
                                bools[i_reg]=False
                            else :
                                if (debug):
                                    print "one electron passed the cut to enter ", regions[i_reg]
                                    print "d0 electron is " , ch_in.d0BeamSpot_electron_elel[ilept]
                        # eo the logic for the electron
    
    
                                    
                # filling the correct tree according the the vector of bools
                for i_reg in range(0,len(regions)-1):
                    if bools[i_reg] == True and bools[i_reg+1] == False:
                        trees[i_reg].Fill()
                        continue
                    else :
                        trees[i_reg+1].Fill()
    
            # eo loop over the event 
    
            # to be FIXED only write SR for data.
            # to be FIXED logic for CR...
            # write the tree
            for i_reg in range(0,len(regions)):
                directories[i_reg].cd()
                trees[i_reg].Write()
#                trees[i_reg].GetCurrentFile().Close()

            # closing files
            new_file.Close()
            ch_in.GetCurrentFile().Close()
    
    # eo loop over dataset
                                                       
        
# eo loop over final state
