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
#date="NoBlindingRerun_30_11_2016"
date="Systematics_29_8_2016"

# debug 
debug=False
fastRun=True

# dictionary to link region Name with bound
region_dict={"PCR": {"lb": 0.00 , "ub": 0.01},
             "DCR": {"lb": 0.01 , "ub": 0.02},
             "SRs": {"lb": 0.02 , "ub": 10},
             "OffZ": {"lb": 0.00 , "ub": 10}
             }

#print region_dict["DCR"]["lb"]

# regions
#regions=["PCR", "DCR", "SRs"]
regions=["DCR"]


# loading the xml
datasetNames = []

#channels=["_ElEl","_MuMu"]
channels=["_MuMu"]
#channels=["_ElEl"]

# loop over the channel (lepton in final state) 
for chan in channels:

    isElEl=False
    isMuMu=False

    treeName="tree"
    # getting the correct xml file depending on the final state
    if "MuMu" in chan:
        isMuMu=True
        tree = ET.ElementTree(file='../config/MuMuV4.xml')
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        tree = ET.ElementTree(file='../config/ElElV4.xml')
        FinalState="At least two electrons"
        print FinalState
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/ElMuV4.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    idataset=0

    

    # loop over the dataset
    for d in datasets:
        title= d.attrib["title"]
#        if d.attrib['add'] == '1' and ("Data" in title  or "DrellYann" in title ) and "QCD" not in title  :
#        if d.attrib['add'] == '1' and "QCD" not in title and "DrellYann" in title  :
        if d.attrib['add'] == '1' and "QCD" not in title and "Data" not in title and "stop" not in title :
            datasetNames.append(str(d.attrib['name']))
            print "\n"+str(d.attrib['name'])
            sampleName=d.attrib['name']

    
            # build the chain
            ch_in = ROOT.TChain(treeName,treeName)
            input_files = pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root"
#            print input_files
            ch_in.Add(input_files)
    
            # max_events
            max_events = ch_in.GetEntries()

    
            # loop over the region
            for i_reg, reg in enumerate(regions):

                # counter 
                n_selected = 0

                # output file
                new_file = ROOT.TFile(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+reg+".root", 'RECREATE')
                mytree=ch_in.CloneTree(0)

                # easier variable
                lb = region_dict[reg]["lb"] 
                ub = region_dict[reg]["ub"]
                if debug :
                    print "lb is ", lb
                    print "ub is ", ub

    
                # bo loop over the event
                for i_event in range(max_events):
            
                    i_entry = ch_in.LoadTree(i_event)
                    ch_in.GetEntry(i_event)

                    # Define one bool per event. If it fails one of the requirement it will be rejected.
                    keep = True

                    # make printout every 10 000 events
                    if i_event % 10000 == 0:
                        print 'Processing event %i of %i' % (i_event, max_events)
        
                    #  skip 99% of the events just to run faster
                    if fastRun and not i_event % 100 == 0:
                        continue
                        
                    if isElEl:
                        nLept=ch_in.nElectrons
                        nLeptPair=ch_in.nElectronPairs
                        leptInvMAss=ch_in.invMass_elel[0]

                    if isMuMu:
                        nLept=ch_in.nMuons
                        nLeptPair=ch_in.nMuonPairs
                        leptInvMAss=ch_in.invMass_mumu[0]

                    # removes contribution from Z if offZ region
                    if reg[i_reg] == "OffZ":
                        if  81.2 <=  leptInvMAss <= 101.2 :
                            keep=False
                        continue
    
    
                    # loop over all the leptons to check the d0
                    for ilept in range (0,nLept):

                        if isElEl:
                            d0_lept = abs(ch_in.d0BeamSpot_electron[ilept])

                        if isMuMu:
                            d0_lept = abs(ch_in.d0BeamSpot_muon[ilept])
                        
                        # if one of the leptons is outside the [lb;ub] region we reject the event
                        if d0_lept < lb or ub < d0_lept:
                            keep=False
                            continue
    

    
    
        
    
        
                                        
                    # fill the tree if the condition is passed
                    if keep == True:
                        n_selected+=1
                        mytree.Fill()
                

                # eo loop over the event 
#                print "n_selected is ", n_selected


                # write the tree
                mytree.Write()

                # closing files
                new_file.Close()
                ch_in.GetCurrentFile().Close()

            # eo loop over the regions

    # eo loop over dataset
        
# eo loop over final state
