#!/usr/bin/env python
"""
This script will make a copy of some specific tree in a root file a fill them only if a certain condition is met.
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
#regions=["PCR", "DCR", "SR1","SR2", "SR3"]
regions=["PCR", "DCR"]
# corresponding bounds
bounds=[0., 0.01, 0.02 , 0.5, 0.1]
# corresponding bool 
bools=[True, True, True, True, True]
# corresponding trees





isElEl=False
isMuMu=True

# loading the xml
tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
root = tree.getroot()
datasets = root.find('datasets')
print "found  "  + str(len(datasets)) + " datasets"
datasetNames = []

chan="_MuMu"

idataset=0
for d in datasets:
    if d.attrib['add'] == '1' :
        datasetNames.append(str(d.attrib['name']))
        print str(d.attrib['name'])
            # one array per dataset [name, title, Eqlumi, N1, N2, N3, SR1, SR2, SR3]    
        sampleName=d.attrib['name']

        # build the chain
        ch_in = ROOT.TChain("doubleMuTree","doubleMuTree")
        input_files = pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root"
        ch_in.Add(input_files)




## max_events

        max_events = ch_in.GetEntries()

#if ops.max_events > 0 and max_events > ops.max_events:

#    max_events = ops.max_events



## branch status

#set_status_of_branches(ch_in, ops.branches_on_file, ops.branches_off_file)



## start new file

        # create one root file for each sample X channel X regions
#        for i_reg in range(0,len(regions)):
#            new_file = ROOT.TFile(sampleName+chan+regions[i_reg]+"Skimmed.root", 'RECREATE')
#            ch_out = ch_in.CloneTree(0)


        new_file = ROOT.TFile(sampleName+chan+"Skimmed.root", 'RECREATE')

        # 
        trees=[]       
        for i_reg in range(0,len(regions)):
            trees.append(ch_in.CloneTree(0))

#            print trees[i_reg]
            trees[i_reg].SetName(ch_in.GetName()+regions[i_reg])

        # bo loop over the event
        for i_event in xrange(max_events):
    
            i_entry = ch_in.LoadTree(i_event)
            ch_in.GetEntry(i_event)
            

            
            # Define one bool per cut on d0
            # This defines 5 inclusive regions and exclusive region can be defined by requiring to to pass cut x and failling cut x+1
            # The exclusive regions the following: Promt Control Region (PCR), Displaced CR (DCR), Singal Region (SR) 1, 2 and 3.

            bools=[True, True, True, True, True]
            
            passed1= True
            passed2= True
            passed3= True

            # make printout every 10 000 events
            if i_event % 10000 == 0:
                print 'Processing event %i of %i' % (i_event, max_events)

            #  skip 99% of the events just to run faster
            if not i_event % 100 == 0:
                continue

             # loop over the 2 highest pt letpon                                                                                                                                
            for ilept in range (0,2):
                
                # make the logic for the muon                                                                                                                                  
                if isMuMu:
#                    LeptonWeight *= iev.sf_muon_mumu[ilept]

                    # if one of the leptons  is smaller than bound, the event fails                          

                    # looping over all the regions
                    for i_reg in range(0,len(regions)):
#                        print "bound is ", bounds[i_reg]
#                        print "bool is ", bools[i_reg]
#                        print "region is ", regions[i_reg]
                        
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
                        #for bound in range(0,len(SRxBounds)): 


                        # if one of the leptons  is smaller than bound, the event fails                                                                                            
                    if abs(ch_in.d0BeamSpot_electron_elel[ilept]) < bound1:
#                        passed1=False
                        if (debug):
                            print "Electron and muon entering N1"
                            print "d0 electron is " , ch_in.d0BeamSpot_electron_elel[ilept]


                    # eo the logic for the electron


            passed_skim = False
            if (ch_in.pt_muon_mumu[0] > 70):
                #print ch_in.pt_muon_mumu[0]
                passed_skim = True
    


            if passed_skim:
                trees[0].Fill()
            if (ch_in.pt_muon_mumu[0] > 150):
                trees[1].Fill()

        # eo loop over the event 

        #    ch_out.Print() 
        for i_reg in range(0,len(regions)):
            trees[i_reg].GetCurrentFile().Write()

        trees[0].GetCurrentFile().Close()
                                                   
        
