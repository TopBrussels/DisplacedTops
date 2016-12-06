############
# Pyroot macro to create d0 histograms from bb+El and bb+Mu Control region. These histograms will be used in the ClosureTestQCDEstimate.py script in order to caclculate Transfer Factors (TFs) to make the closure test.
# You can also get the TF from the DCR to the SR using the overflow bin. TF(DCR->SR) = N signal/ N DCr = int N+1/int 1->N
# April 2016 by qpython@cern.ch 

import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt


# trunc
trunc = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

# folder name
folder="QCDSyst_20_9_2016/"


# output root file
outfile = rt.TFile("rootFiles/d0forTFsNew.root",'RECREATE')


#leptons
#leptons=["el","mu"]
leptons=["mu","el"]

# bool
debug = False


# list of btag WP
wps = ["_Loose", "_Medium", "_Tight"]


# loop over the wp
for wp in wps :


    # define d0 histograms For Closure Test
    electrond0=rt.TH1D("electrond0"+wp,"electrond0"+wp,10,0.01,0.02)
    muond0=rt.TH1D("muond0"+wp,"muond0"+wp,10,0.01,0.02)
    
    
    # define d0 histograms for QCD estimate
    d0_bins = sorted([0.01, 0.02, 0.05, 0.1, 10])
    d0_array = array('d',d0_bins)
    
    
    electrond0Wide=rt.TH1D("electrond0Wide"+wp,"electrond0Wide"+wp, len(d0_array)-1, d0_array)
    muond0Wide=rt.TH1D("muond0Wide"+wp,"muond0Wide"+wp, len(d0_array)-1, d0_array)
    
    
    
    # name of the tree
    treeName="tree"
    
    
    for lepton in leptons:
        isElEl =False
        isMuMu =False
    
    
        if lepton == "el" :
            isElEl = True
        if lepton == "mu" :
            isMuMu = True
    
        chTF = rt.TChain(treeName,treeName)
        
        if (isMuMu):
            chTF.Add(trunc+folder+wp+"_bbMu/DisplacedTop_Run2_TopTree_Study_DataRunD_bbMu.root")
        if (isElEl):
            chTF.Add(trunc+folder+wp+"_bbEl/DisplacedTop_Run2_TopTree_Study_DataRunD_bbEl.root")
         
        nevents=chTF.GetEntries()
        
        ii=0
    
        for iev in chTF:
        
            # define shorter variable 
            if isMuMu:
                d0mu=abs(iev.d0BeamSpot_muon[0])
                muond0.Fill(d0mu,1)
                muond0Wide.Fill(d0mu)
    
                # some debug statement
                if (0.01 < d0mu) and (d0mu < 0.02) and debug :
                        print "d0 is " , d0mu
                        print "weight is ", weight
                        print "PileUpWeight is ", PileUpWeight
                        print "LeptonWeight is ", LeptonWeight
    
            if isElEl :
                d0el=abs(iev.d0BeamSpot_electron[0])
                electrond0.Fill(d0el,1)
                electrond0Wide.Fill(d0el)
        
        
    
    # end of loop
    
    
    # write the histo in the output file
    outfile.cd()
    electrond0.Write()
    electrond0Wide.Write()
    muond0.Write()
    muond0Wide.Write()



outfile.Close()


