import os, sys
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain("tree","tree")
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
ch.Add("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/*.root")
# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

# create canvas
t3=ROOT.TCanvas()


ch.Draw("fabs(d0_electron):fabs(d0_muon)","charge_muon*charge_electron<0","colz")
