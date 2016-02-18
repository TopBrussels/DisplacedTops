import os, sys
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain("doubleMuTree","doubleMuTree")
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
#pathToRootFile="/user/qpython/TopBrussels7X/CMSSW_7_2_1_patch1/src/TopBrussels/HToZZBachelorProjectNtupleMaker/FreyasNtuple"
pathToRootFile="../MergedTrees/12_2_2016/_MuMu/"
ch.Add(pathToRootFile+"DisplacedTop_Run2_TopTree_Study_NP_overlay_stopTobl_m500_Ctau10_MuMu.root")
# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

# book some histograms
h_mupt = ROOT.TH1F("h_mupt","muon p_{T}",100,0,500)
h_mueta = ROOT.TH1F("h_mueta","muon #eta",100,-4,4)
h_mud0 = ROOT.TH1F("h_mud0","muon |d_{0}|",200,0,2);
h_elept = ROOT.TH1F("h_elept","ele p_{T}",100,0,500)
h_eleeta = ROOT.TH1F("h_eleeta","ele #eta",100,-4,4)
h_eled0 = ROOT.TH1F("h_eled0","ele |d_{0}|",200,0,2);


# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()

# for bookkeeping
ii=0
nevents=ch.GetEntries()

# start of loop over events
for iev in ch:
    if ii % 10000 ==0 :
        print ii, "/", nevents
    ii+=1
# comment out the following lines if you are testing, it stops after a certain number of events
#    if ii==10000 :
#        break


# loop over muons - fill in lorentz vector and fill some histograms
    for imu in range(0,iev.nMuons_mumu) :
#        lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
        
        # Fill with the value from the tree! (open root file an search!)
        h_mupt.Fill(iev.pt_muon_mumu[imu])
#        h_mueta.Fill(lvmu.Eta())
#        h_mud0.Fill(abs(iev.d0_muon[imu]))

# loop over electrons - fill in lorentz vector and fill some histograms
    for iele in range(0,iev.nElectrons_mumu) :
        

#        lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
        h_eleeta.Fill(3)
#        h_eleeta.Fill(lve.Eta())
#        h_elept.Fill(lve.Pt())
#        h_eled0.Fill(abs(iev.d0_electron[iele]))

# end of loop


# create canvas
t3=ROOT.TCanvas()

# create sub-pads and cd() to them, draw some histograms
t3.Divide(3,2)
t3.cd(1)
h_mupt.Draw()
t3.cd(2)
h_mueta.Draw()
t3.cd(3)
ROOT.gPad.SetLogy()
h_mud0.Draw()
t3.cd(4)
h_elept.Draw()
t3.cd(5)
h_eleeta.Draw()
t3.cd(6)
ROOT.gPad.SetLogy()
h_eled0.Draw()
t3.Update()
t3.Print("myplot.gif") # TCanvas::Print also can make .pdf or .root/.C plots
