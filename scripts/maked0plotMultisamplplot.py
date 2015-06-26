import os, sys
import ROOT
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

# the analysis structure see TTree/TChain description on root.cern.ch
ch = ROOT.TChain("tree","tree")
# TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
#pathToRootFile="/user/qpython/TopBrussels7X/CMSSW_7_2_1_patch1/src/TopBrussels/HToZZBachelorProjectNtupleMaker/FreyasNtuple"
pathToRootFile="../Craneen/"

samples = ["m500Ctau10","m800Ctau10","m1100Ctau10","m500Ctau100","m500Ctau1000","TTJets_tree.root","WToLNu_tree.root"]

samplesrootfiles=[
    pathToRootFile+"stopTobl_m500_Ctau10_tree.root",
    pathToRootFile+"stopTobl_m800_Ctau10_tree.root",
    pathToRootFile+"stopTobl_m1100_Ctau10_tree.root",
    pathToRootFile+"stopTobl_m500_Ctau100_tree.root",
    pathToRootFile+"stopTobl_m500_Ctau1000_tree.root",
    pathToRootFile+"TTJets_tree.root",
    pathToRootFile+"WToLNu_tree.root",
#    pathToRootFile+"stopTobl_m1100_Ctau10_tree.root",
    ]


# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()

# book some histograms
#"""

#h_eled0_bis = ROOT.TH1F("h_eled0_bis","ele |d_{0}|",200,0,2);
#"""

# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=ROOT.TLorentzVector()
lve=ROOT.TLorentzVector()

# for bookkeeping
ii=0


plotvector=[]


# loop over the samples
for isam in range(0,len(samples)) :

#    print "isam is " isam
    
    ch.Add(samplesrootfiles[isam])

    nevents=ch.GetEntries()

#    h_mud0_work = h_mud0.Clone("h_mud0_"+samples[isam])
#    h_eled0_work = h_eled0.Clone("h_eled0_"+samples[isam])

    h_mupt = ROOT.TH1F("h_mupt","muon p_{T}",100,0,500)
    h_mueta = ROOT.TH1F("h_mueta","muon #eta",100,-4,4)
    h_mud0 = ROOT.TH1F("h_mud0","muon |d_{0}|",200,0,2);
    h_elept = ROOT.TH1F("h_elept","ele p_{T}",100,0,500)
    h_eleeta = ROOT.TH1F("h_eleeta","ele #eta",100,-4,4)
    h_eled0 = ROOT.TH1F("h_eled0","ele |d_{0}|",200,0,2);



    # start of loop over events
    for iev in ch:
        if ii % 10000 ==0 :
            print ii, "/", nevents
        ii+=1
    # comment out the following lines if you are testing, it stops after a certain number of events
    #    if ii==10000 :
    #        break
    
    
    # loop over muons - fill in lorentz vector and fill some histograms
        for imu in range(0,iev.nMuons) :
            lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
            
            h_mupt.Fill(lvmu.Pt())
            h_mueta.Fill(lvmu.Eta())
            h_mud0.Fill(abs(iev.d0_muon[imu]))
    
    # loop over electrons - fill in lorentz vector and fill some histograms
        for iele in range(0,iev.nElectrons) :
    
            lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])


            h_eleeta.Fill(lve.Eta())
            h_elept.Fill(lve.Pt())
            h_eled0.Fill(abs(iev.d0_electron[iele]))
            




# end of loop

#    c1=ROOT.TCanvas()
#    c1.cd()
#    h_eled0.Draw()
#    h_eled0_bis.Draw("SAME")
#    c1.Update()
#    c1.Print("facoplot.gif")

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
    t3.Print("myplot_"+samples[isam]+".gif") # TCanvas::Print also can make .pdf or .root/.C plots
