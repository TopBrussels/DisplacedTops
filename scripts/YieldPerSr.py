import os, sys
from array import array
import ROOT as rt
import tdrstyle, CMS_lumi
##############
# example pyroot loop for yield counting on output trees of Ntupler
# March 2015 by qpython@cern.ch
#

# using Lorentz Vectors (lv) as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()


# for bookkeeping
lumivalue = 1000.0


#samplesnamesfancy = ["Drell-Yan","t#bar{t}+jets","W+jets","W+jets","QCD","QCD","QCD","QCD","QCD","QCD","QCD","QCD","QCD"]
samples = ["dy","ttbar","wjetsplus","wjetsminus","qcd1","qcd2","qcd3","qcd4","qcd5","qcd6","qcd7","qcd8","qcd9"]
samplesxsecs = [2008.4,831.76,11811.4,8677.3,677300000.*0.007,866600000.*0.00044,164300000.*0.00816,21810000.*0.01522,2999000.*0.02424,3529000.*0.158,128500.*0.0406,185900000.*0.00272,3495000.*0.01255]
samplespresel = [2820473.0,25437856.0,699606.0,226439.0,1986513.0,4768929.0,3742583.0,3893676.0, 3468633.,1958930.0,999553.0,1730223.0,1999717.0,1000000.0]
weight=[0,0,0,0,0,0,0,0,0,0,0,0,0] # computable from other array
pathToRootFile="/user/qpython/TopBrussels7X/CMSSW_7_2_1_patch1/src/TopBrussels/HToZZBachelorProjectNtupleMaker/FreyasNtuple"
samplesrootfiles=[pathToRootFile+"/DYJetsToLL_M-50_13TeV-madgraph-pythia8*.root",
                  pathToRootFile+"/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola*.root",
                  pathToRootFile+"/WplusToMuNu_CT10_13TeV-powheg-pythia8*.root",
                  pathToRootFile+"/WminusToMuNu_CT10_13TeV-powheg-pythia8*.root",
                  pathToRootFile+"/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_170toInf_bcToE_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8*.root"]

#Define array containing the number of event passing d_0 cut on electron and muon
N1=[0,0,0,0,0,0,0,0,0,0,0,0,0] # muon and electron with abs(d0) > 0.02
N2=[0,0,0,0,0,0,0,0,0,0,0,0,0] # muon and electron with abs(d0) > 0.05
N3=[0,0,0,0,0,0,0,0,0,0,0,0,0] # muon and electron with abs(d0) > 0.1
SR1=[0,0,0,0,0,0,0,0,0,0,0,0,0] # Exclusive singal region SR1=N1-N2
SR2=[0,0,0,0,0,0,0,0,0,0,0,0,0] # Exclusive singal region SR2=N2-N3
SR3=[0,0,0,0,0,0,0,0,0,0,0,0,0] # Exclusive singal region SR3=N3
Sum_SR1=0
Sum_SR2=0
Sum_SR3=0

#for isam in range (0,1) :
for isam in range(0,len(samples)) :
    # the analysis structure see TTree/TChain description on rt.cern.ch
    ch = rt.TChain("tree","tree")
    # TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
    ch.Add(samplesrootfiles[isam])
# very loud but useful to know what variables are stored in a tree... it prints them all
#    ch.Print()


    nevents=ch.GetEntries()
    if nevents == 0 :
        continue


    weight[isam]=lumivalue*samplesxsecs[isam]/(samplespresel[isam])
    ii=0
    # start of loop over events
    for iev in ch:
        if ii % (nevents/50.) ==0 :
            print samples[isam]," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
        ii+=1


    # loop over muons - fill in lorentz vector and fill some histograms
        for imu in range(0,iev.nMuons) :
        
            for iele in range(0,iev.nElectrons) :
                if iev.charge_muon[imu]*iev.charge_electron[iele]>0 :
                    continue

                lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
                lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
                
                # Define the displaced regions
                if abs(iev.d0_electron[iele])>0.02 and abs(iev.d0_muon[imu])>0.02 :
                    print "Electron and muon entering N1"
                    print "d0 electron is " , iev.d0_electron[iele]
                    print "d0 muon is " , iev.d0_muon[imu]
                    N1[isam]=N1[isam]+1*weight[isam]
                    print N1
                if abs(iev.d0_electron[iele])>0.05 and abs(iev.d0_muon[imu])>0.05 :
                    print "Electron and muon entering N2"
                    print "d0 electron is " , iev.d0_electron[iele]
                    print "d0 muon is " , iev.d0_muon[imu]
                    N2[isam]=N2[isam]+1*weight[isam]
                    print N2
                if abs(iev.d0_electron[iele])>0.1 and abs(iev.d0_muon[imu])>0.1 :
                    print "Electron and muon entering N3"
                    print "d0 electron is " , iev.d0_electron[iele]
                    print "d0 muon is " , iev.d0_muon[imu]
                    N3[isam]=N3[isam]+1*weight[isam]
                    print N3

    # define the non-overlaping singal region (SR) out of the displaced regions
    SR1[isam]=N1[isam]-N2[isam]    
    SR2[isam]=N2[isam]-N3[isam]
    SR3[isam]=N3[isam]

    # compute the sum of background in the 3 SR
    Sum_SR1=Sum_SR1+SR1[isam]
    Sum_SR2=Sum_SR2+SR2[isam]
    Sum_SR3=Sum_SR3+SR3[isam]
    
    # end of event loop

print "the array of weight is", weight

# print the final number per SR    
print "we expect " , Sum_SR1 ,"events in the SR1"               
print "we expect " , Sum_SR2 ,"events in the SR2"               
print "we expect " , Sum_SR3 ,"events in the SR3"               
print "end of the program"


# end of sample loop



