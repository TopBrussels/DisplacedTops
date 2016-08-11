############
#Pyroot script that calculates trigger SF.
# A data sample is used using MET trigger and is compare with ttbar MC.
# August 2016 by qpython@cern.ch 

# basic import
import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt


# list of channels
channels=["ttElEl","ttMuMu"]
#channels=["ttMuMu"] 
#channels=["ttElEl"] 

# list of samples
sampleNames=["DataRunD","TTJets_Dilept"]

# list of indices
indices=["0","1"]

#base of the path to the root file
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

#directory
directory="CMSSW76V4_TTLetp_10_8_2016"

# verbosity
debug = False

# fast run
fastRun = False

# name of the tree in the root file
treeName="tree"

# root output
outfile = rt.TFile("rootFiles/TriggerSF.root",'RECREATE')

# loop over the different channels
for chan in channels :
    print "chan is ", chan

    isttElEl =False
    isttMuMu =False
    
    if chan == "ttElEl" :
        isttElEl = True
        lepton="Electron"
    if chan == "ttMuMu" :
        isttMuMu = True
        lepton="Muon"
            
    if isttElEl:
        tree = ET.ElementTree(file='../config/ttLeptonsV4.xml')    
    if isttMuMu:
        tree = ET.ElementTree(file='../config/ttLeptonsV4.xml')

    root =  tree.getroot()
    datasets = root.find('datasets')
    if (debug):
        print "found  "  + str(len(datasets)) + " datasets"
        

    # TTBar or Data
    for sampleName in sampleNames:
        print "sampleName is ", sampleName

        ch = rt.TChain(treeName,treeName)
        ch.Add(pathTrunc+directory+"/_"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+"_"+chan+".root")

        
        for index in indices :
            print "index is ", index
    
            
            # string of the histo depending on the lepton (channel), the index and the samples name
            histStr=lepton+index+sampleName
    
            # tot histo
            totPtLeptonIndexSample=rt.TH1D("totPt"+histStr,"totPt"+histStr, 8, 0.0, 200)
            totEtaLeptonIndexSample=rt.TH1D("totEta"+histStr,"totEta"+histStr, 6, 0.0, 2.4)
            
            # pass histo clone from tot to ensure bin consistency
            passPtLeptonIndexSample=totPtLeptonIndexSample.Clone("passPt"+histStr)
            passEtaLeptonIndexSample=totEtaLeptonIndexSample.Clone("passEta"+histStr)
    
    
            
    
            #TH2D 
            #passElectronPtVsElectronPt=rt.TH2D("passElectronPtVsElectronPt","passElectronPtVsElectronPt", 8, 0.0, 200, 8, 0.0, 200)
            #passElectronEtaVsElectronEta=rt.TH2D("passElectronEtaVsElectronEta","passElectronEtaVsElectronEta", 6, 0.0, 2.4, 6, 0.0, 2.4)
    
    
            #Clone 
            #MuonPtVsMuonPt=passElectronPtVsElectronPt.Clone("MuonPtVsMuonPt")
            # MuonEtaVsMuonEta=passElectronEtaVsElectronEta.Clone("MuonEtaVsMuonEta")
    
    
                
        
            # start of loop over events 
            ii=0
            for iev in ch:
                ii=ii+1
                if fastRun and 100 < ii :
                    continue
        
                
            # Filling corresponding histograms
                if isttElEl :
    
                    # define all the variables
                    if index=="0":
                        ptIndexlept=iev.pt_electron[0]
                        etaIndexlept=iev.eta_electron[0]
                    
                    if index=="1":
                        ptIndexlept=iev.pt_electron[1]
                        etaIndexlept=iev.eta_electron[1]
    
                if isttMuMu :
                    if index=="0":
                        ptIndexlept=iev.pt_muon[0]
                        etaIndexlept=iev.eta_muon[0]
    
                    if index=="1":
                        ptIndexlept=iev.pt_muon[1]
                        etaIndexlept=iev.eta_muon[1]
        
        
                # fill tot histo
                totPtLeptonIndexSample.Fill(ptIndexlept)
                totEtaLeptonIndexSample.Fill(etaIndexlept)
                
                # fill pass histo
                if (iev.crossTrigged):
                    passPtLeptonIndexSample.Fill(ptIndexlept)
                    passEtaLeptonIndexSample.Fill(etaIndexlept)
    
            # eo loop over the event
    



            # make the efficiency with the properly filled histograms
            effPt=rt.TEfficiency(passPtLeptonIndexSample,totPtLeptonIndexSample)
#            effPt.Draw()
    
            effEta=rt.TEfficiency(passEtaLeptonIndexSample,totEtaLeptonIndexSample)
#            effEta.Draw()
    

            # write the histo on the output file
            outfile.cd()
            totPtLeptonIndexSample.Write()
            passPtLeptonIndexSample.Write()
            effPt.Write()
    
            totEtaLeptonIndexSample.Write()
            passEtaLeptonIndexSample.Write()
            effEta.Write()

        # eo loop over the indices

    
    #eo loop over the sampleName
    

#eo loop over the channnel


outfile.Close()





"""
                # get the lumi from the Data
                if d.attrib['add'] == '1' and "Data" in str(d.attrib['name']):
                    lumivalue=float(d.attrib['EqLumi'])
    
                # check if in the right composite dataset
                if d.attrib['add'] == '1' and FilterString ==  str(d.attrib['title']): 
                    datasetNames.append(str(d.attrib['name']))
                    print str(d.attrib['name'])
                    ch = rt.TChain(treeName,treeName)
                    sampleName=d.attrib['name']
    
    
                    # fix the type of dataset (bgMC, signal or data)                           
                    isBgMC = False
                    isSignal = False
                    isData = False
"""
                    

"""
                
            if isMuMu:
                pt1lept=iev.pt_muon[0]
                pt2lept=iev.pt_muon[1]
    
    #            muond0VsMuond0.Fill(d01mu,d02mu,weight*PileUpWeight*LeptonWeight)
                if (debug):
                    print "d0 is " , d0mu
                    print "weight is ", weight
                    print "PileUpWeight is ", PileUpWeight
                    print "LeptonWeight is ", LeptonWeight
"""   


    
