############
#Pyroot script that calculates trigger SF.
# A data sample is used using MET trigger and is compare with ttbar MC.
# August 2016 by qpython@cern.ch 

# basic import
import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt


#channel 
#channels=["_MuMu"] 
#channels=["_ElEl"] 
channels=["_ElEl","_MuMu"]


#root file  
#base of the path to the root file
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

#date
directory="CMSSW76V4_TTLetp_10_8_2016"
# array with composite dataset and matching string

compositeDatasets= ["WJets", "Diboson", "SingleTop", "TTJets", "DrellYann","stopTobl_m500_Ctau10"] # title in the xml config


# verbosity
debug = False

# fast run
fastRun = False


# SF vs electron pt
electron1Pt=rt.TH1D("electron1Pt","electron1Pt", 8, 0.0, 200)
electron2Pt=electron1Pt.Clone("electron2Pt")


# SF vs electron eta
electron1Eta=rt.TH1D("electron1Eta","electron1Eta", 6, 0.0, 2.4)
electron2Eta=electron1Eta.Clone("electron2Eta")


# 
electronPtVsElectronPt=rt.TH2D("electronPtVsElectronPt","electronPtVsElectronPt", 8, 0.0, 200, 8, 0.0, 200)
electronEtaVsElectronEta=rt.TH2D("electronEtaVsElectronEta","electronEtaVsElectronEta", 6, 0.0, 2.4, 6, 0.0, 2.4)


# clone
MuonPtVsMuonPt=electronPtVsElectronPt.Clone("MuonPtVsMuonPt")
MuonEtaVsMuonEta=electronEtaVsElectronEta.Clone("MuonEtaVsMuonEta")


# remove low d0 part of the histo
#muond0VsMuond0=rt.TH2D("muond0VsMuond0","muond0VsMuond0",35, 0.0015, 0.05, 35, 0.015, 0.05)
#muond0VsElectrond0=rt.TH2D("muond0VsElectrond0","muond0VsElectrond0",35, 0.0015, 0.05, 35, 0.0015, 0.05)


# name of the tree in the root file
treeName="tree"
ch = rt.TChain(treeName,treeName)

outfile = rt.TFile("rootFiles/TriggerSF.root",'RECREATE')

# loop over the different channels
for chan in channels :
    isElEl =False
    isMuMu =False
    
    sampleName="DataRunD"
    
    if chan == "_ElEl" :
        isElEl = True
    if chan == "_MuMu" :
        isMuMu = True
            
    if isElEl:
        tree = ET.ElementTree(file='../config/Yield'+chan+'V4.xml')    
    if isMuMu:
        tree = ET.ElementTree(file='../config/Yield'+chan+'V4.xml')


    root =  tree.getroot()
    datasets = root.find('datasets')
    if (debug):
        print "found  "  + str(len(datasets)) + " datasets"
        
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
                

    ch.Add(pathTrunc+directory+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")

        

    # start of loop over events 
    ii=0
    for iev in ch:
        ii=ii+1
        if fastRun and 100 < ii :
            continue

        
    # Filling corresponding histograms
        if isElEl :
            pt1lept=iev.pt_electron[0]
            pt2lept=iev.pt_electron[1]
            electron1Pt.Fill(pt1lept)
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
     
    # eo loop over the event


#eo loop over the channnel

# write the two histo and close the file
outfile.cd()
electron1Pt.Write()
#muond0VsMuond0.Write()
outfile.Close()

