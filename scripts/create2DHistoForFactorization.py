############
#Pyroot script that reads some tree and fill d0 vs d0 2D histograms. 
#These 2D histograms are used to extend the non-zero prediction of the given background in the SR using the product of the efficiency to pass the d0 cut of the first lepton with the efficiency to pass the d0 cut of the second lepton
# April 2016 by qpython@cern.ch 

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
date="NoDisplacedTriggerNoBlinding"

# array with composite dataset and matching string
dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "DrellYann","Signal"]
compositeDatasets= ["WJets", "Diboson", "SingleTop", "TTJets", "DrellYann","stopTobl_m500_Ctau10"] # title in the xml config


# verbosity
debug = False

# fast run
fastRun = False



electrond0VsElectronsd0Sum=rt.TH2D("electrond0VsElectronsd0Sum","electrond0VsElectronsd0", 50, 0.0, 0.05, 50, 0.0, 0.05)
muond0VsMuond0Sum=rt.TH2D("muond0VsMuond0Sum","muond0VsMuond0",50, 0.0, 0.05, 50, 0.0, 0.05)
muond0VsElectrond0Sum=rt.TH2D("muond0VsElectrond0Sum","muond0VsElectrond0",50, 0.0, 0.05, 50, 0.0, 0.05)


# remove low d0 part of the histo
#muond0VsMuond0=rt.TH2D("muond0VsMuond0","muond0VsMuond0",35, 0.0015, 0.05, 35, 0.015, 0.05)
#muond0VsElectrond0=rt.TH2D("muond0VsElectrond0","muond0VsElectrond0",35, 0.0015, 0.05, 35, 0.0015, 0.05)


# name of the tree in the root file
treeName="tree"

# loop over the different composite data set
i_comp=0
for compositeDataset in compositeDatasets:
    print "\n", "compositeDataset is " , dataSetTitles[i_comp] , ":"


    # define d0 histograms, one per composite dataset
    electrond0VsElectronsd0=rt.TH2D("electrond0VsElectronsd0"+dataSetTitles[i_comp],"electrond0VsElectronsd0", 50, 0.0, 0.05, 50, 0.0, 0.05)
    muond0VsMuond0=rt.TH2D("muond0VsMuond0"+dataSetTitles[i_comp],"muond0VsMuond0",50, 0.0, 0.05, 50, 0.0, 0.05)
    muond0VsElectrond0=rt.TH2D("muond0VsElectrond0"+dataSetTitles[i_comp],"muond0VsElectrond0",50, 0.0, 0.05, 50, 0.0, 0.05)
          
    
    FilterString=compositeDatasets[i_comp]
    outfile = rt.TFile("rootFiles/"+dataSetTitles[i_comp]+"2D.root",'RECREATE')

    # loop over the different channels
    for chan in channels :
        isElEl =False
        isMuMu =False
    
    
        if chan == "_ElEl" :
            isElEl = True
        if chan == "_MuMu" :
            isMuMu = True
    
        if isElEl:
            tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')    
        if isMuMu:
            tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')


        root =  tree.getroot()
        datasets = root.find('datasets')
        if (debug):
            print "found  "  + str(len(datasets)) + " datasets"
        
        # reset the list of dataset
        datasetNames = []


        # loop over the dataset inside the composite dataset (name)
        i_dataset=0
        for d in datasets:

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

                if "Data" not in sampleName and "NP" not in sampleName:
                    isBgMC = True
                if "NP" in sampleName:
                    isSignal = True
                if "Data" in sampleName:
                    isData = True
                

                ch.Add(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")

        
                # get the lumi weight 
                weight= lumivalue / float(d.attrib['EqLumi'])
                if (debug):
                    print "lumivalue is " ,lumivalue
                    print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                    print "weight is " , weight


                # start of loop over events 
                ii=0
                for iev in ch:
                    
                    ii=ii+1
                    if fastRun and 100 < ii :
                        continue

#                    # just pick few points for signal
#                    if isSignal and 900 < ii:
#                        continue
                    
                    
                    # PU weight
                    if isMuMu:
                        PileUpWeight=iev.evt_puSF
                    if isElEl:
                        PileUpWeight=iev.evt_puSF
        

                    # Lepton Weight
                    LeptonWeight=1.0
                    for ilept in range (0,2):
                        if isMuMu:
                            LeptonWeight *= iev.sf_muon[0]
                        if isElEl:
                            LeptonWeight *= iev.sf_electron[0]
                    
        
                    # Filling corresponding histograms
                    if isMuMu:
                        d01mu=abs(iev.d0BeamSpot_muon[0])
                        d02mu=abs(iev.d0BeamSpot_muon[1])
                        muond0VsMuond0.Fill(d01mu,d02mu,weight*PileUpWeight*LeptonWeight)
                        if (debug):
                            print "d0 is " , d0mu
                            print "weight is ", weight
                            print "PileUpWeight is ", PileUpWeight
                            print "LeptonWeight is ", LeptonWeight
                    if isElEl :
                        d01el=abs(iev.d0BeamSpot_electron[0])
                        d02el=abs(iev.d0BeamSpot_electron[1])
                        electrond0VsElectronsd0.Fill(d01el,d02el,weight*PileUpWeight*LeptonWeight)
        
                    # eo loop over the event

            # eo loop over the dataset

        #eo loop over the channnel

    # write the two histo and close the file
    outfile.cd()
    muond0VsMuond0.Write()
    electrond0VsElectronsd0.Write()
    outfile.Close()

    

    i_comp=i_comp+1
    # end of loop over the comp dataset
