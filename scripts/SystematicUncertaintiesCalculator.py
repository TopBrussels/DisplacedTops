##############
# pyroot macro to calcultate the Yield for different Signal region and create a compilable tex file.
# February 2016 by qpython@cern.ch
#############


import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt



# list of systematics
systematicUncertainties=["sf_id_electron","sf_reco_electron","sf_iso_muon","sf_id_muon","evt_puSF"]

# create an dictionary that will connect the name of the systematics to the correct reweighting
my_dict = {}

# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"



# using Lorentz Vectors (lv) as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()

#channel
#channels=["_MuMu"]
channels=["_ElEl"]
#channels=["_ElEl","_MuMu"]

# path to tree
folderName="Systematics_29_8_2016"
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

# debug
debug=False

#
FastRun=True

treeName="tree"


# define the bound of the Signal region
#LowBounds=[0.001,0.02,0.05]
bound1 = 0.01

# template histogram that contains one single bin
hist=rt.TH1D("template","template",1,0,1,)

lumivalue = 3


# loop over the channel (lepton in final statue)
for chan in channels:
    

    isElEl=False
    isMuMu=False
    doubleArray=[]

    # get the xmlfile that corresponds to the channel
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

    SystType="sf_reco_electron"

    # create one root file per systType.
    outfile = rt.TFile("rootFiles/"+"Systematics/"+SystType+".root",'RECREATE')


    # loop over datasets
    for d in datasets:
#        if d.attrib['add'] == '1' :
#        if d.attrib['add'] == '1' and "QCD_" in str(d.attrib['name']):
        if d.attrib['add'] == '1' and "DYJets" not in str(d.attrib['name']) and "QCD" not in str(d.attrib['name']): 
#            print "found dataset to be added..." + str(d.attrib['name'])
            datasetNames.append(str(d.attrib['name']))
            print str(d.attrib['name'])
            # one array per dataset [name, title, Eqlumi, N1, N2, N3, SR1, SR2, SR3]
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


            ch.Add(pathTrunc+folderName+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")


            # you can also use two bins, one for the weighted event and one for the unweighted event

            
            # get number of events
            nevents=ch.GetEntries()
            
#            if nevents == 0 :
#                continue
            
    
            # calculate weight
            if isData :
                lumivalue=float(d.attrib['EqLumi'])
                
#            weight= lumivalue / float(d.attrib['EqLumi'])
            if False:
                print "lumivalue is " ,lumivalue
                print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                print "weight is " , weight

            # loop over the systematics shif:
            for systShift in ["down","up"]:
                
                # one histo per sample X systtype
                N1=rt.TH1D(sampleName+"_"+systShift,"N1",1,0,1)

                
                Yield=0
                ii=0
                # start of loop over events
                for iev in ch:
    
                    #                if ii % (nevents/50.) ==0 :
                    #                print  d.attrib['title']," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
                    ii+=1
                    if (FastRun and ii > 100):
                        continue
                        
    
                    if isMuMu:
                        PileUpWeight=iev.evt_puSF
                    if isElEl:
                        PileUpWeight=iev.evt_puSF
                    
                    LeptonWeight = 1
                    
        
                    # loop over the 2 highest pt letpon
                    for ilept in range (0,2):
    #                    x=str(ilept)
    #                    my_dict[x] = iev.sf_muon[ilept]
    #                    print my_dict
    
    #                    print "my_dict[x] is ", my_dict[x]
    
                        
    
                        # make the logic for the muon
                        if isMuMu:
                            sf_iso_muon = iev.sf_iso_up_muon[ilept] 
                            #print "sf_iso_muon is " ,sf_iso_muon
                            sf_id_muon = iev.sf_id_muon[ilept]
                            #print "sf_id_muon is ", sf_id_muon
    
                            LeptonWeight *= sf_iso_muon * sf_id_muon
    
                            if abs(iev.d0BeamSpot_muon[ilept]) < bound1:
                                passed1=False
                                if (debug):
                                    print "Electron and muon entering N1"
                                    print "d0 muon is " , iev.d0BeamSpot_muon[ilept]
                                
    
                        # eo the logic for the muon 
    
    
                                    
                        # make the logic for the electron
                        if isElEl :
    
                            sf_reco_electron_dict={'down':iev.sf_reco_down_electron[ilept], 'up':iev.sf_reco_up_electron[ilept]}
    
                            sf_reco_electron = sf_reco_electron_dict[systShift] 
                            #print "sf_reco_electron is " ,sf_reco_electron
                            sf_id_electron = iev.sf_id_electron[ilept]
                            #print "sf_id_electron is ", sf_id_electron
                            LeptonWeight *= sf_reco_electron * sf_id_electron
                            #for bound in range(0,len(SRxBounds)): 
    
    
                            # if one of the leptons  is smaller than bound, the event fails
                            if abs(iev.d0BeamSpot_electron[ilept]) < bound1:
                                passed1=False
                                if (debug):
                                    print "Electron and muon entering N1"
                                    print "d0 electron is " , iev.d0BeamSpot_electron[ilept]
    
                        # eo the logic for the electron
    
    
    
                    # Filling the histogram
                    if isData :
                        PileUpWeight=1
                        LeptonWeight=1
                                
    
                    # print info
    #                print "weight  is ", weight 
#                    print " PileUpWeight is ", PileUpWeight
#                    print "LeptonWeight is ", LeptonWeight
                    
    

                    Yield = Yield + PileUpWeight*LeptonWeight
#                    print Yield
    
                    N1.Fill(0.5,PileUpWeight*LeptonWeight)
                        
                
                # eo event loop
                outfile.cd()
                N1.Write()
                print "sampleName is " , sampleName
                print "systShift is ", systShift
                print "Yield is ", Yield
            # end of systematics shift
            idataset=idataset+1    
    
    


    # end of sample loop


# end of the channel loop

