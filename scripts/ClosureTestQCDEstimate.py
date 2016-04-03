##############
# Pyroot macro to calcultate the Yield for different Signal region and create a compilable tex file.
# February 2016 by qpython@cern.ch
#############


import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt



# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"



# using Lorentz Vectors (lv) as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()

#channel
channels=["_MuMu"]
#channels=["_ElEl"]
#channels=["_ElEl","_MuMu"]

# path to tree
date="1_4_2016"
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

# debug
debug=False

# define various bound on the first and second leptonn
boundsLept1 = [0.012, 0.015, 0.018]
boundsLept2 = [0.012, 0.015, 0.018]
#boundLept2 = 0.015


# template histogram that contains one single bin
hist=rt.TH1D("template","template",1,0,1,)

lumivalue = 3


# loop over the bound of the first lepton
for boundLept1 in boundsLept1 :

    # loop over the bound of the first lepton
    for boundLept2 in boundsLept2 :
    
        print "boundLept1 is" , boundLept1
        print "boundLept2 is" , boundLept2
        
        # loop over the channel (lepton in final statue)
        for chan in channels:
        
            # loop over the low bounds
        #    for ilb in range (1,len(LowBounds)+1):
        #        dict = ('bgMCSum+')
        #        bgMCSum
            
            # Histogram containing the number of events
            NonQCDBase=rt.TH1D("NonQCDBase"+chan,"NonQCDBase",1,0,1)
            DataBase=rt.TH1D("DataBase"+chan,"DataBase",1,0,1)
            NonQCDTarget=rt.TH1D("NonQCDTarget"+chan,"NonQCDTarget",1,0,1)
            DataTarget=rt.TH1D("DataTarget"+chan,"DataTarget",1,0,1)
            EstimatedQCDTarget=rt.TH1D("EstimatedQCDTarget"+chan,"EstimatedQCDTarget",1,0,1)
        
        
        
            isElEl=False
            isMuMu=False
            doubleArray=[]
        
            # get the xmlfile that corresponds to the channel
            if "MuMu" in chan:
                isMuMu=True
                tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
                treeName="doubleMuTree"
                FinalState="At least two muons"
                print FinalState
            elif "ElEl" in chan:
                isElEl=True
                tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')
                treeName="doubleElTree"
                FinalState="At least two electrons"
                print FinalState
            elif "ElMu" in chan:
                tree = ET.ElementTree(file='../config/Yield_FullSamplesElMuV0.xml')
            else:
                print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
                sys.exit()
        
                
        
            root = tree.getroot()
            datasets = root.find('datasets')
            print "found  "  + str(len(datasets)) + " datasets"
            datasetNames = []
            idataset=0
        
            # index to get the bin conten of the 3 different signal region
        #    iSR1=10
        #    iSR2=12
        #    iSR3=14
        
            # loop over datasets
            for d in datasets:
                if d.attrib['add'] == '1' :
        #        if d.attrib['add'] == '1' and "QCD_" in str(d.attrib['name']):
        #            print "found dataset to be added..." + str(d.attrib['name'])
                    datasetNames.append(str(d.attrib['name']))
                    print str(d.attrib['name'])
                    # one array per dataset [name, title, Eqlumi, N1, N2, N3, SR1, SR2, SR3]
                    ch = rt.TChain(treeName,treeName)
                    sampleName=d.attrib['name']
                    
                    # fix the type of dataset (bgMC, signal or data)
                    isBgMC = False
                    isQCD = False 
                    isData = False
                    
                    if "Data" not in sampleName and "NP" not in sampleName:
                        isBgMC = True
                    if "QCD" in sampleName:
                        isQCD = True
                    if "Data" in sampleName:
                        isData = True
        
        
                    ch.Add(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")
                    Sum_SR1=Sum_SR2=Sum_SR3=0
                    
                    
                    # get number of events
                    nevents=ch.GetEntries()
                    
        #            if nevents == 0 :
        #                continue
                    
            
                    # calculate weight
                    if isData :
                        lumivalue=float(d.attrib['EqLumi'])
                        
                    weight= lumivalue / float(d.attrib['EqLumi'])
                    if (1):
                        print "lumivalue is " ,lumivalue
                        print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                        print "weight is " , weight
                        
                    ii=0
                    # start of loop over events
                    for iev in ch:
        
                        if isMuMu:
                            PileUpWeight=iev.evt_puSF_mumu
                        if isElEl:
                            PileUpWeight=iev.evt_puSF_elel
                        
                        LeptonWeight=1.0
        
                        
        #                if ii % (nevents/50.) ==0 :
        
                        #  skip 99% of the events just to run faster                             
        #                if not ii % 100  == 0:                                              
        #                    continue    
        #                print  d.attrib['title']," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
        
                        ii+=1
                        isInBase= False
                        isInTarget= False
            
                        # get the scale factor
                        for ilept in range (0,2):
                            # bo the logic for the muon
                            if isMuMu:
                                LeptonWeight *= iev.sf_muon_mumu[ilept]                        
        
        
                        # define shorter variable
                        d01=abs(iev.d0BeamSpot_muon_mumu[0])
                        d02=abs(iev.d0BeamSpot_muon_mumu[1])
                            
                        # event in base if both lepton are smaller than a bound but still in DCR
                        if 0.01 < d01 and  d01 < boundLept1 and 0.01 < d02 and d02 < boundLept2:
                            isInBase=True
                            if (debug):
        #                    if (True):
                                print "Leptons entering Base"
                                print "d0 muon[0] is " , iev.d0BeamSpot_muon_mumu[0]
                                print "d0 muon[1] is " , iev.d0BeamSpot_muon_mumu[1]
                                    
        
                        # event in base if both lepton are bigger than previous bound but still in DCR
                        if boundLept1 < d01 and d01 < 0.02 and boundLept2 < d02 and d02 < 0.02:
                            isInTarget=True
                            if (debug):
                                print "Lepton entering Target"
                                print "d0 muon is " , iev.d0BeamSpot_muon_mumu[ilept]
                        # eo the logic for the muon 
        
        
        
                        # Filling the histogram
                        if (isData) :
                            PileUpWeight=1
                            LeptonWeight=1
                                    
                        # filling base histo
                        if (isInBase):
        #                    print "isInbase!!"
                            if (isData):
                                DataBase.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                            if (isBgMC and not isQCD):
                                NonQCDBase.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                                
                        # filling target histo
                        elif (isInTarget):
                            if (isData):
                                DataTarget.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                            if (isBgMC and not isQCD):
                                NonQCDTarget.Fill(0.5,weight*PileUpWeight*LeptonWeight)
    #                    else :
    #                        print "Neither in target neither in base!!"
        
                        
                    
                    # eo event loop
        
            # eo dataset loop
            
    
            EstimatedQCDTarget.Fill(0.5,1000)
        
        
            # make some printout
            print  "NonQCDBase is ", NonQCDBase.GetBinContent(1) 
            print  "DataBase is ", DataBase.GetBinContent(1)  
            print  "NonQCDTarget is ", NonQCDTarget.GetBinContent(1)  
            print  "DataTarget is ", DataTarget.GetBinContent(1) 
            print  "EstimatedQCDTarget is ", EstimatedQCDTarget.GetBinContent(1) 
                                
            
                    # Fill the two D array for clearer output
                    
                    
        
        """
        
                    datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], nevents, NonQCDBase.GetBinContent(1),NonQCDBase.GetBinError(1), DataBase.GetBinContent(1), DataBase.GetBinError(1), NonQCDTarget.GetBinContent(1),NonQCDTarget.GetBinError(1), DataTarget.GetBinContent(1),DataTarget.GetBinError(1), EstimatedQCDTarget.GetBinContent(1),EstimatedQCDTarget.GetBinError(1)]
                    print datasetArray
                    
                    
        #
                    doubleArray.append(datasetArray)
                    print "double array is"
                    print doubleArray
        
                
                    # end of event loop
                    idataset=idataset+1    
            
            
            # print the summary contained it the double array
            for i in range (0,len(doubleArray)):
        #    for i in range (0,lendidataset):
                print "---------"
                print "NEW SAMPLE!!!"
                print "---------"
                for j in range (0,len(doubleArray[i])):
                    print doubleArray[i][j]
                print ""
                print ""
            
            # print the minimum necessary for the Yield
            for i in range (0,len(doubleArray)):
        #    for i in range (0,idataset):
                print "---------"
                print "Sample is ", doubleArray[i][1], "and the yield in the signal regions are:"
                print "---------"
                for j in range (5,len(datasetArray)):
                    print doubleArray[i][j]
                print ""
                print ""
    
        """
    
print "end of the program !!!!"

