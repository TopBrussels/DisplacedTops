import xml.etree.cElementTree as ET
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


#channels=["_ElEl","_MuMu"]
channels=["_MuMu"]

date="12_2_2016"

# double array containing an array of samples. Each sample is an array with different varaibles such as ["fancyname", "name" x-sec,samplepresels,...] 
N1=N2=N3=SR1=SR2=SR3=Sum_SR1=Sum_SR2=Sum_SR3=0
doubleArray=[
]

lumivalue = 2600


pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"



# loop over the channel (lepton in final statue)
for chan in channels:
    

    # get the xmlfile that corresponds to the channel
    if "MuMu" in chan:
        tree = ET.ElementTree(file='../config/TreeProc_FullSamplesMuMuV0.xml')
#        tree = ET.ElementTree(file='../config/FullSamplesMuMuV0.xml')
    elif "ElEl" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElElV0.xml')
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElMuV0.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the ri\
ght directories!!!"
        sys.exit()

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    idataset=0

    # loop over datasets
    for d in datasets:
        if d.attrib['add'] == '1': # and "ZZ" in str(d.attrib['title']):
            print "found dataset to be added..." + str(d.attrib['name'])
            datasetNames.append(str(d.attrib['name']))
            print str(d.attrib['name'])
            # one array per dataset [name, title, Eqlumi, N1, N2, N3, SR1, SR2, SR3]
            datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], 0, 0, 0, 0, 0, 0]
            print datasetArray
            # the analysis structure see TTree/TChain description on rt.cern.ch
            ch = rt.TChain("doubleMuTree","doubleMuTree")
#            TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
#            ch.Add(samplesrootfiles[d])
            sampleName=d.attrib['name']
            ch.Add("/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/12_2_2016/_MuMu/DisplacedTop_Run2_TopTree_Study_"+sampleName+"_MuMu.root")
            
            # very loud but useful to know what variables are stored in a tree... it prints them all
            #    ch.Print()
            
            # get number of events
            nevents=ch.GetEntries()
            if nevents == 0 :
                continue
    
            # calculate weight
#            weight=lumivalue / float(d.attrib['EqLumi'])
            weight=1
            ii=0
            # start of loop over events
            for iev in ch:
                if ii % (nevents/50.) ==0 :
                    print  d.attrib['title']," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
                    ii+=1
    
    
        # loop over muons - fill in lorentz vector and fill some histograms
                    for imu in range(0,iev.nMuons_mumu) :
                

                
#                    if iev.charge_muon[imu]*iev.charge_electron[iele]>0 : # skip event with same charge
#                        continue
#                    if iev.pfIso_electron[iele] > 0.1 : # skip non isolated electron
#                        print "isolation of the electron is " , iev.pfIso_electron[iele] , ". This event will be skiped"
#                        continue
#                    if iev.pfIso_muon[imu] > 0.1 : # skip non isolated muon
#                        print "isolation of the muon is " , iev.pfIso_muon[imu] , ". This event will be skiped"
#                        continue
    
#                    lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
#                    lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
                    
                    # Define the displaced regions
                #if abs(iev.d0_electron[iele])>0.02 and abs(iev.d0_muon[imu])>0.02 :
#                        if abs(iev.d0BeamSpot_muon_mumu[imu]) > 0.02:
                        
                        passed= False
                        if abs(iev.d0BeamSpot_muon_mumu[imu]) > 0.0000001:
                        
                            print "Electron and muon entering N1"
                    #print "d0 electron is " , iev.d0_electron[iele]
                            print "d0 muon is " , iev.d0BeamSpot_muon_mumu[imu]
                            passed=True
                            
                if (passed==True):
                    N1=N1+1*weight
                    print N1


#                if abs(iev.d0_electron[iele])>0.05 and abs(iev.d0_muon[imu])>0.05 :
#                    print "Electron and muon entering N2"
#                    print "d0 electron is " , iev.d0_electron[iele]
#                    print "d0 muon is " , iev.d0_muon[imu]
#                    N2[d]=N2[d]+1*weight
#                    print N2
#                if abs(iev.d0_electron[iele])>0.1 and abs(iev.d0_muon[imu])>0.1 :
#                    print "Electron and muon entering N3"
#                    print "d0 electron is " , iev.d0_electron[iele]
#                    print "d0 muon is " , iev.d0_muon[imu]
#                    N3[d]=N3[d]+1*weight
#                    print N3
#
    
            # define the non-overlaping singal region (SR) out of the displaced -regions
            SR1=N1-N2
            SR2=N2-N3
            SR3=N3
    
            # Fill the two D array for clearer output
            datasetArray[0]= d.attrib['name']
            datasetArray[1]= d.attrib['title']
            datasetArray[2]= float(d.attrib['EqLumi'])
            datasetArray[3]=N1
            datasetArray[4]=N2
            datasetArray[5]=N3
            datasetArray[6]=SR1
            datasetArray[7]=SR2
            datasetArray[8]=SR3

            doubleArray.append(datasetArray)
            print "double array is"
            print doubleArray

        
            # compute the sum of background in the 3 SR
            Sum_SR1=Sum_SR1+SR1
            Sum_SR2=Sum_SR2+SR2
            Sum_SR3=Sum_SR3+SR3

        
            # end of event loop
            idataset=idataset+1    
    
    
    # print the final number per SR    
    print "we expect " , Sum_SR1 ,"events in the SR1"               
    print "we expect " , Sum_SR2 ,"events in the SR2"               
    print "we expect " , Sum_SR3 ,"events in the SR3"               
    
    # print the summary contained it the double array
#    for i in range (0,len(doubleArray)):
    for i in range (0,idataset):
        print "---------"
        print "NEW SAMPLE!!!"
        print "---------"
        for j in range (0,len(doubleArray[i])):
            print doubleArray[i][j]
        print ""
        print ""
    
    # print the minimum necessary for the Yield
#    for i in range (0,len(doubleArray)):
    for i in range (0,idataset):
        print "---------"
        print "Sample is ", doubleArray[i][1], "and the yield in the signal regions are:"
        print "---------"
        for j in (6,7,8):
            print doubleArray[i][j]
        print ""
        print ""
    
    
    print "end of the program !!!!"
    

    # end of sample loop
    
    

