##############
# Pyroot macro to calcultate the Yield for different Signal region and create a compilable tex file.
# February 2016 by qpython@cern.ch
#############

from tabulate import tabulate
from uncertainties import ufloat
import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt
import math


# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"


# testing uncertainties propagation
#x = ufloat(1, 0.25)
#print x**2
#print x-x
#print (x**2+ x +1000).derivatives[x]


# using Lorentz Vectors (lv) as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()

#channel
#channels=["_MuMu"]
channels=["_ElEl"]
#channels=["_ElEl","_MuMu"]

# path to antiIso tree
date="1_4_2016"
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"


# root file wit d0 distribution
inputFile = rt.TFile("d0forTFs.root",'READ')

h_electrond0=inputFile.Get("electrond0")
h_muond0=inputFile.Get("muond0")


# debug
debug=False

# define various bound on the first and second leptonn
boundsLept1 = [0.012, 0.015, 0.018]
boundsLept2 = [0.012, 0.015, 0.018]

# define vector of Transfer Factors
TFs = [[0,0,0],[0,0,0],[0,0,0]]


# random lumivalue
lumivalue = 3

# double array for table writting
doubleArray = []

i_lept1 = 0
# loop over the bound of the first lepton
for boundLept1 in boundsLept1 :

    i_lept2 = 0
    # loop over the bound of the first lepton
    for boundLept2 in boundsLept2 :
    
        print "boundLept1 is" , boundLept1
        print "boundLept2 is" , boundLept2
        
        # loop over the channel (lepton in final state)
        for chan in channels:
        
            # loop over the low bounds
        #    for ilb in range (1,len(LowBounds)+1):
        #        dict = ('bgMCSum+')
        #        bgMCSum
            
            # Histogram containing the number of events
            NonQCDBase=rt.TH1D("NonQCDBase"+chan+str(boundLept1)+str(boundLept2),"NonQCDBase",1,0,1)
            DataBase=rt.TH1D("DataBase"+chan+str(boundLept1)+str(boundLept2),"DataBase",1,0,1)
            QCDBase=rt.TH1D("QCDBase"+chan+str(boundLept1)+str(boundLept2),"QCDBase",1,0,1)
            NonQCDTarget=rt.TH1D("NonQCDTarget"+chan+str(boundLept1)+str(boundLept2),"NonQCDTarget",1,0,1)
            DataTarget=rt.TH1D("DataTarget"+chan+str(boundLept1)+str(boundLept2),"DataTarget",1,0,1)
            EstimatedQCDTarget=rt.TH1D("EstimatedQCDTarget"+chan+str(boundLept1)+str(boundLept2),"EstimatedQCDTarget",1,0,1)
        
            # 
            NBase=rt.TH1D("NBase"+chan+str(boundLept1)+str(boundLept2),"NBase",1,0,1)
            NTarget=rt.TH1D("NTarget"+chan+str(boundLept1)+str(boundLept2),"NTarget",1,0,1)
        
        
            isElEl=False
            isMuMu=False
            isElMu=False        

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
                treeName="tree"
                FinalState="One muon and one electron"
                print FinalState
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
#                    print str(d.attrib['name'])
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
                    if (0):
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


                        # define shorter variable depending on the channel
                        if (isElEl):
                            d01=abs(iev.d0BeamSpot_electron_elel[0])
                            d02=abs(iev.d0BeamSpot_electron_elel[1])

                        if (isMuMu):
                            d01=abs(iev.d0BeamSpot_muon_mumu[0])
                            d02=abs(iev.d0BeamSpot_muon_mumu[1])

                        if (isElMu):
                            d01=abs(iev.d0BeamSpot_electron_elel[0])
                            d02=abs(iev.d0BeamSpot_muon_mumu[0])

                            
                        # skip events outside range of interest
                        if d01 < 0.01 or 0.02 < d01 :
                            continue
                        if d02 < 0.01 or 0.02 < d02 :
                            continue

            
                        # get the lepton scale factor depending on the channel
                        if isElEl:
                            for ilept in range (0,2):
                                LeptonWeight *= iev.sf_electron_elel[ilept]         
               
                        if isMuMu:
                            for ilept in range (0,2):
                                LeptonWeight *= iev.sf_muon_mumu[ilept]

                        if isElMu:
                            LeptonWeight= iev.sf_muon_mumu[0]*iev.sf_electron_elel[0]

                                

                            
                        

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
        
        
        
                        # Filling the histogram
                        if (isData) :
                            PileUpWeight=1
                            LeptonWeight=1
                                    
                        # filling base histo
                        if (isInBase):
        #                    print "isInbase!!"
                            if (isData):
                                DataBase.Fill(0.5,weight*PileUpWeight*LeptonWeight)
#                                QCDBase.Fill(0.5,-1)
                            if (isBgMC and not isQCD):
                                NonQCDBase.Fill(0.5,weight*PileUpWeight*LeptonWeight)
#                                QCDBase.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                                
                                
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


            # Get the TFs
            

            
            # get the correct histograms depending on the channel
            if (isElEl):
                lepton1=h_electrond0
                lepton2=h_electrond0

            if (isMuMu):
                lepton1=h_muond0
                lepton2=h_muond0

            if (isElMu):
                lepton1=h_electrond0
                lepton2=h_muond0
            
                
            # convert cut into bin Number
            ibin1 = lepton1.FindBin(boundLept1) - 1
            ibin2 = lepton2.FindBin(boundLept2) - 1


            # example to calculate error (looks like it is just the sqrt of the integral..)
#            Base = lepton1.Integral(ibin1+1,10)
#            Base_err=rt.Double()
#            lepton1.IntegralAndError(ibin1+1,10,Base_err)
#            print "Base is " , Base
#            print "Base_err is ", Base_err

            TF = 0

            TF1 = lepton1.Integral(ibin1+1,10)/ lepton1.Integral(1,ibin1)
            TF1_err = TF1 * (1/math.sqrt(lepton1.Integral(ibin1+1,10)) + 1/math.sqrt(lepton1.Integral(1,ibin1)) )
            TF1_ = ufloat (TF1,TF1_err)
#            print "TF1_err/TF1 is ", TF1_err/TF1


            TF2 = lepton2.Integral(ibin2+1,10)/ lepton2.Integral(1,ibin2)
            TF2_err =  TF2 * (1/math.sqrt(lepton2.Integral(ibin2+1,10)) + 1/math.sqrt(lepton2.Integral(1,ibin2)) )
            TF2_ = ufloat (TF2, TF2_err)
            
            TF = TF1 * TF2
            TF_ = TF1_ * TF2_
            
            print "TF is ", TF
            print "TF_ is ", TF_

            print "boundLept1 , ibin1 , TF1 , boundLept2, ibin2 , TF2 , TF are ..."
            print boundLept1 , ibin1 , TF1 , boundLept2, ibin2 , TF2 , TF

            
            NQCDBase = DataBase.GetBinContent(1)-NonQCDBase.GetBinContent(1)
            NQCDBase_ = ufloat (DataBase.GetBinContent(1),DataBase.GetBinError(1)) - ufloat(NonQCDBase.GetBinContent(1),NonQCDBase.GetBinError(1))
#            print "NEstimateQCDBase is " , NEstimateQCDBase
#            print "NEstimateQCDBase_ is " , NEstimateQCDBase_
            EstimatedQCDTarget.Fill(0.5,TF * NQCDBase )
            EstimatedQCDTarget_ = NQCDBase_ * TF1_ * TF2_
            DirectQCDTarget_ = ufloat (DataTarget.GetBinContent(1),DataTarget.GetBinError(1)) - ufloat(NonQCDTarget.GetBinContent(1),NonQCDTarget.GetBinError(1))

            NQCDBase_err = DataBase.GetBinError(1)+NonQCDBase.GetBinError(1)
            CombinedError = EstimatedQCDTarget.GetBinContent(1) * (TF1_err/TF1 + TF2_err/TF2 + NQCDBase_err/NQCDBase )
            EstimatedQCDTarget.SetBinError(1,CombinedError )

#            EstimatedQCDTarget.Fill(0.5,TF * QCDBase.GetBinContent(1) )
        
        
            # make some printout
            print  "NonQCDBase is ", NonQCDBase.GetBinContent(1) 
            print  "DataBase is ", DataBase.GetBinContent(1)  
            print  "NonQCDTarget is ", NonQCDTarget.GetBinContent(1)  
            print  "DataTarget is ", DataTarget.GetBinContent(1) , " +/- " , DataTarget.GetBinError(1)
            print  "EstimatedQCDTarget is ", EstimatedQCDTarget.GetBinContent(1) , " +/- " , EstimatedQCDTarget.GetBinError(1)
                                


            
                    # Fill the two D array for clearer output
            singleArray = [str(boundLept1)+" ; "+str(boundLept2),DataTarget.GetBinContent(1),  DataTarget.GetBinError(1), EstimatedQCDTarget.GetBinContent(1), EstimatedQCDTarget.GetBinError(1) ] 


            singleArray_ = [str(boundLept1)+" ; "+str(boundLept2), DirectQCDTarget_, EstimatedQCDTarget_ ] 


            
            doubleArray.append(singleArray_)
#            doubleArray[i_lept1].append(singleArray)
            print "doubleArray is " , doubleArray
    
        i_lept2=i_lept2+1
    i_lept1=i_lept1+1



            
# print the summary contained it the double array
for i in range (0,len(doubleArray)):
    #    for i in range (0,lendidataset):
    print "---------"
    print "NEW bound pair!!!"
    print "---------"
    for j in range (0,len(doubleArray[i])):
        print doubleArray[i][j]


# get the info for the table
headers=["bounds","DirectCount","Error","EstimatedCount","Error"]
print tabulate(doubleArray, headers, tablefmt="latex")

# writing results in a tex file                                                                   
outputFile = "ClosureTestTable"+chan+".tex"
fout = open (outputFile, "w")
fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
fout.write("\\begin{table}"+newLine)
fout.write("\\caption{ " + "QCD Closure Test"+chan+ "}"+newLine)

# the actual tabular
fout.write(tabulate(doubleArray, headers, tablefmt="latex"))


# end of table                                                                   
fout.write("\\end{table}"+newLine)
fout.write("\\end{document}"+newLine)
fout.close()


# write a line for each value of the bound on the the second letpon

        
"""
        
                    datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], nevents, NonQCDBase.GetBinContent(1),NonQCDBase.GetBinError(1), DataBase.GetBinContent(1), DataBase.GetBinError(1), NonQCDTarget.GetBinContent(1),NonQCDTarget.GetBinError(1), DataTarget.GetBinContent(1),DataTarget.GetBinError(1), EstimatedQCDTarget.GetBinContent(1),EstimatedQCDTarget.GetBinError(1)]
                    print datasetArray
                    
                    
            
            
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

