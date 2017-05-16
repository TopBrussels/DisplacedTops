
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
import numpy as np
import json

import facoLib as fl


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
#channels=["_ElEl"]
channels=["_ElEl","_MuMu"]

# path to tree
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"
# 
folderName="Systematics_29_8_2016/"


# root file wit d0 distribution
inputFile = rt.TFile("rootFiles/d0forTFsNew.root",'READ')

h_electrond0 = inputFile.Get("electrond0")
h_muond0 = inputFile.Get("muond0")



# empty dict of hist
histo_dict = {}


# dictrionary of the average TFs
N_QCD_dict = {}
N_QCD_err_dict = {}
N_QCD_syst_dict = {}


#
wps = ["Loose", "Medium", "Tight"]
lepts = ["electron" , "muon"]
binTypes = ["", "Wide"]

# Fill the dictionary of histo
# loop over wps
for wp in wps  :
    # loop over leptons
    for lept in lepts:
        # loop over binnig type
        for binType in binTypes:
            hist_str = lept+"d0"+binType+"_"+wp
            histo_dict[hist_str] = inputFile.Get(hist_str)
            
                




# output root file                                                                       
outfile = rt.TFile("rootFiles/TFS.root",'RECREATE')


#
# convert number of the SR into its definition (SR1-> LSR, SR2 -> MSR)
dictSRX_toXSR = {1:"L", 2: "M" , 3 :"T"}


# Calculate all the TFs (DCR->LSR,MSR,TSR) 
# loop over leptons
for lept in lepts:
    
    SRx_Header = ["Region"]
    SRx_doubleArray = []

    # histo
    fancyHist = rt.TH1D("TF"+lept,"TF"+lept, 15, 0.5, 15.5)

    print "\n"
    print "lepton is ", lept

    # Normalisation region yiel
    yield_DCR = 0
    yield_DCR_err = 0
    
    if lept is "electron":
        yield_DCR = 0.1
        yield_DCR_err = 2.5
    if lept is "muon" :
        yield_DCR = 5.406
        yield_DCR_err = 2.61

    yield_DCR_ = ufloat (yield_DCR, yield_DCR_err)



    # loop over the 3 SRs
    for i_SR, SR in enumerate(range (1, 4)):
        

        SRx_SingleArray = ["SR"+str(i_SR+1)]
        

        TFs_valueAndErrorArray_ = []

        TFs_list = []
        TFs_err_list = []


        # loop over the wps
        for i_wp, wp in enumerate(wps) :
            
            
            
            print "\n"
            print "Working point is ", wp

            hist_str = lept+"d0Wide"+"_"+wp
        
            # yield in LIDCR with error
            N_DCR = histo_dict[hist_str].GetBinContent(1)
            N_DCR_err = histo_dict[hist_str].GetBinError(1)
            N_DCR_ = ufloat (N_DCR, N_DCR_err)

            # yield in SRX with error
            N_SR = histo_dict[hist_str].GetBinContent(SR+1)
            N_SR_err = histo_dict[hist_str].GetBinError(SR+1)
            N_SR_ = ufloat (N_SR, N_SR_err) 

            # get TF
            TF_ = (N_SR_/N_DCR_)**2
            TF = TF_.nominal_value
            TFs_list.append(TF)
            TF_err = TF_.std_dev
            TFs_err_list.append(TF_err)
            
            SRx_Header.append(wp)
            SRx_SingleArray.append(TF_)
            TFs_valueAndErrorArray_.append(TF_)
            
            
            print "TF (DCR->SR"+str(SR)+") is " ,  TF_

            # put everything in a histogram
            ibin = 2 + i_SR * 5 + i_wp

            # fill 
            fancyHist.SetBinContent(ibin, TF)
            fancyHist.SetBinError(ibin, TF_err)

            # set bin label
            fancyHist.GetXaxis().SetBinLabel(ibin, wp)
        
        # eo over the WP
    
    
        # calculate average TF for each SR
        TF_mean = np.mean(TFs_list)
        TF_max_err = max(TFs_err_list) # take the max error
        TF_mean_ = ufloat (TF_mean, TF_max_err)
        print "average with max error is ", TF_mean_
        
        SRx_SingleArray.append(TF_mean_)
        SRx_Header.append("mean")


        # get the QCD yield
        N_QCD_ = yield_DCR_ * TF_mean_
        print "yield DCR is ", yield_DCR, "and TF_mean_ is ", TF_mean_
        SRx_SingleArray.append(N_QCD_)
        SRx_Header.append("N_QCD")


#        print "List for the current SR (" , i_SR , ") is ", TFs_valueAndErrorArray_
        
        # properly calculated average if we assume that the TF from all wp are  uncorelated..
        # which is not the case
        #        TF_mean_ = np.mean(TFs_valueAndErrorArray_)
        #        print "average is " , TF_mean_
    
        # fill
        fancyHist.SetBinContent(5 + i_SR * 5,TF_mean_.nominal_value)
        fancyHist.SetBinError(5 + i_SR * 5, TF_mean_.std_dev)
    
        # set bin label
        fancyHist.GetXaxis().SetBinLabel(1 + i_SR * 5, dictSRX_toXSR[SR]+"SR" )
        fancyHist.GetXaxis().SetBinLabel(5 + i_SR * 5, "Average" )

        # save the value in the dictionary
        N_QCD_dict["SR" + str(SR) + "_" + lept] = N_QCD_.nominal_value
        N_QCD_err_dict["SR" + str(SR) + "_" + lept] = N_QCD_.std_dev



        
        # calculate the systematic uncertainty
        TF_sys = (max(TFs_list) - min(TFs_list))
        TF_sys_rel = TF_sys/ TF_mean
        SRx_SingleArray.append(fl.floatToPercent(TF_sys_rel))
        SRx_Header.append("Systematic")

        # put it in a dict
        N_QCD_syst_dict["SR" + str(SR) + "_" + lept] = TF_sys_rel


        # put the single array in the double array
        SRx_doubleArray.append(SRx_SingleArray)

    # eo loop over the SRs

    # style
    fancyHist.SetMaximum(1)
    fancyHist.SetMinimum(0.0001)
    fancyHist.SetYTitle("TF(DCR->[L,M,T]SR)")
    fancyHist.SetStats(False)

#    fancyHist.SetTitle("")
    

    # canvas
    outfile.cd()
    c1 = rt.TCanvas(lept)
    c1.cd()
    fancyHist.Draw()
    c1.SetLogy()
    fancyHist.Write()
    c1.Write()
    c1.SaveAs("plots/TFs_" + lept + ".pdf")


    # make a table
    caption = "Values of the transfer factors and the yields in the " + str(lept) + "s final state."
    fl.makeTable("QCDSumaryTable_" + lept, SRx_doubleArray, SRx_Header, True, caption , False, True)



#eo lopp over the lept


outfile.Close()


# save the dict
with open ("jsonFiles/N_QCD.json", 'w') as f:
    json.dump(N_QCD_dict, f)

with open ("jsonFiles/N_QCD_err.json", 'w') as f:
    json.dump(N_QCD_err_dict, f)

with open ("jsonFiles/N_QCD_syst.json", 'w') as f:
    json.dump(N_QCD_syst_dict, f)




# checking input histogram
"""
boundsLept1 = []

print "printing histo .."
for ibin in range (1,h_electrond0.GetNbinsX()+2):
    print "ibin is  ", ibin , " and corresponding lowEdge is " , h_electrond0.GetBinLowEdge(ibin)
    print "appending lowEdge value in a vector \n"
    boundsLept1.append(h_electrond0.GetBinLowEdge(ibin))

print "vector of bound is"
print boundsLept1
print "\n making a loop with this array"
for boundLept1 in boundsLept1 :
    ibin1 = h_electrond0.FindBin(boundLept1) -1 
    
    print "boundLept1 is ", boundLept1 , " annd ibin1  is ", ibin1 
    


boundsLept1P = [0.010,0.011,0.012, 0.013,0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.020]
print "using a hand defined array .."
print boundsLept1P
for boundLept1P in boundsLept1P :
    ibin1 = h_electrond0.FindBin(boundLept1P) -1

    print "boundLept1P is ", boundLept1P , " and ibin1  is ", ibin1
"""













# Manage the number of print outs
debug=False

# define various bound on the first and second leptonn
binsLept1 = [2, 5, 8]
binsLept2 = [2, 5, 8]


# define used histo for the Closure test histograms names
h_muond0 = histo_dict["muond0_Tight"]
h_electrond0 = histo_dict["electrond0_Tight"]

# random lumivalue
lumivalue = 3

# double array for table writting
doubleArray = []

i_lept1 = 0
# loop over the bound of the first lepton
for binLept1 in binsLept1 :

    i_lept2 = 0
    # loop over the bound of the first lepton
    for binLept2 in binsLept2 :
    
        print "binLept1 is" , binLept1
        print "binLept2 is" , binLept2
        
        # loop over the channel (lepton in final state)
        for chan in channels:
 
        
            # loop over the low bounds
        #    for ilb in range (1,len(LowBounds)+1):
        #        dict = ('bgMCSum+')
        #        bgMCSum
            
            # Histogram containing the number of events
            NonQCDBase=rt.TH1D("NonQCDBase"+chan+str(binLept1)+str(binLept2),"NonQCDBase",1,0,1)
            DataBase=rt.TH1D("DataBase"+chan+str(binLept1)+str(binLept2),"DataBase",1,0,1)
            QCDBase=rt.TH1D("QCDBase"+chan+str(binLept1)+str(binLept2),"QCDBase",1,0,1)
            NonQCDTarget=rt.TH1D("NonQCDTarget"+chan+str(binLept1)+str(binLept2),"NonQCDTarget",1,0,1)
            DataTarget=rt.TH1D("DataTarget"+chan+str(binLept1)+str(binLept2),"DataTarget",1,0,1)
            QCDTarget=rt.TH1D("QCDTarget"+chan+str(binLept1)+str(binLept2),"QCDTarget",1,0,1)
            EstimatedQCDTarget=rt.TH1D("EstimatedQCDTarget"+chan+str(binLept1)+str(binLept2),"EstimatedQCDTarget",1,0,1)
        
        
        
            isElEl=False
            isMuMu=False
            isElMu=False        

            # get the xmlfile that corresponds to the channel
            if "MuMu" in chan:
                isMuMu=True
                tree = ET.ElementTree(file='../config/MuMuV4.xml')
                treeName="tree"
                FinalState="At least two muons"
                print FinalState
            elif "ElEl" in chan:
                isElEl=True
                tree = ET.ElementTree(file='../config/ElElV4.xml')
                treeName="tree"
                FinalState="At least two electrons"
                print FinalState
            elif "ElMu" in chan:
                tree = ET.ElementTree(file='../config/ElMuV4.xml')
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
        
        
                    ch.Add(pathTrunc+folderName+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")
                    print pathTrunc+folderName+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root"
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
                            PileUpWeight=iev.evt_puSF
                        if isElEl:
                            PileUpWeight=iev.evt_puSF
                        
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
                            d01=abs(iev.d0BeamSpot_electron[0])
                            d02=abs(iev.d0BeamSpot_electron[1])

                        if (isMuMu):
                            d01=abs(iev.d0BeamSpot_muon[0])
                            d02=abs(iev.d0BeamSpot_muon[1])

                        if (isElMu):
                            d01=abs(iev.d0BeamSpot_electron[0])
                            d02=abs(iev.d0BeamSpot_muon[0])

                            
                        # skip events outside range of interest
                        if d01 < 0.01 or 0.02 < d01 :
                            continue
                        if d02 < 0.01 or 0.02 < d02 :
                            continue

            
                        # get the lepton scale factor depending on the channel
                        if isElEl:
                            for ilept in range (0,2):
                                LeptonWeight *= iev.sf_electron[ilept]         
               
                        if isMuMu:
                            for ilept in range (0,2):
                                LeptonWeight *= iev.sf_muon[ilept]

                        if isElMu:
                            LeptonWeight= iev.sf_muon[0]*iev.sf_electron[0]


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
            
                
                        # convert bin into bound
                        boundLept1= lepton1.GetBinLowEdge(binLept1+1)
                        boundLept2= lepton2.GetBinLowEdge(binLept2+1)

                            
                        

                        # event in base if both lepton are smaller than a bound but still in DCR
                        if 0.01 < d01 and  d01 < boundLept1 and 0.01 < d02 and d02 < boundLept2:
                            isInBase=True
                            if (debug):
        #                    if (True):
                                print "Leptons entering Base"
                                print "d0 muon[0] is " , iev.d0BeamSpot_muon[0]
                                print "d0 muon[1] is " , iev.d0BeamSpot_muon[1]
                                    
        
                        # event in base if both lepton are bigger than previous bound but still in DCR
                        if boundLept1 < d01 and d01 < 0.02 and boundLept2 < d02 and d02 < 0.02:
                            isInTarget=True
                            if (debug):
                                print "Lepton entering Target"
                                print "d0 muon is " , iev.d0BeamSpot_muon[ilept]
        
        
        
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
            

            

            # example to calculate error (looks like it is just the sqrt of the integral..)
#            Base = lepton1.Integral(binLept1+1,10)
#            Base_err=rt.Double()
#            lepton1.IntegralAndError(binLept1+1,10,Base_err)
#            print "Base is " , Base
#            print "Base_err is ", Base_err

            TF = 0

            TF1 = lepton1.Integral(binLept1+1,10)/ lepton1.Integral(1,binLept1)
            TF1_err = TF1 * (1/math.sqrt(lepton1.Integral(binLept1+1,10)) + 1/math.sqrt(lepton1.Integral(1,binLept1)) )
            TF1_ = ufloat (TF1,TF1_err)
#            print "TF1_err/TF1 is ", TF1_err/TF1


            TF2 = lepton2.Integral(binLept2+1,10)/ lepton2.Integral(1,binLept2)
            TF2_err =  TF2 * (1/math.sqrt(lepton2.Integral(binLept2+1,10)) + 1/math.sqrt(lepton2.Integral(1,binLept2)) )
            TF2_ = ufloat (TF2, TF2_err)
            
            TF = TF1 * TF2
            TF_ = TF1_ * TF2_
            
            print "TF is ", TF
            print "TF_ is ", TF_

            print "boundLept1 , binLept1 , TF1 , boundLept2, binLept2 , TF2 , TF are ..."
            print boundLept1 , binLept1 , TF1 , boundLept2, binLept2 , TF2 , TF

            
            NQCDBase = DataBase.GetBinContent(1) - NonQCDBase.GetBinContent(1)
            NQCDBase_err = DataBase.GetBinError(1)+ NonQCDBase.GetBinError(1)
            NQCDBase_ = ufloat (DataBase.GetBinContent(1),DataBase.GetBinError(1)) - ufloat(NonQCDBase.GetBinContent(1),NonQCDBase.GetBinError(1))

            EstimatedQCDTarget.Fill(0.5,TF * NQCDBase )
            EstimatedQCDTarget_ = NQCDBase_ * TF1_ * TF2_
            DirectQCDTarget_ = ufloat (DataTarget.GetBinContent(1),DataTarget.GetBinError(1)) - ufloat(NonQCDTarget.GetBinContent(1),NonQCDTarget.GetBinError(1))

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
            singleArray = [str(binLept1)+" ; "+str(binLept2),DataTarget.GetBinContent(1),  DataTarget.GetBinError(1), EstimatedQCDTarget.GetBinContent(1), EstimatedQCDTarget.GetBinError(1) ] 


            singleArray_ = [str(binLept1)+" ; "+str(binLept2), DirectQCDTarget_, EstimatedQCDTarget_ ] 


            
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
#headers=["bounds","DirectCount","Error","EstimatedCount","Error"]
headers=["bounds","Direct Count", "Estimated Count"]
print tabulate(doubleArray, headers, tablefmt="latex")

# writing results in a tex file                                                                   
outputFile = "tables/ClosureTestTable"+chan+".tex"
fout = open (outputFile, "w")
fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
fout.write("\\begin{table}"+newLine)
fout.write("\\caption{ " + "QCD Closure Test"+chan.replace("_"," ")+ "}"+newLine)

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

