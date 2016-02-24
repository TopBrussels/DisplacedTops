##############
# pyroot macro to calcultate the Yield for different Signal region and create a compilable tex file.
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
#channels=["_MuMu"]
#channels=["_ElEl"]
channels=["_ElEl","_MuMu"]

# path to tree
date="12_2_2016"
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

# debug
debug=False

# define the bound of the Signal region
#SRxBounds=[0.001,0.002,0.1]
bound1 = 0.001
bound2 = 0.02
bound3 = 0.05


lumivalue = 3


# loop over the channel (lepton in final statue)
for chan in channels:
    
    bgMCSum1=rt.TH1D("bgMCSum1"+chan,"bgMCSum1",1,0,1)    
    bgMCSum2=bgMCSum1.Clone("bgMCSum2"+chan)
    bgMCSum3=bgMCSum1.Clone("bgMCSum3"+chan)
    isElEl=False
    isMuMu=False
    doubleArray=[]

    # get the xmlfile that corresponds to the channel
    if "MuMu" in chan:
        isMuMu=True
        tree = ET.ElementTree(file='../config/TreeProc_FullSamplesMuMuV0.xml')
        treeName="doubleMuTree"
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        tree = ET.ElementTree(file='../config/TreeProc_FullSamplesElElV0.xml')
        treeName="doubleElTree"
        FinalState="At least two electrons"
        print FinalState
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElMuV0.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

        

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    idataset=0

    # index to get the bin conten of the 3 different signal region
    iSR1=10
    iSR2=12
    iSR3=14

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
            isSignal = False
            isData = False
            
            if "Data" not in sampleName and "NP" not in sampleName:
                isBgMC = True
            if "NP" in sampleName:
                isSignal = True
            if "Data" in sampleName:
                isData = True


            ch.Add(pathTrunc+date+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")
            Sum_SR1=Sum_SR2=Sum_SR3=0
            
            # define 3 histograms containing inclusive yields
            N1=rt.TH1D("N1"+sampleName+chan,"N1",1,0,1)
            N2 = N1.Clone("N2"+sampleName+chan)
            N3 = N1.Clone("N3"+sampleName+chan)

            # define 3 histograms containing exclusive yields
            SR1 = N1.Clone("SR1"+sampleName+chan)
            SR2 = N1.Clone("SR2"+sampleName+chan)
            SR3 = N1.Clone("SR3"+sampleName+chan)



            # you can also use two bins, one for the weighted event and one for the unweighted event
#            N1=rt.TH1D("N1"+sampleName,"N1",2,0,2)

            
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
#                print  d.attrib['title']," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
                ii+=1
                passed1= True
                passed2= True
                passed3= True
    
                # loop over the 2 highest pt letpon
                for ilept in range (0,2):

                    # make the logic for the muon
                    if isMuMu:
                        LeptonWeight *= iev.sf_muon_mumu[ilept]                        

                    # if one of the leptons  is smaller than bound, the event fails
                        if abs(iev.d0BeamSpot_muon_mumu[ilept]) < bound1:
                            passed1=False
                            if (debug):
                                print "Electron and muon entering N1"
                                print "d0 muon is " , iev.d0BeamSpot_muon_mumu[ilept]
                            

                        if abs(iev.d0BeamSpot_muon_mumu[ilept]) < bound2:
                            passed2=False
                            if (debug):
                                print "Electron and muon entering N2"
                                print "d0 muon is " , iev.d0BeamSpot_muon_mumu[ilept]

                        if abs(iev.d0BeamSpot_muon_mumu[ilept]) < bound3:
                            passed3=False
                            if (debug):
                                print "Electron and muon entering N3"
                                print "d0 muon is " , iev.d0BeamSpot_muon_mumu[ilept]
                    # eo the logic for the muon 


                                
                    # make the logic for the electron
                    if isElEl :
                        LeptonWeight *= iev.sf_electron_elel[ilept]
                        #for bound in range(0,len(SRxBounds)):                                                                                                                                          


                        # if one of the leptons  is smaller than bound, the event fails
                        if abs(iev.d0BeamSpot_electron_elel[ilept]) < bound1:
                            passed1=False
                            if (debug):
                                print "Electron and muon entering N1"
                                print "d0 electron is " , iev.d0BeamSpot_electron_elel[ilept]


                        if abs(iev.d0BeamSpot_electron_elel[ilept]) < bound2:
                            passed2=False
                            if (debug):
                                print "Electron and muon entering N2"
                                print "d0 electron is " , iev.d0BeamSpot_electron_elel[ilept]

                        if abs(iev.d0BeamSpot_electron_elel[ilept]) < bound3:
                            passed3=False
                            if (debug):
                                print "Electron and muon entering N3"
                                print "d0 electron is " , iev.d0BeamSpot_electron_elel[ilept]
                    # eo the logic for the electron



                # Filling the histogram
                if isData :
                    PileUpWeight=1
                    LeptonWeight=1
                            

                if (passed1==True):
                    N1.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
                if (passed2==True):
                    N2.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
                if (passed3==True):
                    N3.Fill(0.5,weight*PileUpWeight*LeptonWeight)


                if (passed1==True) and (passed2==False) and (passed3==False):
                    SR1.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    if isBgMC:
                        bgMCSum1.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
                elif (passed1==True) and (passed2==True) and (passed3==False):
                    SR2.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    if isBgMC:
                        bgMCSum2.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
                elif (passed1==True) and (passed2==True) and (passed3==True):
                    SR3.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    if isBgMC:
                        bgMCSum3.Fill(0.5,weight*PileUpWeight*LeptonWeight)
            
            # eo event loop


    
            # Fill the two D array for clearer output
            datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], nevents, N1.GetBinContent(1),N1.GetBinError(1), N2.GetBinContent(1),N2.GetBinError(1), N3.GetBinContent(1),N3.GetBinError(1), SR1.GetBinContent(1),SR1.GetBinError(1), SR2.GetBinContent(1),SR2.GetBinError(1), SR3.GetBinContent(1),SR3.GetBinError(1)]
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
        for j in (7,8,9):
            print doubleArray[i][j]
        print ""
        print ""

    # print the final number per SR    
    for i in range (0,idataset):
        if "Data" not in doubleArray[i][0]  and "NP" not in doubleArray[i][0]:
            Sum_SR1=Sum_SR1+doubleArray[i][iSR1]
            Sum_SR2=Sum_SR2+doubleArray[i][iSR2]
            Sum_SR3=Sum_SR3+doubleArray[i][iSR3]
            

    print "we expect " , Sum_SR1 ,"events in the SR1"               
    print "we expect " , Sum_SR2 ,"events in the SR2"               
    print "we expect " , Sum_SR3 ,"events in the SR3"               

    

    # writing results in a tex file
    outputFile = "YieldTable"+chan+".tex"
    fout = open (outputFile, "w")
    fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
    fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
    fout.write("\\begin{table}"+newLine)
    fout.write("\\caption{ " + FinalState + "}"+newLine)
    fout.write("\\begin{tabular}{lrrr}"+newLine+hLine)

    line = "Event Source &  \multicolumn{3}{c}{Event Yield $\pm$ 1$\sigma$ (stat.)} "
    line = line +endLine+newLine+hLine+hLine
    fout.write(line)
    line = " &  d0 $>$ "+str(bound1)+ " & d0 $>$ " + str(bound2) + " & d0 $>$ " + str(bound3)
    line = line + endLine + hLine
    fout.write(line)
    
    # write a line for each background sample
    bgMCcounter = 0
    for i in range (0,len(doubleArray)):
        if "Data" not in doubleArray[i][0] and "NP" not in doubleArray[i][0]:
            bgMCcounter = bgMCcounter + 1
            rawlabel = doubleArray[i][0]
            label = rawlabel.replace("","").replace("#tau","$\\tau$").replace("\Nu","$\\nu$").replace("\rightarrow","${\\rightarrow}$").replace(" ","\\ ").replace("_","")
            line = label + " & " + str(round(doubleArray[i][iSR1],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR1+1],3)) + " & " + str(round(doubleArray[i][iSR2],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR2+1],3)) + " &  " + str(round(doubleArray[i][iSR3],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR3+1],3))
            line = line + endLine + newLine + hLine
            fout.write(line)

    #write a line with the sum of the backgrounds
    if bgMCcounter is not 0:
        line = hLine+"Total expected background & " + str(round(bgMCSum1.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum1.GetBinError(1),3)) + " & " + str(round(bgMCSum2.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum2.GetBinError(1),3)) +" & " + str(round(bgMCSum3.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum3.GetBinError(1),3))
        line = line + endLine + newLine + hLine
        fout.write(line)

    #write a line for each data sample
    for i in range (0,len(doubleArray)):
        if "Data" not in doubleArray[i][0]:
            continue
        label = "Observation"
        myYield= doubleArray[i][iSR1]
        fout.write(label + " & " + str(myYield) + endLine + newLine) 

    # write a line for each signal sample
    for i in range (0,len(doubleArray)):
        if "NP" not in doubleArray[i][0]:
            continue
        rawlabel = doubleArray[i][0]
        label = rawlabel.replace("_","")
#        myYield= doubleArray[i][iSR1]
        line = label + " & " + str(round(doubleArray[i][iSR1],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR1+1],3)) + " & " + str(round(doubleArray[i][iSR2],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR2+1],3)) + " &  " + str(round(doubleArray[i][iSR3],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR3+1],3))
        line = line + endLine + newLine + hLine
        fout.write(line)


        
    # end of tabular
    fout.write("\\end{tabular}"+newLine)
    fout.write("\\end{table}"+newLine)
    fout.write("\\end{document}"+newLine)
    fout.close()



    # end of sample loop


print "end of the program !!!!"
# end of the channel loop
