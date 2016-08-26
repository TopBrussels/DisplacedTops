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
channels=["_MuMu"]
#channels=["_ElEl"]
#channels=["_ElEl","_MuMu"]

# path to tree
folderName="CMSSW76V4_NewCutFlow"
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
        tree = ET.ElementTree(file='../config/Yield_FullSamplesMuMuV0.xml')
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        tree = ET.ElementTree(file='../config/Yield_FullSamplesElElV0.xml')
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


            ch.Add(pathTrunc+folderName+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+".root")


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

                    sf_iso_muon = sf_iso_muon_up[ilept] 
                    print "sf_iso_muon is " ,sf_iso_muon
                    sf_id_muon = sf_id_muon[ilept]
                    print "sf_id_muon is ", sf_id_muon
                    

                    # make the logic for the muon
                    if isMuMu:
                        LeptonWeight *= sf_iso_muon * sf_id_muon

                        if abs(iev.d0BeamSpot_muon[ilept]) < bound1:
                            passed1=False
                            if (debug):
                                print "Electron and muon entering N1"
                                print "d0 muon is " , iev.d0BeamSpot_muon[ilept]
                            

                    # eo the logic for the muon 


                                
                    # make the logic for the electron
                    if isElEl :
                        LeptonWeight *= 1
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
                print " PileUpWeight is ", PileUpWeight
                print "LeptonWeight is ", LeptonWeight
                
                if (passed1==True):
                    N1.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
            
            # eo event loop


    
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
    outputFile = "tables/YieldTable"+chan+".tex"
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
