import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt
import tdrstyle, CMS_lumi
##############
# example pyroot loop for yield counting on output trees of Ntupler
# March 2015 by qpython@cern.ch
#


# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"



# using Lorentz Vectors (lv) as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()


#channels=["_ElEl","_MuMu"]
channels=["_MuMu"]

date="12_2_2016"

# double array containing an array of samples. Each sample is an array with different varaibles such as ["fancyname", "name" x-sec,samplepresels,...] 
#doubleArray=[['Data', 'Data', '2629.405', 893018L, 346646, 33, 0, 346613, 33, 0], ['NP_overlay_stopTobl_m500_Ctau10', 'stopTobl_m500_Ctau10', '170941.984262', 9312L, 9103, 7789, 3322, 1314, 4467, 3322], ['WJetsToLNu', 'W\\rightarrow l+\\Nu', '392.612052979', 13L, 8, 0, 0, 8, 0, 0], ['WWToLNuQQ', 'Diboson', '139928.395704', 23L, 14, 0, 0, 14, 0, 0], ['WWTo2l2Nu', 'Diboson', '162587.288553', 42970L, 15089, 0, 0, 15089, 0, 0], ['ZZ', 'Diboson', '59650.1845912', 13757L, 4718, 0, 0, 4718, 0, 0], ['SingleTop_tW', 'SingletTop', '28089.8876404', 3573L, 1283, 0, 0, 1283, 0, 0], ['SingleTop_tbarW', 'SingletTop', '28073.0337079', 3585L, 1278, 0, 0, 1278, 0, 0], ['TTJets_Madgrap', 'TTJets_Madgraph', '12281.3443782', 33704L, 11864, 2, 0, 11862, 2, 0], ['QCD_MuEnriched_50to80', 'QCDMuEnriched', '11.5493477661', 1L, 0, 0, 0, 0, 0, 0], ['QCD_MuEnriched_80to120', 'QCDMuEnriched', '36.6018728008', 6L, 3, 0, 0, 3, 0, 0], ['QCD_MuEnriched_120to170', 'QCDMuEnriched', '137.512206973', 20L, 8, 2, 0, 6, 2, 0], ['QCD_MuEnriched_170to300', 'QCDMuEnriched', '455.559903971', 34L, 23, 6, 0, 17, 6, 0], ['QCD_MuEnriched_300to470', 'QCDMuEnriched', '4904.063158', 52L, 23, 3, 1, 20, 2, 1], ['QCD_MuEnriched_470to600', 'QCDMuEnriched', '24402.503292', 35L, 14, 0, 0, 14, 0, 0], ['QCD_MuEnriched_600to800', 'QCDMuEnriched', '79034.0011142', 38L, 13, 0, 0, 13, 0, 0], ['QCD_MuEnriched_800to1000', 'QCDMuEnriched', '421108.780958', 40L, 17, 1, 0, 16, 1, 0], ['QCD_MuEnriched_1000toInf', 'QCDMuEnriched', '1079300.96335', 49L, 23, 0, 0, 23, 0, 0], ['DYJetsToLL_M-50toInf_Madgraph', 'Z/\\gamma^{*}\\rightarrow ll', '1494.44466574', 539047L, 193576, 0, 0, 193576, 0, 0]]

doubleArray=[]

lumivalue = 3
debug=False

pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"
FinalState="At least two muons"

# define the bound of the Signal region
#SRxBounds=[0.001,0.002,0.1]
bound1 = 0.001
bound2 = 0.01
bound3 = 0.1

#bgMCSum=rt.TH1D("bgMCSum","bgMCSum",1,0,1)
bgMCSum=0


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
    iSR1=10
    iSR2=12
    iSR3=14

    # loop over datasets
    for d in datasets:
        if d.attrib['add'] == '1'  and "QCD_MuEnriched_120to170" in str(d.attrib['name']):
            print "found dataset to be added..." + str(d.attrib['name'])
            datasetNames.append(str(d.attrib['name']))
            print str(d.attrib['name'])
            # one array per dataset [name, title, Eqlumi, N1, N2, N3, SR1, SR2, SR3]
            ch = rt.TChain("doubleMuTree","doubleMuTree")
            sampleName=d.attrib['name']
            ch.Add("/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"+date+"/_MuMu/DisplacedTop_Run2_TopTree_Study_"+sampleName+"_MuMu.root")
            Sum_SR1=Sum_SR2=Sum_SR3=0
            N1=rt.TH1D("N1"+sampleName,"N1",1,0,1)
            # ideal way to do this (fix):
            # N2 = N1.Clone("N2"+samplename)
            N2=rt.TH1D("N2"+sampleName,"N2",1,0,1)
            N3=rt.TH1D("N3"+sampleName,"N3",1,0,1)
            SR1=rt.TH1D("SR1"+sampleName,"SR1",1,0,1)
            SR2=rt.TH1D("SR2"+sampleName,"SR2",1,0,1)
            SR3=rt.TH1D("SR3"+sampleName,"SR3",1,0,1)

            # you can also use two bins, one for the weighted event and one for the unweighted event
#            N1=rt.TH1D("N1"+sampleName,"N1",2,0,2)
            # N2 = N1.Clone("N2"+samplename)
            # ...
            
            # get number of events
            nevents=ch.GetEntries()
            
            if nevents == 0 :
                continue
            
    
            # calculate weight
            if "Data" in sampleName:
                lumivalue=float(d.attrib['EqLumi'])
                
            weight= lumivalue / float(d.attrib['EqLumi'])
            if (1):
                print "lumivalue is " ,lumivalue
                print " float(d.attrib['EqLumi']) is ",  float(d.attrib['EqLumi'])
                print "weight is " , weight
                
            ii=0
            # start of loop over events
            for iev in ch:

                PileUpWeight=iev.evt_puSF_mumu
                MuonWeight=1.0

                
#                if ii % (nevents/50.) ==0 :
#                print  d.attrib['title']," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
                ii+=1
                passed= True
                passed2= True
                passed3= True
    
        # loop over muons - fill in lorentz vector and fill some histograms
#                for imu in range(0,iev.nMuons_mumu) :
                # only the two leading muons
                SFWeight = 1.0

                for imu in range (0,2):
                    
                    
                    MuonWeight *= iev.sf_muon_mumu[imu]
    
#                    lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
#                    lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
                    

                        
#                    for bound in range(0,len(SRxBounds)):

                    
                    # if one leptons smaller than bound, event fails
                    if abs(iev.d0BeamSpot_muon_mumu[imu]) < bound1:
                        passed=False
                        if (debug):
                            print "Electron and muon entering N1"
                            print "d0 muon is " , iev.d0BeamSpot_muon_mumu[imu]
                            

                    if abs(iev.d0BeamSpot_muon_mumu[imu]) < bound2:
                        passed2=False
                        if (debug):
                            print "Electron and muon entering N2"
                            print "d0 muon is " , iev.d0BeamSpot_muon_mumu[imu]

                    if abs(iev.d0BeamSpot_muon_mumu[imu]) < bound3:
                        passed3=False
                        if (debug):
                            print "Electron and muon entering N2"
                            print "d0 muon is " , iev.d0BeamSpot_muon_mumu[imu]

                if "Data" in sampleName:
                    PileUpWeight=1
                    MuonWeight=1
                            

                if (passed==True):
                    N1.Fill(0.5,weight*PileUpWeight*MuonWeight)

                if (passed2==True):
                    N2.Fill(0.5,weight*PileUpWeight*MuonWeight)

                if (passed3==True):
                    N3.Fill(0.5,weight*PileUpWeight*MuonWeight)


                if (passed==True) and (passed2==False) and (passed3==False):
                    SR1.Fill(0.5,weight*PileUpWeight*MuonWeight)
                    
                elif (passed==True) and (passed2==True) and (passed3==False):
                     SR2.Fill(0.5,weight*PileUpWeight*MuonWeight)

                elif (passed==True) and (passed2==True) and (passed3==True):
                     SR3.Fill(0.5,weight*PileUpWeight*MuonWeight)







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
#            SR1=N1-N2
#            SR2=N2-N3
#            SR3=N3
    
            # Fill the two D array for clearer output
            datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], nevents, N1.GetBinContent(1),N1.GetBinError(1), N2.GetBinContent(1),N2.GetBinError(1), N3.GetBinContent(1),N3.GetBinError(1), SR1.GetBinContent(1),SR1.GetBinError(1), SR2.GetBinContent(1),SR2.GetBinError(1), SR3.GetBinContent(1),SR3.GetBinError(1)]
            print datasetArray
            
            
            
#            datasetArray[0]= d.attrib['name']
#            datasetArray[1]= d.attrib['title']
#            datasetArray[2]= float(d.attrib['EqLumi'])
#            datasetArray[3]=nevents
#            datasetArray[4]=N1
#            datasetArray[5]=N2
#            datasetArray[6]=N3
#            datasetArray[7]=SR1
#            datasetArray[8]=SR2
#            datasetArray[9]=SR3
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
        Sum_SR1=Sum_SR1+doubleArray[i][iSR1]
        Sum_SR2=Sum_SR2+doubleArray[i][iSR2]
        Sum_SR3=Sum_SR3+doubleArray[i][iSR3]
            

    print "we expect " , Sum_SR1 ,"events in the SR1"               
    print "we expect " , Sum_SR2 ,"events in the SR2"               
    print "we expect " , Sum_SR3 ,"events in the SR3"               

#    bgMCSum.Add()
    

    # writing results in a tex file
    outputFile = "YieldTable_MuMu.tex"
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
            myYield = doubleArray[i][iSR1]
            bgMCSum = bgMCSum+myYield
            line = label + " & " + str(round(myYield,3)) + " $\pm$ "+ str(round(doubleArray[i][iSR1+1],3)) + " & " + str(round(doubleArray[i][iSR2],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR2+1],3)) + " &  " + str(round(doubleArray[i][iSR3],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR3+1],3))
            line = line + endLine + newLine + hLine
            fout.write(line)

    #write a line with the sum of the backgrounds
    if bgMCcounter is not 0:
        line = hLine+"Total expected background & " + str(round(bgMCSum,3)) + " $\pm$ " + ""
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
        myYield= doubleArray[i][iSR1]
        fout.write(label + " & " + str(myYield) + endLine + newLine) 


        
    # end of tabular
    fout.write("\\end{tabular}"+newLine)
    fout.write("\\end{table}"+newLine)
    fout.write("\\end{document}"+newLine)
    fout.close()

        

    print "end of the program !!!!"


    
    

    # end of sample loop
    
    

