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

#channel
channels=["_MuMu"]
#channels=["_ElEl"]
#channels=["_ElEl","_MuMu"]

# path to tree
folderName="NoBlindingRerun_30_11_2016"
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

# debug
debug=False

treeName="tree"


lumivalue = 3


# loop over the channel (lepton in final statue)
for chan in channels:
    
    bgMCSum1=rt.TH1D("bgMCSum1"+chan,"bgMCSum1",1,0,1)



    isElEl=False
    isMuMu=False
    doubleArray=[]

    # get the xmlfile that corresponds to the channel
    if "MuMu" in chan:
        isMuMu=True
        myFile='../config/MuMuV4.xml'
        FinalState="At least two muons"
        print FinalState
    elif "ElEl" in chan:
        isElEl=True
        myFile='../config/ElElV4.xml'
        FinalState="At least two electrons"
        print FinalState
    elif "ElMu" in chan:
        myFile='../config/ElMuV0.xml'
    else:
        print "No tree will be loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

    tree = ET.ElementTree(file=myFile)

        

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []

    # loop over datasets
    for idataset, d in enumerate(datasets):
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
            
            if "Data" not in sampleName and "stop" not in sampleName:
                isBgMC = True
            if "stop" in sampleName:
                isSignal = True
            if "Data" in sampleName:
                isData = True


            ch.Add(pathTrunc+folderName+"/"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+chan+"DCR.root")

            N1=rt.TH1D("N1"+sampleName+chan,"N1",1,0,1)




            # get number of events
            nevents=ch.GetEntries()
            
    
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

                if isMuMu:
                    PileUpWeight=iev.evt_puSF
                if isElEl:
                    PileUpWeight=iev.evt_puSF
                
                LeptonWeight=1.0
                

    
                # loop over the 2 highest pt letpon
                for ilept in range (0,2):

                    # make the logic for the muon
                    if isMuMu:
                        LeptonWeight *= iev.sf_muon[ilept]                        


                    # make the logic for the electron
                    if isElEl :
                        LeptonWeight *= iev.sf_electron[ilept]
                        #for bound in range(0,len(SRxBounds)):                                                                                                                                          

                if isData :
                    PileUpWeight=1
                    LeptonWeight=1
                            

                N1.Fill(0.5,weight*PileUpWeight*LeptonWeight)


                if isBgMC :
                    bgMCSum1.Fill(0.5,weight*PileUpWeight*LeptonWeight)
                    
            
            # eo event loop


    
            # Fill the two D array for clearer output
            datasetArray = [ d.attrib['name'], d.attrib['title'], d.attrib['EqLumi'], nevents, N1.GetBinContent(1),N1.GetBinError(1)]
            print datasetArray
            
            
#
            doubleArray.append(datasetArray)
            print "double array is"
            print doubleArray

        
            # end of event loop

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
    line = " &  d0 $>$ "+str("bound1")
#    line = " &  d0 $>$ "+str("bound1")+ " & d0 $>$ " + str("bound2") + " & d0 $>$ " + str("bound3")
    line = line + endLine + hLine
    fout.write(line)
    
    # write a line for each background sample
    bgMCcounter = 0
    for i in range (0,len(doubleArray)):
        if "Data" not in doubleArray[i][0] and "NP" not in doubleArray[i][0]:
            bgMCcounter = bgMCcounter + 1
            rawlabel = doubleArray[i][0]
            label = rawlabel.replace("","").replace("#tau","$\\tau$").replace("\Nu","$\\nu$").replace("\rightarrow","${\\rightarrow}$").replace(" ","\\ ").replace("_","")
            line = label + " & " + str(round(doubleArray[i][4],3))
            #line = label + " & " + str(round(doubleArray[i][iSR1],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR1+1],3)) + " & " + str(round(doubleArray[i][iSR2],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR2+1],3)) + " &  " + str(round(doubleArray[i][iSR3],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR3+1],3))
            line = line + endLine + newLine + hLine
            fout.write(line)

    #write a line with the sum of the backgrounds
    if bgMCcounter is not 0:

        line = hLine+"Total expected background & " + str(round(bgMCSum1.GetBinContent(1),3))
        #line = hLine+"Total expected background & " + str(round(bgMCSum1.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum1.GetBinError(1),3)) + " & " + str(round(bgMCSum2.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum2.GetBinError(1),3)) +" & " + str(round(bgMCSum3.GetBinContent(1),3)) + " $\pm$ " + str(round(bgMCSum3.GetBinError(1),3))
        line = line + endLine + newLine + hLine
        fout.write(line)

    #write a line for each data sample
    for i in range (0,len(doubleArray)):
        if "Data" not in doubleArray[i][0]:
            continue
        label = "Observation"
        myYield= doubleArray[i][4]
#        myYield= doubleArray[i][iSR1]
        fout.write(label + " & " + str(myYield) + endLine + newLine) 

    # write a line for each signal sample
    for i in range (0,len(doubleArray)):
        if "NP" not in doubleArray[i][0]:
            continue
        rawlabel = doubleArray[i][0]
        label = rawlabel.replace("_","")
#        myYield= doubleArray[i][iSR1]
        line = label + " & " + str(round(doubleArray[i][4],3))
#        line = label + " & " + str(round(doubleArray[i][iSR1],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR1+1],3)) + " & " + str(round(doubleArray[i][iSR2],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR2+1],3)) + " &  " + str(round(doubleArray[i][iSR3],3)) + " $\pm$ "+ str(round(doubleArray[i][iSR3+1],3))
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
