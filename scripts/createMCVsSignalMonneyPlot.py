#############
"""
This script will use the information produced by the create2DHistoForFactorization.py script and make a single
plot with the Sum of the background vs a signal sample in the 2D d0 plane. 

"""
#############




#basic import
from array import array
from ROOT import *
import ROOT as rt

# list of channels
channels=["_ElEl","_MuMu"]


# list of root file names
dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "DrellYann"]
dataSetColours=[38, 5, 46, 872, 30, 45]



inputFiles=[]

for sample in dataSetTitles:
    inputFiles.append(TFile("rootFiles/"+sample+"2D.root"))



# loop over the channel (ee or mumu) 
i_chan = 0
for chan in channels:
    if "ElEl" in chan:
        histo_index=1
        print "In ElEl final state!! \n"
    if "MuMu" in chan:
        histo_index=0
        print "In MuMu final state!! \n"

    # declare one 2D histo per channel (final state)
    electrond0VsElectronsd0Sum=rt.TH2D("electrond0VsElectronsd0Sum","electrond0VsElectronsd0", 50, 0.0, 0.05, 50, 0.0, 0.05)
    muond0VsMuond0Sum=rt.TH2D("muond0VsMuond0Sum","muond0VsMuond0",50, 0.0, 0.05, 50, 0.0, 0.05)
    muond0VsElectrond0Sum=rt.TH2D("muond0VsElectrond0Sum","muond0VsElectrond0",50, 0.0, 0.05, 50, 0.0, 0.05)
    

    canvas=rt.TCanvas("c2")
    canvas.cd()


    # loop over the samples
    i_sam=0
    for sample in dataSetTitles:


        histnames = []
        rootkeys = inputFiles[i_sam].GetListOfKeys()
        for key in rootkeys:
            histnames.append(key.GetName())
    
        histo=inputFiles[i_sam].Get(histnames[histo_index])
        print histo

        histo.SetMarkerColor(dataSetColours[i_sam])
        if (i_sam == 0):
            histo.Draw("")
        else :
            histo.Draw("SAME")
        
        i_sam=i_sam+1
    # eo loop over the samples


    canvas.Print("plot.pdf")


    
    
        
