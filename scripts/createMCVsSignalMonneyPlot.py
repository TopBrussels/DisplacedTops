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
import CMS_lumi, tdrstyle


#set the tdr style
tdrstyle.setTDRStyle()


#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 10
if( iPos==0 ): CMS_lumi.relPosX = 0.12

H_ref = 600; 
W_ref = 800; 
W = W_ref
H  = H_ref


iPeriod = 4

# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.04*W_ref


# canvas
canvas = rt.TCanvas("c2","c2",50,50,W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)


# list of channels
channels=["_ElEl","_MuMu"]


# list of root file names
dataSetTitles=["Signal","WJets", "Diboson", "SingleTop", "TTJets", "DrellYann"]
dataSetColours=[1, 38, 5, 46, 872, 30, 45]



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

        # set the colour
        histo.SetMarkerColor(dataSetColours[i_sam])


        # draw signal and then all bkgd 
        if (i_sam == 0):
            histo.SetMarkerStyle(1)
            histo.Draw("")
        else :
            histo.SetMarkerStyle(4)
            histo.Draw("SAME")
        

            
        i_sam=i_sam+1
    # eo loop over the samples

    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    
    # update
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    frame = canvas.GetFrame()
    frame.Draw()


    canvas.Print("plots/MonneyPlot"+chan+".pdf")
# eo loop over the channels


    
    
        
