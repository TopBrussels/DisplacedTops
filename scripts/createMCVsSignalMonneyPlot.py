#############
"""
This script will use the information produced by the create2DHistoForFactorization.py script and make a single
plot with the Sum of the background vs a signal sample in the 2D d0 plane. 

"""
#############




#basic import
from array import array
import array
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
        leptonStr="electron"
        print "In ElEl final state!! \n"
    if "MuMu" in chan:
        histo_index=0
        leptonStr="muon"
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
        histo.SetXTitle(leptonStr+"1 d_{0} [cm]")
        histo.SetYTitle(leptonStr+"2 d_{0} [cm]")
        histo.GetYaxis().SetTitleOffset(1.4)

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


    #set the colors and size for the legend
    histLineColor = rt.kOrange+7
    histFillColor = rt.kOrange-2
    markerSize  = 1.0
    
    latex = rt.TLatex()
    n_ = 2
    
    x1_l = 0.92
    y1_l = 0.60

    dx_l = 0.30
    dy_l = 0.18
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l
    
    legend =  rt.TPad("legend_0","legend_0",x0_l,y0_l,x1_l, y1_l )
    #legend.SetFillColor( rt.kGray )
    legend.Draw()
    legend.cd()

    ar_l = dy_l/dx_l
    #gap_ = 0.09/ar_l
    gap_ = 1./(n_+1)
    bwx_ = 0.12
    bwy_ = gap_/1.5

    x_l = [1.2*bwx_]
    #y_l = [1-(1-0.10)/ar_l]
    y_l = [1-gap_]
    ex_l = [0]
    ey_l = [0.04/ar_l]
    
   #array must be converted 
    x_l = array.array("f",x_l)
    ex_l = array.array("f",ex_l)
    y_l = array.array("f",y_l)
    ey_l = array.array("f",ey_l)
    
    gr_l =  rt.TGraphErrors(1, x_l, y_l, ex_l, ey_l)
    
    rt.gStyle.SetEndErrorSize(0)
    gr_l.SetMarkerSize(0.9)
    gr_l.Draw("0P")
    
    latex.SetTextFont(42)
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    
    latex.SetTextSize(0.25)    
    latex.SetTextAlign(12) 
    
#    box_ = rt.TBox()
    xx_ = x_l[0]
    yy_ = y_l[0]
    latex.DrawLatex(xx_+1.*bwx_,yy_,"RPV stops, c#tau 1 cm")

    yy_ -= gap_
#    box_.SetLineStyle( rt.kSolid )
#    box_.SetLineWidth( 1 )
    # box_.SetLineColor( kBlack )
#    box_.SetLineColor( histLineColor )
#    box_.SetFillColor( histFillColor )
#    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
#    box_.SetFillStyle(0)
#    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
    #Draw Z->ee text
    latex.DrawLatex(xx_+1.*bwx_,yy_,"SM background ")

    #update the canvas to draw the legend
    canvas.Update()

    raw_input("Press Enter to end")

    canvas.Print("plots/MonneyPlot"+chan+".pdf")
# eo loop over the channels


    
    
        
