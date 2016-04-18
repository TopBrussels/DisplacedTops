##############                                                                         
# pyroot macro to calcultate the Yield for different all the NonQCD sample. Compares yield when just cutting on the d0 and when using the parametrisation to extend the prediction
# April 2016 by qpython@cern.ch 
#############  


import os, sys
from array import array
from ROOT import *
import ROOT as rt
from uncertainties import ufloat
import numpy as n


#channel 
channels=["_ElEl","_MuMu"]

# samples
#dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "ZToll", "NonQCD"]
dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "ZToll"]
inputFiles=[]

for sample in dataSetTitles:
    inputFiles.append(TFile("rootFiles/"+sample+"2D.root"))



# hadd all non QCD to make the sum
"""cmd"""


# loop over the channel (ee or mumu)
i_chan = 0
for chan in channels:
    
    if "ElEl" in chan:
        histo_index=1
        print "In ElEl final state!! \n"
    if "MuMu" in chan:
        histo_index=0
        print "In MuMu final state!! \n"
    

    # loop over samples
    i_sam=0
    for sample in dataSetTitles:
        
        # declare vector necessary for the 4 graps
        xValues=[]

        yCut=[]
        yFact=[]
        yFactUp=[]
        yFactDown=[]

        
        histnames = []
        rootkeys = inputFiles[i_sam].GetListOfKeys()
        for key in rootkeys:
            histnames.append(key.GetName())
    
        histo=inputFiles[i_sam].Get(histnames[histo_index])
        print histo
    
        
        # set some usefull varibles
    
        #Nbin
        NbinsX=histo.GetXaxis().GetNbins()
        NbinsY=histo.GetYaxis().GetNbins()
        if not (NbinsX == NbinsY):
            print "the number of bin in the x-axis (", NbinsX, ") and the number of the bin in the y-axis (", NbinsY, ") is different!!!"
        else:
            print "NbinsX=NbinsY=",NbinsX
        
    
        # NTot
        NTot_error=rt.Double()
        NTot=histo.IntegralAndError(0,histo.GetXaxis().GetNbins()+2,0,histo.GetYaxis().GetNbins()+2,NTot_error);
        NTot_ = ufloat (NTot, NTot_error)
    
        print "NTot_ is ", NTot_ , "\n"
        #
        
    
    
        # loop over bins
    #    ibin = 0
        for ibin in range (1,NbinsX+1):
            print "ibin is ", ibin
    
    
            # cut on X and Y and count
            Ncut_error=rt.Double()
            Ncut=histo.IntegralAndError(ibin,histo.GetXaxis().GetNbins()+2,ibin,histo.GetYaxis().GetNbins()+2,Ncut_error);
            Ncut_ = ufloat (Ncut, Ncut_error)
            print "Ncut_ is " , Ncut_
    
    
            ##### 

            # cut on X and count
            Ncutx_error=rt.Double()
            Ncutx = histo.IntegralAndError(ibin,histo.GetXaxis().GetNbins()+2,0,histo.GetYaxis().GetNbins()+2,Ncutx_error);
            Ncutx_ = ufloat (Ncutx, Ncutx_error)
            print "Ncutx is ", Ncutx
            
            #  efficiency of first lepton (X)
            eff_x = Ncutx/NTot
            print "eff_x is ", eff_x
    
    

            # cut on Y and count
            Ncuty_error=rt.Double()
            Ncuty = histo.IntegralAndError(0,histo.GetXaxis().GetNbins()+2,ibin,histo.GetYaxis().GetNbins()+2,Ncuty_error);
            Ncuty_ = ufloat(Ncuty, Ncuty_error)
            print "Ncuty is ", Ncuty
    
            #  efficiency of first lepton (Y)
            eff_y = Ncuty/NTot
            print "eff_y is ", eff_y
    
    
            # get the yield with the factorisation Nfact = NTot * effx * effy = Ncutx * Ncuty / NTot
            Nfact_ = Ncutx_ * Ncuty_ / NTot_
            print "Nfact_ is ", Nfact_
    
            #Filling the vector for the 4 Graphs
            xValues.append(ibin/100.)

            yCut.append(Ncut)
            yFact.append(Nfact_.nominal_value)
            yFactUp.append(Nfact_.nominal_value + Nfact_.std_dev)
            yFactDown.append(Nfact_.nominal_value - Nfact_.std_dev)

    
            ####
    
            
            ibin = ibin + 1
            # eo loop over bins


        # convert into array of double to allow compability wiht TGraph
        xValuesDouble = array("d",xValues)

        yCutDouble = array ("d",yCut)
        yFactDouble = array("d",yFact)
        yFactUpDouble = array("d",yFactUp)
        yFactDownDouble = array("d",yFactDown)            
        print xValuesDouble

        sampleColour=3

        # defining the Graphs
        gCut=rt.TGraph(len(xValuesDouble), xValuesDouble,yCutDouble)
        gCut.SetTitle(sample)
        gCut.SetLineColor(kGray)
        gCut.SetMarkerStyle(21)
        gCut.SetMarkerColor(kGray)
#        gCut.SetLineWidth(3.5)
        

        # Central value 
        gFact=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactDouble)
        gFact.SetTitle(sample)
        gFact.SetLineColor(sampleColour)

        # Up Value
        gFactUp=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactUpDouble)
        gFactUp.SetTitle(sample)
        gFactUp.SetLineColor(sampleColour)

        # Down Value
        gFactDown=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactDownDouble)
        gFactDown.SetTitle(sample)
        gFactDown.SetLineColor(sampleColour)
        
        

        # define a Canvas
        canv=rt.TCanvas("c1"+sample+chan)

        # draw all the graphs
        gCut.Draw("pa")
        gFact.Draw("l0")
        gFactUp.Draw("l0")
        gFactDown.Draw("l0")
        gCut.Draw("p")


        canv.SetLogy()
        canv.Print("param"+sample+chan+".gif")


        i_sam=i_sam+1
        # eo loop over 
    

    i_chan=i_chan+1
    #eo loop over channels

    
