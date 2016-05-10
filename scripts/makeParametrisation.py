##############                                                                         
# pyroot macro to calcultate the Yield for different all the NonQCD sample. Compares yield when just cutting on the d0 and when using the parametrisation to extend the prediction
# Also produces single letpon efficiencies on the same canva
# April 2016 by qpython@cern.ch 
#############  

from tabulate import tabulate
import os, sys
from array import array
from ROOT import *
import ROOT as rt
from uncertainties import ufloat
import numpy as n




# usefull variables for writing a tex file 
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"


#channel 
channels=["_ElEl","_MuMu"]

# samples
dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "ZToll", "NonQCD"]
dataSetTitles=["WJets", "Diboson", "SingleTop", "TTJets", "ZToll"]
dataSetColours=[38, 5, 46, 872, 30, 45]


inputFiles=[]

for sample in dataSetTitles:
    inputFiles.append(TFile("rootFiles/"+sample+"2D.root"))

# bool to select upper Edge of the integrals
includeOverFlowBin=True



# loop over the channel (ee or mumu)
i_chan = 0
for chan in channels:
    
    if "ElEl" in chan:
        histo_index=1
        print "In ElEl final state!! \n"
    if "MuMu" in chan:
        histo_index=0
        print "In MuMu final state!! \n"

    # array for the yield 
    yieldArray=[]
    Sum = 0
    Sum_error = 0
    Sum_ = ufloat (Sum, Sum_error)
    

    # loop over samples
    i_sam=0
    for sample in dataSetTitles:
        
        # declare vector necessary for the 4 graps
        xValues=[]
        xValues_error=[]

        yCut=[]
        yCut_error=[]


        yEff_y=[]
        yEff_x=[]
        yFact=[]
        yFactUp=[]
        yFactDown=[]

        
        histnames = []
        rootkeys = inputFiles[i_sam].GetListOfKeys()
        for key in rootkeys:
            histnames.append(key.GetName())
    
        histo=inputFiles[i_sam].Get(histnames[histo_index])
        print histo


        # check correlation
        histo.GetCorrelationFactor()
        print "sample is ", sample , " ,channel is ", chan , " ,CorrelationFactor is ", histo.GetCorrelationFactor()

    
        
        # set some usefull varibles
        #Nbin
        NbinsX=histo.GetXaxis().GetNbins()
        NbinsY=histo.GetYaxis().GetNbins()
        if not (NbinsX == NbinsY):
            print "the number of bin in the x-axis (", NbinsX, ") and the number of the bin in the y-axis (", NbinsY, ") is different!!!"
        else:
            print "NbinsX=NbinsY=",NbinsX

        # Upper edge for the integral
        MaxX=NbinsX
        MaxY=NbinsY
        if (includeOverFlowBin):
            MaxX=NbinsX+1
            MaxY=NbinsY+1
        
    
        # NTot
        NTot_error=rt.Double()
        NTot=histo.IntegralAndError(0, MaxX, 0, MaxY, NTot_error);
        NTot_ = ufloat (NTot, NTot_error)
    
        print "NTot_ is ", NTot_ , "\n"
        #
        
    
    
        # loop over bins
    #    ibin = 0
        for ibin in range (0,NbinsX+1):
#            print "ibin is ", ibin
    
    
            # cut on X and Y and count
            Ncut_error=rt.Double()
            Ncut=histo.IntegralAndError(ibin, MaxX, ibin, MaxY, Ncut_error);
            Ncut_ = ufloat (Ncut, Ncut_error)
#            print "Ncut_ is " , Ncut_
    
    
            ##### 

            # cut on X and count
            Ncutx_error=rt.Double()
            Ncutx = histo.IntegralAndError(ibin, MaxX, 0, MaxY, Ncutx_error);
            Ncutx_ = ufloat (Ncutx, Ncutx_error)
#            print "Ncutx is ", Ncutx
            
            #  efficiency of first lepton (X)
            eff_x = Ncutx/NTot
#            print "eff_x is ", eff_x
    
    

            # cut on Y and count
            Ncuty_error=rt.Double()
            Ncuty = histo.IntegralAndError(0, MaxX, ibin, MaxY, Ncuty_error);
            Ncuty_ = ufloat(Ncuty, Ncuty_error)
#            print "Ncuty is ", Ncuty
    
            #  efficiency of first lepton (Y)
            eff_y = Ncuty/NTot
#            print "eff_y is ", eff_y
    
    
            # get the yield with the factorisation Nfact = NTot * effx * effy = Ncutx * Ncuty / NTot
            Nfact_ = Ncutx_ * Ncuty_ / NTot_
#            print "Nfact_ is ", Nfact_
    
            #Filling the vector for the 4 Graphs
            xValues.append(ibin/1000.)
            xValues_error.append(0.)

            yCut.append(Ncut)
            yCut_error.append(Ncut_error)

            yFact.append(Nfact_.nominal_value)
            yFactUp.append(Nfact_.nominal_value + Nfact_.std_dev)
            yFactDown.append(Nfact_.nominal_value - Nfact_.std_dev)

            # vector for efficiencies
            yEff_x.append(eff_x)
            yEff_y.append(eff_y)

            # get the yield for the SR
            if ibin == 20 and not dataSetTitles[i_sam] == "NonQCD":
                singleArray=[sample,Nfact_]
                yieldArray.append(singleArray)
                Sum_=Sum_+Nfact_
            

            ####
    
            
            ibin = ibin + 1
            # eo loop over bins


        # convert into array of double to allow compability wiht TGraph
        xValuesDouble = array("d",xValues)
        xValues_errorDouble = array("d",xValues_error)

        yCutDouble = array ("d",yCut)
        yCut_errorDouble = array ("d",yCut_error)

        yFactDouble = array("d",yFact)
        yFactUpDouble = array("d",yFactUp)
        yFactDownDouble = array("d",yFactDown)            
#        print xValuesDouble

        yEff_xDouble = array ("d",yEff_x)
        yEff_yDouble = array ("d",yEff_y)
        
        # defining the Graphs

        # cut and count
        gCut=rt.TGraphErrors(len(xValuesDouble), xValuesDouble,yCutDouble, xValues_errorDouble,yCut_errorDouble)
        gCut.SetTitle(sample+chan)
        gCut.SetLineColor(kGray)
        gCut.SetMarkerStyle(21)
        gCut.SetMarkerColor(kGray)
        gCut.SetLineWidth(1)
        

        # Factorised Central value 
        gFact=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactDouble)
        gFact.SetTitle(sample+chan)
        gFact.SetLineColor(dataSetColours[i_sam])
        gFact.SetLineWidth(3)

        # Factorised Up Value
        gFactUp=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactUpDouble)
        gFactUp.SetTitle(sample+chan)
        gFactUp.SetLineColor(dataSetColours[i_sam])
        gFactUp.SetLineStyle(2)
        gFactUp.SetLineWidth(3)

        # Factorised Down Value
        gFactDown=rt.TGraph(len(xValuesDouble), xValuesDouble,yFactDownDouble)
        gFactDown.SetTitle(sample+chan)
        gFactDown.SetLineColor(dataSetColours[i_sam])
        gFactDown.SetLineStyle(2)
        gFactDown.SetLineWidth(3)
        
        

        # define a Canvas
        canv=rt.TCanvas("c1"+sample+chan)

        # draw all the graphs
        gCut.Draw("pa")
        gFact.Draw("l0")
        gFactUp.Draw("l0")
        gFactDown.Draw("l0")
        gCut.Draw("p")

        # make the legend box
        leg = rt.TLegend(0.5,0.7,0.9,0.85)
        leg.SetFillColor(kWhite)
        leg.SetBorderSize(0)
        
        # creates a line
        line = rt.TLine(0.1,0.2,0.3,0.4)
        line.SetLineStyle(2)

        # add the entries
        leg.AddEntry(gCut,sample+"MC from cut-and-count method","pl")
        leg.AddEntry(gFact,sample+"MC from factorized method (#pm 1 #sigma)","l")
        leg.AddEntry(line,"edges of signal/prompt regions","l")
        leg.Draw("same")
        canv.SetLogy()

        # fiddle with y range for appropriate display
        gCut.GetHistogram().SetMaximum(2*gCut.GetHistogram().GetMaximum())
        gCut.GetHistogram().SetMinimum(0.000005);
        
        # axis labels 
        gCut.GetHistogram().GetXaxis().SetTitle("d_{0} > x cut value [cm]")
        gCut.GetHistogram().GetYaxis().SetTitle("predicted events after d_{0} cut")

        # add lines that define transition between regions
        lowval =gCut.GetHistogram().GetMinimum()
        highval =gCut.GetHistogram().GetMaximum()
        line.DrawLine(0.01,lowval,0.01,highval)
        line.DrawLine(0.02,lowval,0.02,highval)
        leg.Draw("same")


        # save the canva
        canv.Modified()
#        canv.Print("plots/param"+sample+chan+".gif")
        canv.Print("plots/param"+sample+chan+".pdf")



        ##############################################
        # bo new graph with the single lepton efficiencies
        #

        # eff of first lepton
        gEff_x=rt.TGraph(len(xValuesDouble), xValuesDouble,yEff_xDouble)
        gEff_x.SetTitle(sample+chan)
        gEff_x.SetLineColor(dataSetColours[i_sam])
        gEff_x.SetLineStyle(2)
        gEff_x.SetLineWidth(3)


        # eff of second lepton
        gEff_y=rt.TGraph(len(xValuesDouble), xValuesDouble,yEff_yDouble)
        gEff_y.SetTitle(sample+chan)
        gEff_y.SetLineColor(dataSetColours[i_sam])
        gEff_y.SetLineWidth(3)
        
         # define a Canvas
        c2=rt.TCanvas("c2"+sample+chan)

        # draw the graphs
        gEff_x.Draw("l0a")
        gEff_y.Draw("l0")


        # make the legend box
        leg2 = rt.TLegend(0.5,0.7,0.9,0.85)
        leg2.SetFillColor(kWhite)
        leg2.SetBorderSize(0)
        
        # add entries to the legend
        leg2.AddEntry(gEff_x,"efficiency of the first lepton","l")
        leg2.AddEntry(gEff_y,"efficiency of the second lepton","l")
        leg2.Draw()

        gEff_x.GetHistogram().SetMinimum(0.0001);
        
        # set axis title
        gEff_x.GetHistogram().GetXaxis().SetTitle("d_{0} > x cut value [cm]")
        gEff_x.GetHistogram().GetYaxis().SetTitle("efficiency")
        
        # save canva
        c2.SetLogy()
        c2.Modified()
        c2.Print("plots/eff"+sample+chan+".pdf")

        #
        # eo new graph with the single lepton efficiencies
        ##############################################


        i_sam=i_sam+1
        # eo loop over sample

    # put the sum of all nonQCD sample in the array
    singleArray=["Sum of NonQCD background",Sum_]
    yieldArray.append(singleArray)



    # writing results in a tex file 
    fileName="parametrisationTable"+chan
    outputFile = "tables/"+fileName+".tex"
    fout = open (outputFile, "w")
    fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
    fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
    fout.write("\\begin{table}"+newLine)
    fout.write("\\caption{ " + "Yield estimated with the factorisation method in the SR (both lepton with d0 $<$ 0.02 cm) in the "+chan.replace("_"," ")+ " channel." "}"+newLine)

    # the actual tabular
    headers=["background source","Yield +/- uncertainty"]
    fout.write(tabulate(yieldArray, headers, tablefmt="latex"))

    # end of table                                                                                           
    fout.write("\\end{table}"+newLine)
    fout.write("\\end{document}"+newLine)
    fout.close()

    # compile tex into pdf
    cmd="pdflatex tables/"+fileName+".tex"
    os.system(cmd)
    
    # mv table
    cmd="mv "+fileName+".pdf"+" tables/"
    os.system(cmd)
    
    # clean the mess
    cmd="rm "+fileName+".aux"
    os.system(cmd)
    cmd="rm "+fileName+".log"
    os.system(cmd)
    

    i_chan=i_chan+1
    #eo loop over channels

