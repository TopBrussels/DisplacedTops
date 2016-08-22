############
#Pyroot script that calculates trigger SF.
# A data sample is used using MET trigger and is compare with ttbar MC.
# August 2016 by qpython@cern.ch 

# basic import
import xml.etree.cElementTree as ET
import os, sys
from array import array
import ROOT as rt


# list of channels
channels=["ttElEl","ttMuMu"]
#channels=["ttMuMu"] 
#channels=["ttElEl"] 

# list of samples
sampleNames=["DataRunD","TTJets_Dilept"]

# list of lepton indices (first pt or second pt)
indices=["0","1"]





#base of the path to the root file
pathTrunc="/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/"

#directory
directory="CMSSW76V4_TTLetp_11_8_2016"

# verbosity
debug = False

# fast run
fastRun = False

# name of the tree in the root file
treeName="tree"

# root output
outfile = rt.TFile("rootFiles/TriggerSF.root",'RECREATE')



# Make the plot comparison
c_ptComp=rt.TCanvas("Pt")
c_ptComp.Divide(2,2)

c_etaComp=rt.TCanvas("Eta")
c_etaComp.Divide(2,2)


canvaIndex=1
# loop over the different channels
for chan in channels :
    print "chan is ", chan

    isttElEl =False
    isttMuMu =False
    
    if chan == "ttElEl" :
        isttElEl = True
        lepton="Electron"
    if chan == "ttMuMu" :
        isttMuMu = True
        lepton="Muon"
            
    if isttElEl:
        tree = ET.ElementTree(file='../config/ttLeptonsV4.xml')    
    if isttMuMu:
        tree = ET.ElementTree(file='../config/ttLeptonsV4.xml')

    root =  tree.getroot()
    datasets = root.find('datasets')
    if (debug):
        print "found  "  + str(len(datasets)) + " datasets"



    # loop over the indices
    for index in indices :
        print "index is ", index        

        # list of all the efficiencies
        ptEffs=[]
        etaEffs=[]

        # TTBar or Data
        for sampleName in sampleNames:
            print "sampleName is ", sampleName

            
            ch = rt.TChain(treeName,treeName)
            ch.Add(pathTrunc+directory+"/_"+chan+"/DisplacedTop_Run2_TopTree_Study_"+sampleName+"_"+chan+".root")

            # string of the histo depending on the lepton (channel), the index and the samples name
            histStr=lepton+index
    
            # tot histo
            totPtLeptonIndexSample=rt.TH1D("totPt"+histStr+sampleName,"totPt"+histStr+sampleName, 8, 0.0, 200)
            totEtaLeptonIndexSample=rt.TH1D("totEta"+histStr+sampleName,"totEta"+histStr+sampleName, 6, 0.0, 2.4)
            
            # pass histo clone from tot to ensure bin consistency
            passPtLeptonIndexSample=totPtLeptonIndexSample.Clone("passPt"+histStr+sampleName)
            passEtaLeptonIndexSample=totEtaLeptonIndexSample.Clone("passEta"+histStr+sampleName)
    
        
            #TH2D 
            #passElectronPtVsElectronPt=rt.TH2D("passElectronPtVsElectronPt","passElectronPtVsElectronPt", 8, 0.0, 200, 8, 0.0, 200)
            #passElectronEtaVsElectronEta=rt.TH2D("passElectronEtaVsElectronEta","passElectronEtaVsElectronEta", 6, 0.0, 2.4, 6, 0.0, 2.4)
    
    
            #Clone 
            #MuonPtVsMuonPt=passElectronPtVsElectronPt.Clone("MuonPtVsMuonPt")
            # MuonEtaVsMuonEta=passElectronEtaVsElectronEta.Clone("MuonEtaVsMuonEta")
                    
        
            # start of loop over events 
            ii=0
            for iev in ch:
                ii=ii+1
                if fastRun and 100 < ii :
                    continue
        
                
            # Filling corresponding histograms
                if isttElEl :
    
                    # define all the variables
                    if index=="0":
                        ptIndexlept=iev.pt_electron[0]
                        etaIndexlept=iev.eta_electron[0]
                    
                    if index=="1":
                        ptIndexlept=iev.pt_electron[1]
                        etaIndexlept=iev.eta_electron[1]
    
                if isttMuMu :
                    if index=="0":
                        ptIndexlept=iev.pt_muon[0]
                        etaIndexlept=iev.eta_muon[0]
    
                    if index=="1":
                        ptIndexlept=iev.pt_muon[1]
                        etaIndexlept=iev.eta_muon[1]
        
        
                # fill tot histo
                totPtLeptonIndexSample.Fill(ptIndexlept)
                totEtaLeptonIndexSample.Fill(etaIndexlept)

                
                # fill pass histo
                if (iev.crossTrigged):
                    passPtLeptonIndexSample.Fill(ptIndexlept)
                    passEtaLeptonIndexSample.Fill(etaIndexlept)
    
            # eo loop over the event
    



            # make the efficiency with the properly filled histograms

                    
            # pt eff

            effPt=rt.TEfficiency(passPtLeptonIndexSample,totPtLeptonIndexSample)
            ptEffs.append(effPt)
            

            """
            # draw Eff on the canva
#            if sampleName != "TTJets_Dilept":
            canvPt=rt.TCanvas("Pt"+histStr)
            canvPt.cd()
            rt.gPad.SetGridy()
            effPt.Draw("AP")
            rt.gPad.Update()
            Pg_effPt = effPt.GetPaintedGraph()
            effPt.SetTitle("trigger efficiency;"+lepton+" "+index+" p_{T} [GeV]; #epsilon")
            Pg_effPt.SetMaximum(1.05)
            Pg_effPt.SetMinimum(0)
            

            # make legend 
            leg = rt.TLegend(0.7,0.4,0.9,0.5)
            leg.SetBorderSize(1)

#            else :
#                effPt.Draw("SAME")


            # add entries
            leg.AddEntry(effPt, sampleName , "l")

 #           if sampleName == "TTJets_Dilept":
            leg.Draw()            
            rt.gPad.Update()
            canvPt.Print("plots/EffPt"+histStr+sampleName+".pdf")
            print "saving histo"
            """




            # do it the dirty way
#            ratioPtLeptonIndexSample=passPtLeptonIndexSample.Clone("ratioPt"+histStr)
#            ratioPtLeptonIndexSample.Divide(totPtLeptonIndexSample)
#            ratioPtLeptonIndexSample.Write()


    
            # eta eff
            effEta=rt.TEfficiency(passEtaLeptonIndexSample,totEtaLeptonIndexSample)
            etaEffs.append(effEta)

            
            """
            # canva
            canvEta=rt.TCanvas("Eta"+histStr)
            canvEta.cd()
            rt.gPad.SetGridy()
            effEta.Draw("AP")
            rt.gPad.Update()
            Pg_effEta = effEta.GetPaintedGraph()
            effEta.SetTitle("trigger efficiency;"+lepton+" "+index+" #eta; #epsilon")
            Pg_effEta.SetMaximum(1.05)
            Pg_effEta.SetMinimum(0)


            # make legend
            leg = rt.TLegend(0.1,0.1,0.3,0.3)
            leg.SetBorderSize(1)


            # add entries
            leg.AddEntry(effEta, sampleName , "l")

            # draw the legend
            leg.Draw()

            # save plot
            rt.gPad.Update()
            canvEta.Print("plots/EffEta"+histStr+sampleName+".pdf")
            print "saving histo"
            """


            # write the histo and the canvas on the output root file
            outfile.cd()
            totPtLeptonIndexSample.Write()
            passPtLeptonIndexSample.Write()
#            canvPt.Write()

    
            totEtaLeptonIndexSample.Write()
            passEtaLeptonIndexSample.Write()
#            canvEta.Write()


        # do the SF using the efficiency of Data and TTJets

        if sampleName == "Data":
            gr_num=effPt.CreateGraph()
            gr_num.Write()



        

        #eo loop over the sampleName


        # make the plots comparison for each index, lepton pair
        c_pt=rt.TCanvas(histStr)
        c_pt.cd()
        
#        c_ptComp.cd(canvaIndex)
        ptEffs[0].Draw("AP")
        ptEffs[1].SetLineColor(8)
        ptEffs[1].Draw("Same")
        rt.gPad.SetGridy()

        rt.gPad.Update()
        Pg_effPt = ptEffs[0].GetPaintedGraph()
        ptEffs[0].SetTitle("trigger efficiency;"+lepton+" "+index+" p_{T} [GeV]; #epsilon")
        Pg_effPt.SetMaximum(1.05)
        Pg_effPt.SetMinimum(0)
            

        # make legend 
        leg = rt.TLegend(0.7,0.4,0.9,0.5)
        leg.SetBorderSize(1)


        # add entries
        leg.AddEntry(ptEffs[0], "Data" , "l")
        leg.AddEntry(ptEffs[1], sampleName , "l")

        
        leg.Draw()            
        rt.gPad.Update()
        c_ptComp.Update()


        c_pt.Update
        c_pt.Print("plots/EffPtComp"+histStr+".pdf")


        c_ptComp.Draw()
        print "YES???"
        print "canvaIndex is ", canvaIndex
        
        canvaIndex+=1


        # eta plots
        c_eta=rt.TCanvas(histStr)
        c_eta.cd()
        
#        c_etaComp.cd(canvaIndex)
        etaEffs[0].Draw("APE")
        etaEffs[1].SetLineColor(8)
        etaEffs[1].Draw("Same")
        rt.gPad.SetGridy()

        rt.gPad.Update()
        Pg_effEta = etaEffs[0].GetPaintedGraph()
        etaEffs[0].SetTitle("trigger efficiency;"+lepton+" "+index+" #eta; #epsilon")
        Pg_effEta.SetMaximum(1.05)
        Pg_effEta.SetMinimum(0)
            

        # make legend 
        leg = rt.TLegend(0.1,0.3,0.1,0.3)
        leg.SetBorderSize(1)


        # add entries
        leg.AddEntry(etaEffs[0], "Data" , "l")
        leg.AddEntry(etaEffs[1], sampleName , "l")

    
        leg.Draw()            
        rt.gPad.Update()
        c_etaComp.Update()


        c_eta.Update
        c_eta.Print("plots/EffEtaComp"+histStr+".pdf")


        c_etaComp.Draw()
        



        ### make the SF for each index and lepton pair
        
        # data pt eff hist from TEff
        Data_pt = totPtLeptonIndexSample.Clone("Data_pt"+histStr)
        for i_bin in range (0,ptEffs[0].GetPassedHistogram().GetNbinsX()):
            Data_pt.SetBinContent(i_bin,ptEffs[0].GetEfficiency(i_bin))
            print i_bin, " " , ptEffs[0].GetEfficiency(i_bin)
#        Data_pt.Sumw2()

        # MC pt eff hist from TEff   
        MC_pt = totPtLeptonIndexSample.Clone("MC_pt"+histStr)
        for i_bin in range (0,ptEffs[1].GetPassedHistogram().GetNbinsX()):
            MC_pt.SetBinContent(i_bin,ptEffs[1].GetEfficiency(i_bin))
            print i_bin, " " , ptEffs[1].GetEfficiency(i_bin)
 #       MC_pt.Sumw2()

        # SF pt hist
        SF_pt = Data_pt.Clone("SF_pt"+histStr)
        # Data/MC
        SF_pt.Divide(MC_pt)
        

        # Draw TH1D
        c_sf=rt.TCanvas("c_sf"+histStr)
        c_sf.cd()
        SF_pt.Draw("PE")
        SF_pt.SetMaximum(3)
        SF_pt.SetMinimum(0)
        c_sf.Print("plots/SFPt"+histStr+".pdf")
        



        # data eta eff hist from TEff                                                                                                                                             
        Data_eta = totEtaLeptonIndexSample.Clone("Data_eta"+histStr)
        for i_bin in range (0,etaEffs[0].GetPassedHistogram().GetNbinsX()):
            Data_eta.SetBinContent(i_bin,etaEffs[0].GetEfficiency(i_bin))
            print i_bin, " " , etaEffs[0].GetEfficiency(i_bin)
#        Data_eta.Sumw2()                                                                                                                                                         

        # MC eta eff hist from TEff                                                                                                                                               
        MC_eta = totEtaLeptonIndexSample.Clone("MC_eta"+histStr)
        for i_bin in range (0,etaEffs[1].GetPassedHistogram().GetNbinsX()):
            MC_eta.SetBinContent(i_bin,etaEffs[1].GetEfficiency(i_bin))
            print i_bin, " " , etaEffs[1].GetEfficiency(i_bin)
 #       MC_eta.Sumw2()                                                                                                                                                           

        # SF eta hist                                                                                                                                                             
        SF_eta = Data_eta.Clone("SF_eta"+histStr)
        # Data/MC                                                                                                                                                                
        SF_eta.Divide(MC_eta)


        # Draw TH1D                                                                                                                                                              
        c_sf=rt.TCanvas("c_sf"+histStr)
        c_sf.cd()
        SF_eta.Draw("PE")
        SF_eta.SetMaximum(3)
        SF_eta.SetMinimum(0)
        c_sf.Print("plots/SFEta"+histStr+".pdf")

        
        
        





    # eo loop over the indices


# save the comparison plot!
#c_ptComp.Print("plots/EffPtComp.pdf")

#eo loop over the channnel

outfile.Close()
