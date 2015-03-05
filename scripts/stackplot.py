import os, sys
from array import array
import ROOT as rt
import tdrstyle, CMS_lumi
##############
# example pyroot loop for histogram making on output trees of Ntupler
# January 2015 by freya.blekman@cern.ch
#

tdrstyle.setTDRStyle()

d0range =[0.0,0.02,0.04,0.06,0.08,0.1,2.0]
CMS_lumi.lumi_13TeV = "1 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "PHYS14 samples, work in progress"
iPos = 11

outfile = rt.TFile("output_stackplot.root","RECREATE")
outfile.cd()
# book some histograms
h_mud0 = rt.TH1F("h_mud0","muon |d_{0}| (cm)",100,0,2) #5,array('f',d0range) )
h_eled0 = rt.TH1F("h_eled0","ele |d_{0}| (cm)",100,0,2) #5,array('f',d0range) )
h_mud0.SetXTitle("muon |d_{0}| (cm)")
h_mud0.SetYTitle("ele/mu pairs")
h_eled0.SetXTitle("electron |d_{0}| (cm)")
h_eled0.SetYTitle("ele/mu pairs")
h_isoele = rt.TH1F("h_isoele","electron PF iso (0.3)",40,0,2)
h_isomu = rt.TH1F("h_isomu","muon PF iso (0.3)",40,-1000,0)
h_isoele.SetXTitle("electron PF iso (0.3)")
h_isoele.SetYTitle("ele/mu pairs")
h_dr = rt.TH1F("h_dr","#Delta R(e, #mu)",21,0,7)
h_dr.SetXTitle("#Delta R(e, #mu)")
h_dr.SetYTitle("ele/mu pairs")






# using lorentz vectors as easy to calculate angles, pseudorapidity, etc
lvmu=rt.TLorentzVector()
lve=rt.TLorentzVector()




#done_QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8
#done_QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8

# for bookkeeping
lumivalue = 1000.0
samplesnamesfancy = ["Drell-Yan","t#bar{t}+jets","W+jets","W+jets","QCD","QCD","QCD","QCD","QCD","QCD","QCD","QCD","QCD"]
samplesaddtoleg=[1,1,1,0,1,0,0,0,0,0,0,0,0]
samples = ["dy","ttbar","wjetsplus","wjetsminus","qcd1","qcd2","qcd3","qcd4","qcd5","qcd6","qcd7","qcd8","qcd9"]
samplesxsecs = [2008.4,831.76,11811.4,8677.3,677300000.*0.007,866600000.*0.00044,164300000.*0.00816,21810000.*0.01522,2999000.*0.02424,3529000.*0.158,128500.*0.0406,185900000.*0.00272,3495000.*0.01255]
samplespresel = [2820473.0,25437856.0,699606.0,226439.0,1986513.0,4768929.0,3742583.0,3893676.0, 3468633.,1958930.0,999553.0,1730223.0,1999717.0,1000000.0]

pathToRootFile="/user/qpython/TopBrussels7X/CMSSW_7_2_1_patch1/src/TopBrussels/HToZZBachelorProjectNtupleMaker/FreyasNtuple"
samplesrootfiles=[pathToRootFile+"/DYJetsToLL_M-50_13TeV-madgraph-pythia8*.root",
                  pathToRootFile+"/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola*.root",
                  pathToRootFile+"/WplusToMuNu_CT10_13TeV-powheg-pythia8*.root",
                  pathToRootFile+"/WminusToMuNu_CT10_13TeV-powheg-pythia8*.root",
                  pathToRootFile+"/QCD_Pt-20to30_EMEnriched_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-30to50_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-50to80_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-80to120_MuEnrichedPt5_PionKaonDecay_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_170toInf_bcToE_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8*.root",
                  pathToRootFile+"/QCD_Pt_80to170_bcToE_Tune4C_13TeV_pythia8*.root"]
samplescolors=[4,2,rt.kGreen,rt.kGreen,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow,rt.kYellow]

stackplot_mu = rt.THStack("stackplot_mu",h_mud0.GetTitle())
stackplot_ele = rt.THStack("stackplot_ele",h_eled0.GetTitle())
stackplot_dr = rt.THStack("stackplot_dr",h_dr.GetTitle())
stackplot_eleiso = rt.THStack("stackplot_eleiso",h_isoele.GetTitle())
stackplot_muiso = rt.THStack("stackplot_muiso",h_isomu.GetTitle())

legend = rt.TLegend(0.6,0.7,0.9,0.9)
legend.SetFillColor(0)
legend.SetFillStyle(0)

for isam in range(0,len(samples)) :
    # the analysis structure see TTree/TChain description on rt.cern.ch
    ch = rt.TChain("tree","tree")
    # TChain accepts wildcards and ~infinite numbers of input files! So no need to merge files!
    ch.Add(samplesrootfiles[isam])
# very loud but useful to know what variables are stored in a tree... it prints them all
#ch.Print()


    nevents=ch.GetEntries()
    if nevents == 0 :
        continue


    h_mud0_work = h_mud0.Clone("h_mud0_"+samples[isam])
    h_mud0_work.SetFillColor(samplescolors[isam])
    h_mud0_work.SetLineColor(rt.kBlack)
    h_eled0_work = h_eled0.Clone("h_eled0_"+samples[isam])
    h_eled0_work.SetFillColor(samplescolors[isam])
    h_eled0_work.SetLineColor(rt.kBlack)
    h_dr_work = h_dr.Clone("h_dr_"+samples[isam])
    h_dr_work.SetFillColor(samplescolors[isam])
    h_isoele_work = h_isoele.Clone("h_isoele_"+samples[isam])
    h_isoele_work.SetFillColor(samplescolors[isam])
    h_isomu_work = h_isomu.Clone("h_isomu_"+samples[isam])
    h_isomu_work.SetFillColor(samplescolors[isam])


    weight=lumivalue*samplesxsecs[isam]/(samplespresel[isam])
    ii=0
    # start of loop over events
    for iev in ch:
        if ii % (nevents/50.) ==0 :
            print samples[isam]," ", ii, "/", nevents, " ,", (100*ii)/nevents, "%"
        ii+=1


    # loop over muons - fill in lorentz vector and fill some histograms
        for imu in range(0,iev.nMuons) :
        
            for iele in range(0,iev.nElectrons) :
                if iev.charge_muon[imu]*iev.charge_electron[iele]>0 :
                    continue
                #                print "iso: ",iev.pfIso_electron[iele] , " " ,  iev.pfIso_muon[imu]
                h_isoele_work.Fill(iev.pfIso_electron[iele],weight)
                h_isomu_work.Fill(iev.pfIso_muon[imu],weight)
                    #                if iev.pfIso_electron[iele]>0.15 :
                    #                    continue

                lvmu.SetPxPyPzE(iev.pX_muon[imu],iev.pY_muon[imu],iev.pZ_muon[imu],iev.E_muon[imu])
                lve.SetPxPyPzE(iev.pX_electron[iele],iev.pY_electron[iele],iev.pZ_electron[iele],iev.E_electron[iele])
                    #                if abs(iev.d0_electron[iele])>0.02 :
                h_mud0_work.Fill(abs(iev.d0_muon[imu]),weight)
                    #                if abs(iev.d0_muon[imu])>0.02 :
                h_eled0_work.Fill(abs(iev.d0_electron[iele]),weight)
#                h_mud0_work.Fill(abs(iev.d0_muon[imu]),weight)
                        
                if abs(iev.d0_electron[iele])>0.02 and abs(iev.d0_muon[imu])>0.02 :
                    h_dr_work.Fill(lvmu.DeltaR(lve),weight)

                     # end of event loop
    # now fill stack plots
    stackplot_ele.Add(h_eled0_work)
    stackplot_mu.Add(h_mud0_work)
    stackplot_dr.Add(h_dr_work)
    stackplot_eleiso.Add(h_isoele_work)
    stackplot_muiso.Add(h_isomu_work)
    if samplesaddtoleg[isam] == 1 :
        legend.AddEntry(h_eled0_work,samplesnamesfancy[isam],"f")
# end of sample loop



tc = rt.TCanvas()
tc.cd()

stackplot_mu.Draw()
#stackplot_mu.SetMinimum(0.0001)
stackplot_mu.GetHistogram().SetXTitle(h_mud0.GetTitle())
stackplot_mu.GetHistogram().SetYTitle("events")
tc.SetLogy()
CMS_lumi.CMS_lumi(tc, 4, iPos)
tc.cd()
tc.Update()
tc.RedrawAxis()
legend.Draw("same")
tc.Update()
tc.Print("muon_d0_stackplot.C")

tc1 = rt.TCanvas()
tc1.cd()

stackplot_ele.Draw()
#stackplot_ele.SetMinimum(0.0001)
stackplot_ele.GetHistogram().SetXTitle(h_eled0.GetTitle())
stackplot_ele.GetHistogram().SetYTitle("events")

tc1.SetLogy()

CMS_lumi.CMS_lumi(tc1, 4, iPos)

tc1.cd()
tc1.Update()
tc1.RedrawAxis()
legend.Draw("same")
tc1.Update()
tc1.Print("electron_d0_stackplot.C")

tc2 = rt.TCanvas()
tc2.cd()

stackplot_eleiso.Draw()
#stackplot_eleiso.SetMinimum(0.0001)
stackplot_eleiso.GetHistogram().SetXTitle(h_isoele.GetTitle())
stackplot_eleiso.GetHistogram().SetYTitle("events")

tc2.SetLogy()

CMS_lumi.CMS_lumi(tc2, 4, iPos)

tc2.cd()
tc2.Update()
tc2.RedrawAxis()
legend.Draw("same")
tc2.Update()
tc2.Print("electron_iso_stackplot.C")


tc12 = rt.TCanvas()
tc12.cd()

stackplot_muiso.Draw()
#stackplot_muiso.SetMinimum(0.0001)
stackplot_muiso.GetHistogram().SetXTitle(h_isomu.GetTitle())
stackplot_muiso.GetHistogram().SetYTitle("events")

tc12.SetLogy()

CMS_lumi.CMS_lumi(tc12, 4, iPos)

tc12.cd()
tc12.Update()
tc12.RedrawAxis()
legend.Draw("same")
tc12.Update()
tc12.Print("muon_iso_stackplot.C")


tc3 = rt.TCanvas()
tc3.cd()

stackplot_dr.Draw()
#stackplot_dr.SetMinimum(0.0001)
stackplot_dr.GetHistogram().SetXTitle(h_dr.GetTitle())
stackplot_dr.GetHistogram().SetYTitle("events")

tc3.SetLogy()

CMS_lumi.CMS_lumi(tc3, 4, iPos)

tc3.cd()
tc3.Update()
tc3.RedrawAxis()
legend.Draw("same")
tc3.Update()
tc3.Print("electron_muon_dr_stackplot.C")




outfile.Write()
