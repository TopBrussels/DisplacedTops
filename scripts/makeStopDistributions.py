"""
This- script will produce fancy plots for some gen-level quantities of the stops.

"""



# basic import
import os, sys
import ROOT
import math



# path to root files
pathToRootFile="../MACRO_Output_MuEl/"


# vector of samples containing the mass and the ctau and the colour
samples=[[200,1,2],[500,10,3],[800,100,4]]


# vector of histos
histos=[]


# loop over some samples 
for sample in samples :

    Mass=sample[0]
    Ctau=sample[1]
    print "new sample with Mass = " , Mass, "and Ctau = ", Ctau
    
    
    ch  = ROOT.TChain("preCutTree","preCutTree")
    ch.Add(pathToRootFile+"DisplacedTop_Run2_TopTree_Study_stopTobl_m" + str(Mass) + "_Ctau" + str(Ctau) + "_MuEl_NoBlinding_1.root")


#        ch.Add(pathToRootFile+"DisplacedTop_Run2_TopTree_Study_NP_overlay_stopTobl_m500_Ctau10_ElEl_NoBlinding_ZPeak_1.root")
        


        # booking histograms
    h_stop_pt = ROOT.TH1F("h_stop_pt","stop p_{T}",70,0,1400)
    h_stop_eta = ROOT.TH1F("h_stop_eta","stop #eta",100,-5,5)
    h_stop_ctau = ROOT.TH1F("h_stop_d0","stop c#tau",100,0,100)
    h_stop_phi = ROOT.TH1F("h_stop_phi","stop #phi",32,-3.2,3.2)
    


        

    ii=0
        # loop over the events
    for iev in ch:
        #            print ii
        if 100 < ii :
            continue

        # loop over the gen particles
        for i_gen in range(0,iev.nMcParticles_pc) :
            #                print i_gen
            
            # check if a stop
            if abs(iev.type_mcParticle_pc[i_gen]) ==  1000006:
                #                    Print "found a stop"

                # calculate the flight distance (c*tau)
                ctau = math.sqrt (iev.v0_mcParticle_pc[i_gen]**2 + iev.vz_mcParticle_pc[i_gen]**2)

                    # fill the histo
                h_stop_pt.Fill(iev.pt_mcParticle_pc[i_gen])
                h_stop_eta.Fill(iev.eta_mcParticle_pc[i_gen])
                h_stop_ctau.Fill(ctau)
                h_stop_phi.Fill(iev.phi_mcParticle_pc[i_gen])
                

                


        ii+=1
    # eo loop over the events     

    
    histos.append([h_stop_pt, h_stop_eta, h_stop_ctau, h_stop_phi])
    print "histos is ",histos
# eo over the sample



# loop over the samples
for i_sample in range (0,len(samples)):


    # loop over the histo ~ the variables
    for i_histo in range (0,len(histos[i_sample])):

        # one canvas per variable
        canvas = ROOT.TCanvas("canvasOverlay"+str(i_histo))
        canvas.cd()


        histos[i_sample][i_histo].SetLineColor(samples[i_sample][2])

        if i_sample == 0:
            histos[i_sample][i_histo].Draw()
        else :
            histos[i_sample][i_histo].Draw("SAME")
    
    # eo loop over the histo ~ the variables



    # save one canvas

#eo loop over the samples 



canvasOverlay.SaveAs("overlay.pdf")
    

canvas = ROOT.TCanvas("canvas")
canvas.Divide(2,2)
canvas.cd(1)
h_stop_pt.Draw()

canvas.cd(2)
h_stop_eta.Draw()

canvas.cd(3)
h_stop_ctau.Draw()

canvas.cd(4)
h_stop_phi.Draw()

canvas.SaveAs("yo.pdf")
