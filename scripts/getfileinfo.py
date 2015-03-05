import os, sys, subprocess
import ROOT as rt

for filename in os.listdir("./mergedsamples/") :
    #    print filename
    file = rt.TFile("mergedsamples/"+filename,"READ")
    #    file.ls()
    histo = file.Get("cutflow")
    print "sample ",filename," ",histo.GetBinContent(1)
#    for ibin in range(1,histo.GetNbinsX()+1) :
#        print histo.GetBinContent(ibin)

