from ROOT import TChain
import ROOT
import sys
import glob
import xml.etree.cElementTree as ET
import os
from datetime import datetime


# Define time variable
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
date = dd+"_"+mm+"_"+yyyy
#date = "24_11_2016"

#directory name
directory = "NoBlindingRerun"+"_"+date


#usging argument to filter
filterSample = sys.argv[1]



# config
sampleVersion="V4"


# pick the list of channels to be merged
#channels = ["_MuMu","_ElEl","_bbMu","_bbEl","_ttElEl","_ttMuMu"]
#channels= ["_bbMu","_bbEl"]
#channels = ["_bbMu"]
#channels = ["_bbEl"]
channels = ["_MuMu","_ElEl"]
#channels = ["_ElEl"]
#channels = ["_MuMu"]
#channels = ["_ttElEl","_ttMuMu"] 
#channels = ["_ttMuMu"] 
#channels = ["_Fail"]





for chan in channels:


    # set the btag WP filter
    btagFilters=[""]
    
    if "_bb" in chan:
        btagFilters=["_Loose", "_Medium", "_Tight"]
        
    
    

    # loop over the btag WP
    for btagFilter in btagFilters:

        print "btagFilter is " , btagFilter


        #Define path where ntuples are stored
        pathNonMerged = "MACRO_Output"+chan+"/"  
        pathMerged = "MergedTrees/"+directory+"/" + btagFilter + chan +"/"
        
        # create dir if not existing yet
        if not os.path.exists(pathMerged):
            os.makedirs(pathMerged)
        
        # remove "_" from the chann string
        chanStrip=chan.replace("_","")
    
        # get the correct xml file !!!    
        if "ElEl" == chanStrip or "MuMu" == chanStrip or "ElMu" == chanStrip or "bbEl" == chanStrip or "bbMu" == chanStrip:
            tree = ET.ElementTree(file='config/'+chanStrip+sampleVersion+'.xml')
        elif "ttElEl" in chanStrip or "ttMuMu" in chanStrip:
            tree = ET.ElementTree(file='config/'+'ttLeptons'+sampleVersion+'.xml')
        else:
            print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
            print "you have tried to load the config ",  "config/'+chanStrip+sampleVersion+'.xml"
            sys.exit()
    
        #tree = ET.ElementTree(file='config/FullMcBkgdSamplesV9.xml')
        #tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')
        #tree = ET.ElementTree(file='config/DataSamples.xml')
    
        # get the list of dataset
        root = tree.getroot()
        datasets = root.find('datasets')
        print "found  "  + str(len(datasets)) + " datasets"
        datasetNames = []
        
        print ""
        # loop over the datasets to be added and fill the "topTrees" vector
        for d in datasets:
            if d.attrib['add'] == '1':
                print "found dataset to be added..." + str(d.attrib['name'])
    
                # select a subset of the existing root file
                if filterSample in str(d.attrib['name']) :
                    datasetNames.append(str(d.attrib['name']))
                    print str(d.attrib['name'])
        
        
        listOfZombie= []
        
        # loop over data set to search root files
        for n in datasetNames:
            filenames = glob.glob(pathNonMerged + "/*" + n + btagFilter + chan + "*.root")
            
#            print filenames
            
     #       ext=""
     #       if "NoBlinding" in filenames[1]:
      #          ext="_NoBlinding"
      #          print ext
                
            hadd = "hadd " + pathMerged + "DisplacedTop_Run2_TopTree_Study_" + n + filterSample + chan + ".root"
    
            if (len(filenames) == 0):
                print "no root files found in directory" , pathNonMerged ,  " for dataset " , n , " !!"
            else :
                # loop over root files
                for f in filenames:
                    file=ROOT.TFile(f,"read")
                    # check if the file is a zombie
                    if (file.IsZombie()):
                        print "File" , f, "is a Zombie.... Skipping"
                        listOfZombie.append(f)
                    else:
                        print f
                        hadd = hadd + " " + f
                print "Merging ntuples for " + n
                os.system(hadd)
            
        print "\n\n"
        
        # print the list of zombies
        print "The total number of zombie file is ", len(listOfZombie)
        if (len(listOfZombie) > 0):
            outfile = open (pathMerged+"/Zombie"+chan+".txt", 'a')
            print "And the list of the zombie is put in "+pathMerged+"/Zombie"+chan+".txt "
            for zombie in listOfZombie:
                print >> outfile, zombie
        
        mergeData=False
        
        if (mergeData):
        # combining all the Data in one
            dataList=glob.glob(pathMerged+"*Data*.root")
        
            cmd = "hadd " + pathMerged + "DisplacedTop_Run2_TopTree_Study_Data"+chan + "*.root"
            for data in dataList:
                cmd = cmd + " " + data
            os.system(cmd)
                
    
        renameData = False
        # rename the data sample
        if (renameData):
            mvcmd ="mv "+pathMerged+"DisplacedTop_Run2_TopTree_Study_DataRunD"+chan+".root"+" "+ pathMerged+"DisplacedTop_Run2_TopTree_Study_Data"+chan+".root"
            os.system(mvcmd)
    

