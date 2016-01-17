from ROOT import TChain
import ROOT
import glob
import xml.etree.cElementTree as ET
import os
from datetime import datetime


# Define time variable
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
# pick one of the two above
date = dd+"_"+mm+"_"+yyyy


#date = "1_12_2015"
#channel = ["_MuEl","_MuMu","_ElEl"]
channel = ["_MuMu"]


for chan in channel:
    
    #Define path where ntuples are stored
    pathNonMerged = "MACRO_Output"+chan+"/"  #needs to be changed for different lepton channel
    pathMerged = "MergedTrees/"+date+"/"+chan+"/"
    
    if not os.path.exists(pathMerged):
        os.makedirs(pathMerged)
    
    
    # get filenames from the xml!!!
    #tree = ET.ElementTree(file='config/FullMcBkgdSamplesV9.xml')
    tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')
    #tree = ET.ElementTree(file='config/DataSamples.xml')
    
    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    
    # loop over the datasets to be added and fill the "topTrees" vector
    for d in datasets:
        if d.attrib['add'] == '1':
            print "found dataset to be added..." + str(d.attrib['name'])
            # select a subset of the existing root file
            if "" in str(d.attrib['name']) :
                datasetNames.append(str(d.attrib['name']))
                print str(d.attrib['name'])
    
    
    listOfZombie= []
    
    for n in datasetNames:
        filenames = glob.glob(pathNonMerged + "/*" + n + chan + "*.root")
        hadd = "hadd " + pathMerged + "DisplacedTop_Run2_TopTree_Study_" + n + chan + ".root"
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
    if (len(listOfZombie) > 0):
        print "The total number of zombie file is ", len(listOfZombie)
        print "And the list of the zombie is :"
        for zombie in listOfZombie:
            print zombie
    
    mergeData=False
    
    if (mergeData):
    # combining all the Data in one
        dataList=glob.glob(pathMerged+"*Data*.root")
    
        cmd = "hadd " + pathMerged + "DisplacedTop_Run2_TopTree_Study_Data"+chan + ".root"
        for data in dataList:
            cmd = cmd + " " + data
        os.system(cmd)
            
