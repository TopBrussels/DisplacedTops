import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
from datetime import datetime
# libray to copy files
import shutil 


# Define time variable                                                           
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
hh = str(now.hour)
mn= str(now.minute)

# pick one of the two above                                                      
#date = dd+"_"+mm+"_"+yyyy+"_"+hh+"h"+mn+"min"
date = dd+"_"+mm+"_"+yyyy

channels = ["MuMu","ElEl"] 


for chan in channels:
    print "\nSearching list of sample used for ", chan, " channel!"
    # getting the appropriate xml file
    if "MuMu" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesMuMuV9.xml')
    elif "ElEl" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElElV9.xml')
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElMuV9.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()
    #tree = ET.ElementTree(file='../config/FullMcBkgdSamplesV9.xml')
    #tree = ET.ElementTree(file='../config/DataSamples.xml')
    #tree = ET.ElementTree(file='../config/DisplacedTopsSignal.xml')
    #tree = ET.ElementTree(file='../config/FullSamplesElElV9.xml')
    #tree = ET.ElementTree(file='../config/FullSamplesMuMuV9.xml')
    
    root = tree.getroot()
    datasets = root.find('datasets')
    
    
    print "found  "  + str(len(datasets)) + " datasets"
    
    # create new dir if not already existing
    if not os.path.exists("SubmitScripts/"+chan):
        os.makedirs("SubmitScripts/"+chan)
        
    if not os.path.exists("SubmitScripts/"+chan+"/"+date):
        os.makedirs("SubmitScripts/"+chan+"/"+date)
    
    
    # vector containing all the root file for a given dataset
    topTrees = []
    
    
    # loop over all the dataset with add="1"
    for d in datasets:
        if d.attrib['add'] == '1':
            print "found dataset to be added..." + str(d.attrib['name'])
            commandString = "./TreeMaker "+str(d.attrib['name'])+" "+str(d.attrib['title'])+" "+str(d.attrib['add'])+" "+str(d.attrib['color'])+" "+str(d.attrib['ls'])+" "+str(d.attrib['lw'])+" "+str(d.attrib['normf'])+" "+str(d.attrib['EqLumi'])+" "+str(d.attrib['xsection'])+" "+str(d.attrib['PreselEff'])
            topTrees = glob.glob(d.attrib['filenames'])
    
            N_file = 1
            # loop over all the root files and make one job per root file
            for f in range(0,len(topTrees)):
                filename="SubmitScripts/"+chan+"/"+date+"/submit_"+str(d.attrib['name'])+"_"+str(N_file)+".sh"
                # copy a skeleton file that set up the code environment, the wall time and the queue
                shutil.copyfile("submitSkeleton.sh", filename)
                # append to command to be run at the end of the skeleton
                outfile = open (filename, 'a')
    #            print  "dcap://maite.iihe.ac.be:"+topTrees[f]
                print >> outfile, commandString, "dcap://maite.iihe.ac.be:"+topTrees[f], " ", chan , " " , str(N_file) , " 0" , " 2000000" 
                N_file=N_file+1
                
# moving the newly created dir
os.chdir("SubmitScripts/"+chan+"/"+date)
