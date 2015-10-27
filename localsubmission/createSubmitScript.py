import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
# libray to copy files
import shutil 



# getting the appropriate xml file
tree = ET.ElementTree(file='../config/FullMcBkgdSamplesV8.xml')
#tree = ET.ElementTree(file='../config/DataSamples.xml')


root = tree.getroot()
datasets = root.find('datasets')

print "found  "  + str(len(datasets)) + " datasets"

# create new dir if not already existing
if not os.path.exists("SubmitScripts"):
    os.makedirs("SubmitScripts")


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
            filename="SubmitScripts/submit_"+str(d.attrib['name'])+str(N_file)+".sh"
            # copy a skeleton file that set up the code environment, the wall time and the queue
            shutil.copyfile("submitSkeleton.sh", filename)
            # append to command to be run at the end of the skeleton
            outfile = open (filename, 'a')
            print >> outfile, commandString, topTrees[f], " 0" , " 2000000"
            N_file=N_file+1
            
