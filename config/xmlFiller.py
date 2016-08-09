#!/usr/bin/env python

# python script for filling eqlumi into config files

import shutil 
import os
import sys
from ROOT import TChain
import glob
import xml.etree.cElementTree as ET

inputFile = sys.argv[1]
#inputFile1=split(inputFile,.)

output = "Filled_"+inputFile


LumiTot=0

tree = ET.ElementTree(file=inputFile)
root = tree.getroot()
datasets = root.find('datasets')
for d in datasets:
    if d.attrib['add'] == 1:
        continue
    # calculate the sum of the lumi for the different data run 
    if 'data' in d.attrib['title'].lower() and 'run' in d.attrib['name'].lower():
        print "Eqlumi is " , d.attrib['EqLumi']
        LumiTot=LumiTot+float(d.attrib['EqLumi'])
        d.set('add','0')
        continue
    # run only on the merged data dataset with the sum of the lumi
    if 'data' in d.attrib['title'].lower() and not 'run' in d.attrib['name'].lower():
        d.set('EqLumi',str(LumiTot)) 
        d.set('add','1')
        continue
    xsec = float(d.attrib['xsection'])
    files = glob.glob(d.attrib['filenames'])
    root_files = []
    for f in files:
        root_files.append('dcap://maite.iihe.ac.be'+f)
    chain = TChain('eventTree')
    for rf in root_files:
        chain.Add(rf)
    nEntries = chain.GetEntries()
    equivLumi = nEntries/xsec
    d.set('EqLumi',str(equivLumi))
    print d.attrib['name']
print 'filled xml with eqlumis!'
tree.write(output)


# removes the <data> from the file
#with open(tempOutputFileName,"r") as input:
#    with open(outputFileName,"wb") as output:
#        for line in input:
#            if "data>" not in line:
#                output.write(line)


