from ROOT import TChain
from glob import glob
import xml.etree.cElementTree as ET
import os

# get filenames from the xml!!!
#tree = ET.ElementTree(file='config/FullMcBkgdSamplesV7.xml')
tree = ET.ElementTree(file='config/FullMcBkgdSamplesV8.xml')
#tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')
#tree = ET.ElementTree(file='config/DataSamples.xml')
#tree = ET.ElementTree(file='config/FullMcBkgdSamples.xml')


root = tree.getroot()
datasets = root.find('datasets')

topTrees = []
xsections = []

# loop over the datasets to be added and fill the "topTrees" vector
for d in datasets:
    if d.attrib['add'] == '1':
        topTrees.append(d.attrib['filenames'])
        xsections.append(float(d.attrib['xsection']))



# loop over the "topTrees" vector
for n_sample in range(0,len(topTrees)):

    path = topTrees[n_sample]
    print path
    files = glob(path)
    root_files = []
    for f in files:
    	root_files.append('dcap://maite.iihe.ac.be' + f)
    #print root_files
    chain = TChain('eventTree')
    for rf in root_files:
    	chain.Add(rf)
    print 'added files'
    nEntries = chain.GetEntries();
    print "\n"
    print topTrees[n_sample], " contains ", nEntries, " events!"
    if (xsections[n_sample] != 1):
        print " xsections is " ,  xsections[n_sample]
        print "Eqlumi is N/xsect --->" , nEntries , "/" , xsections[n_sample] , " = " , nEntries/xsections[n_sample] 
    print "***************************"
    print "****End of sample *********"    
    print "***************************"
    print "\n"


