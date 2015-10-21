from ROOT import TChain
from glob import glob
import xml.etree.cElementTree as ET
import os

# get filenams for the xml!!!
tree = ET.ElementTree(file='config/FullMcBkgdSamplesV8.xml')

root = tree.getroot()
datasets = root.find('datasets')

topTrees = []

for d in datasets:
    if d.attrib['add'] == '1':
        topTrees.append(d.attrib['filenames'])

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
    print "***************************"
    print "****End of sample *********"    
    print "***************************"
    print "\n"


