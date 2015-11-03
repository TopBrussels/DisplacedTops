import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob

#tree = ET.ElementTree(file='config/Run2DiLepton_TOPTREES.xml')
#tree = ET.ElementTree(file='config/TestFaco.xml')
#tree = ET.ElementTree(file='config/FullMcBkgdSamples.xml')
#tree = ET.ElementTree(file='config/FullMcBkgdSamplesV7.xml')
#tree = ET.ElementTree(file='config/FullMcBkgdSamplesV8.xml')
#tree = ET.ElementTree(file='config/DataSamples.xml')
#tree = ET.ElementTree(file='config/QCDSamplesV6.xml')
tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')

root = tree.getroot()
datasets = root.find('datasets')

print "found  "  + str(len(datasets)) + " datasets"

procsDone = 0
procsStarted = 0
numCores = 8
args = []
execCommands = []
topTrees = []
jobSize = 2000000
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        files = ["./TreeMaker", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff']]
        topTrees = glob.glob(d.attrib['filenames'])
        for f in glob.glob(d.attrib['filenames']):
            files.append("dcap://maite.iihe.ac.be"+f)
        args.append(files)
outfiles = []
fileNames = []
processes = []
tempList = []
if not os.path.exists("Terminal_Output"):
    os.makedirs("Terminal_Output")

for row in args:
    print "checking args..."
    if row[3] == '1':
        title = row[1]

        totalEvents = float(row[8])*float(row[9])
#        tempList = list(row)
#        tempList.extend(["", ""])
        if (totalEvents > jobSize):
            endEvent = 0
            while (totalEvents-(endEvent*jobSize) > 0):
                startStr = str(endEvent*jobSize)
                endStr = str((endEvent+1)*jobSize)
                tempList = list(row)
                tempList.extend(["", ""])
                tempList[len(tempList)-2] = startStr
                tempList[len(tempList)-1] = endStr
                tempList[1] = title+"_"+str(endEvent+1)
                fileNames.append("Terminal_Output/"+tempList[1]+".out")
                execCommands.append(tempList)
                print 'Job {} Created'.format(tempList[1])
                #popen = subprocess.Popen(execCommands, stdout=outfile)
                #processes.append(popen)
                endEvent += 1
        else:
            tempList = list(row)
            tempList.extend(["", ""])
            tempList[len(tempList)-2] = "0"
            tempList[len(tempList)-1] = str(jobSize)
            fileNames.append("Terminal_Output/"+tempList[1]+".out")
            execCommands.append(tempList)
            #popen = subprocess.Popen(execCommands, stdout=outfile)
            #processes.append(popen)
#        popen.wait()
#        for i in row:
#            print i
for i, row in enumerate(execCommands):
    outfile = open(fileNames[i], 'w')
    print "file name  = " + str(fileNames[i])
    print 'File {} opened'.format(fileNames[i])
    outfiles.append(outfile)
    row.insert(0, "nohup")

    popen = subprocess.Popen(row, stdout=outfiles[i])
    # faco get name
#    print "row is ", row
#    print "stdout is ", outfiles[i] 

    print 'Job {} begun'.format(row[2])
    processes.append(popen)
    procsStarted += 1
    print 'Jobs {} of {} started.  Timestamp: {}'.format(procsStarted, len(execCommands), time.ctime())
    while (procsStarted-procsDone) >= (numCores/2):
        time.sleep(60)
        procsDone = 0
        for proc in processes:
            if proc.poll() != None:
                procsDone+= 1
        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())
while (procsDone != len(execCommands)):  #This loop controls the status output for the last 4 jobs that are still running when the above for loop terminates
        time.sleep(60)
        procsDone = 0
        for proc in processes:
            if proc.poll() != None:
                procsDone+= 1
        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())



#while procsDone < len(processes):
#    time.sleep(60)
#    procsDone = 0
#    counter = 0
#    for proc in processes:
#        counter += 1
#        if proc.poll() != None:
#            procsDone += 1
#            print 'Job {} complete.'.format(counter)
#        else:
#            print 'Job {} still running.'.format(counter)
#    print '{} jobs of {} completed.  Timestamp: {}'.format(procsDone, len(processes), time.ctime())

