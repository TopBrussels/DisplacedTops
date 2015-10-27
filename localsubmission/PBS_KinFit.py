# python script to run the KinFits (TopTreeAnalysis) on PBS.

# for command line options
from optparse import OptionParser

# interacting with the os
from subprocess import Popen, PIPE, STDOUT
import sys
import os
import os.path

# working with time
import time
from time import strftime
from datetime import datetime

# import packages for multi-threading
import Queue
import threading

# ------------- #
#  LogHandler   #
# ------------- #

class logHandler:

    def __init__ (self, fileName):
        self.logFile = fileName
        if not self.logFile == "" and os.path.exists(self.logFile):
            os.remove(self.logFile)

    def output(self, string):
        global f
        if not self.logFile == "":
            f = open(self.logFile, "a")
            f.write("\n["+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"] "+string+"\n")
            f.close()
        else:
            print "\n["+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"] "+string

# --------------- #
#   MailHandler   #
# --------------- #

# importing smtp lib
import smtplib


class MailHandler:

    def __init__(self, recepient):
        self.smtpServer = "mach.vub.ac.be"
        #self.smtpServer = "localhost"
        self.senderAddress = "PBSKinFitter@mtop.iihe.ac.be"
        #+Popen('hostname', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().strip()
        #self.toAnnounce = [ "top-brussels-datasets@cern.ch" ]
        self.toAnnounce = recepient.split(',')

    def sendMail(self, subject,msg):
        toAddrs = ""
        for to in range(0,len(self.toAnnounce)):
            toAddrs = toAddrs+self.toAnnounce[to]+", "
        m = "From: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\n\r\n" % (self.senderAddress, toAddrs, subject)
        server = smtplib.SMTP(self.smtpServer)
        server.sendmail(self.senderAddress, toAddrs.split(), m+msg)
        server.quit()

# --------------- #
#  KinFitHandler  #
# --------------- #


class KinFitHandler:

    def __init__(self, nJob, name, systematic, crossSection, intLumi, pathPNFS, inputFileNr):  # , monsterFile):

        self.nJob = nJob
        self.pbsFile = ""
        self.pbsLog = ""
        self.pbsID = ""
        self.taskName = ""

        self.inputName = name
        self.systematicOption = systematic
        self.inputXS = crossSection
        self.inputLumi = intLumi
        self.inputPNFSDir = pathPNFS
        self.inputFileNr = inputFileNr
        #self.monsterFile = monsterFile        # if empty, run kinFits, otherwise run on the monster
        #self.monsterFileName = ""
        #if not self.monsterFile == "":
        #    splitted = self.monsterFile.split("/")
        #    self.monsterFileName = splitted[len(splitted)-1]

        self.xmlCFG = ""
        self.workingDir = ""
        self.resultsDir = ""
        self.plotsMacroDir = ""
        self.log = ""

    def setlog(self, log):
        self.log = log

    def setupWorkingDir(self):

        global options
        global timestamp
        global userName

        self.resultsDir = "./Results/RESULTS_"+options.TaskName+"_"+timestamp+"/"
        self.plotsMacroDir = "./Results/PLOTSMACRO_"+options.TaskName+"_"+timestamp+"/"

        if not os.path.exists(self.resultsDir) and self.nJob == 0:
            os.mkdir(self.resultsDir)

        if not os.path.exists(self.plotsMacroDir) and self.nJob == 0:
            os.mkdir(self.plotsMacroDir)

        if not options.local:
            self.workingDir = "/localgrid/"+userName+"/"+options.TaskName+"_"+timestamp+"_job_"+str(self.nJob)
            if not os.path.exists("/localgrid/"+userName+"/LightTree_PBSScript/"+options.TaskName+"_"+timestamp):
                os.mkdir("/localgrid/"+userName+"/LightTree_PBSScript/"+options.TaskName+"_"+timestamp)
        else:
            print "NOT YET IMPLEMENTED LOCAL RUNNING"
            sys.exit(1)

        cmd = "cp -vfr "+options.WorkingDir+" "+self.workingDir
        Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        # self.log.output("File:  |"+self.monsterFile+"|  FileName:  |"+self.monsterFileName+"|")

        # if not self.monsterFile == "":
        #    self.log.output("Copying file:  "+self.monsterFile+"  to  "+self.workingDir)
        #    self.log.output("  --> Output of command:  "+Popen("cp -vfr "+self.monsterFile+" "+self.workingDir+"/", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

    def createXMLCFG(self):

        self.log.output("Applying command line option:  "+self.systematicOption)
        self.log.output("Creating xml config files for TopTreeAnalysis")
        self.xmlCFG = options.WorkingDir+"/myAnomCouplconfig_"+str(self.nJob)+".xml"

        xmlFile = open(self.xmlCFG,"w")
        xmlFile.write("<?xml version=\"1.0\"?>\n")
        xmlFile.write("<datasets>\n")

        dataSetXML = "<d name=\""+self.inputName+"\" title=\""+self.inputName+"\" add=\"1\" color=\"1\" ls=\"1\" lw=\"2\" normf=\"1\" "
        dataSetXML += "xsection=\""+self.inputXS+"\" EqLumi=\""+self.inputLumi+"\" filenames=\"TopTree_Skimmed_*.root\"/>"
        self.log.output(dataSetXML)
        xmlFile.write(dataSetXML+"\n")

        xmlFile.write("</datasets>\n\n")
        xmlFile.write("<analysis>\n")

        xmlFile.write("<a type=\"Collections\" PVCollection=\"PrimaryVertex\" JetType=\"2\" JetCollection=\"PFJets_selectedPatJetsPF2PAT\" METType=\"2\" METCollection=\"PFMET_patType1CorrectedPFMetPF2PAT\" MuonCollection=\"Muons_selectedPatMuonsPF2PAT\" ElectronCollection=\"Electrons_selectedPatElectronsPF2PAT\" loadGenJetCollection=\"1\" GenJetCollection=\"GenJets_ak5GenJetsNoNu\" loadGenEventCollection=\"1\" GenEventCollection=\"GenEvent\" loadNPGenEventCollection=\"0\" NPGenEventCollection=\"NPGenEvent\" loadMCParticles=\"1\" MCParticlesCollection=\"MCParticles\" TrackMETCollection=\"\" loadTrackMET=\"0\"/>\n")

        xmlFile.write("<a type=\"Selection\" PVertexNdofCut=\"4\" PVertexZCut=\"24.\" PVertexRhoCut=\"2.\" MuonPtCutSR=\"20.\" MuonEtaCutSR=\"2.1\" MuonRelIsoCutSR=\"0.05\" MuonNHitsCutSR=\"10\" MuonD0CutSR=\"0.02\" MuonDRJetsCut=\"0.3\" MuonPtCutVetoSR=\"10.\" MuonEtaCutVetoSR=\"2.5\" MuonRelIsoCutVetoSR=\"0.2\" ElectronPtCut=\"15.\" ElectronEtaCut=\"2.5\" ElectronRelIsoCut=\"0.2\" JetsPtCutSR=\"30.\" JetsEtaCutSR=\"2.4\" applyJetID=\"1\" JetEMFCut=\"0.01\" n90HitsCut=\"1\" fHPDCut=\"0.98\" NofJets=\"4\" NofJetBins=\"2\"/>\n")
        xmlFile.write("<a type=\"Conditions\" isMC=\"1\" MCRound=\"0\" Vars_ByFile=\"0\" VarsFile=\"m0_100_m12_100\" IntToCut=\"4\" Verbose=\"2\" Luminosity=\"9999999\" JES=\"1.\" nPseudoExp=\"0\" nPseudoSession=\"0\" runonTTrees=\"0\" doABCD=\"1\" doVJEstim=\"1\" doVJEstPE=\"1\" doTtJEstim=\"1\" doTemplComp=\"0\" doSystematics=\"0\"/>\n")
        xmlFile.write("<a type=\"CRForTtbarEstimation\" BtagAlgo_ttjEst=\"0\" BtagDiscriCut_ttjEst=\"4.38\" MuonPtCutCR=\"30.\" MuonEtaCutCR=\"2.1\" MuonRelIsoCutCR=\"0.1\" JetsPtCutCR=\"30.\" JetsEtaCutCR=\"2.4\" MblCut=\"160.\" DRBBCut=\"2.3\" HTBBCut=\"500.\" NREvtFraction=\"0.75\"/>\n")
        xmlFile.write("<a type=\"CRForABCDEstimation\" NXbinsABCD=\"200\" NYbinsABCD=\"200\" XbinMinABCD=\"0\" XbinMaxABCD=\"20\" YbinMinABCD=\"0\" YbinMaxABCD=\"20\" cutXmin=\"0\" cutX0=\"0.1\" cutX1=\"0.2\" cutXmax=\"20.\" cutYmin=\"0.\" cutY0=\"3.\" cutY1=\"4.\" cutYmax=\"20.\" region=\"1\"/>\n")
        xmlFile.write("<a type=\"ParamForVJetEstimation\" BtagAlgo_vjEst=\"0\" NofBtagWorkingPoint_vjEst=\"1\" BtagWorkingPoint_vjEst=\"2.03,3.20\" MinMethod=\"Minuit2\" MinOption=\"Combined\" useMJLE=\"0\" useUnBinMLE=\"1\" NVJetPE=\"500\" TagEffInit=\"0.794,0.128,0.097-0.70,0.043,0.02-0.63,0.05,0.010/0.807,0.134,0.124-0.70,0.043,0.02-0.63,0.05,0.010\" NVlikeInit=\"14./4.\" NTTlikeInit=\"6./8.\" EffEbsel=\"0.0515,0.4170,0.5281/0.0187,0.2604,0.7049\" VJEstFixParam=\"0,1,2\" NofIterationsVJestShapeEstim=\"40\"/>\n")
        xmlFile.write("<a type=\"Observables\" runOnObsByString=\"0\" listOfObsInts=\"2\" listOfObsStrings=\"ET1oET2,ET1oET3\" binning=\"../config/Binning.root\" bins=\"20\"/>\n")
        xmlFile.write("<a type=\"CrossSection\" MCExpFilename=\"../config/MCFile.root\" LuminosityError=\"0.1\" TriggerEff=\"1\" TriggerEffError=\"0.05\" SkimEff=\"1.\" SkimEffError=\"0\" MuonSelEff=\"0.43\" MuonSelEffError=\"0.003\" SecondLeptonVetoEff=\"0.4833\" SecondLeptonVetoEffError=\"0.01\" JetSelEff=\"0.8206\" JetSelEffError=\"0.04\" NofSingleTopEvts=\"0\" NofSingleTopEvtsError=\"0\"/>\n")
        xmlFile.write("<a type=\"Search\" doBkgEstim=\"1\" doDumpPseudoExpInfoInTTree=\"1\" DumpTreeName=\"dumpTreeFile.root\" FractionHWEvts=\"0.1,0.2,0.3\"/>\n")
        xmlFile.write("</analysis>")
        xmlFile.close()

    def createPBSCFG (self):

        global timestamp
        global userName

        self.pbsFile = "/localgrid/"+userName+"/"+options.TaskName+"_"+timestamp+"_job_"+str(self.nJob)+".pbs"
        self.pbsLog = "/localgrid/"+userName+"/"+options.TaskName+"_"+timestamp+"_job_"+str(self.nJob)+".out"
        self.taskName = options.TaskName+"_"+timestamp+"_Job"+str(self.nJob)

        pbs = open(self.pbsFile, "w")

        pbs.write("#! /bin/bash\n")
        pbs.write("#PBS -r n\n")
        pbs.write("#PBS -N "+self.taskName+"\n")
        pbs.write("#PBS -j oe\n")
        pbs.write("#PBS -k oe\n")  # dit zorgt voor realtime log uitspuwen
        # if self.monsterFile == "":
        pbs.write("#PBS -l walltime=20:00:00\n")
        # else:
        #    pbs.write("#PBS -l walltime=06:00:00\n")

        pbs.write("echo \"dumping some info off the worker node\"\n")
        pbs.write("hostname\n")
        pbs.write("df -h\n")
        pbs.write("uptime\n")
        pbs.write("free\n")
        pbs.write("ls -l /scratch/\n")
        pbs.write("export ROOTSYS=/localgrid/"+userName+"/root\n")
        pbs.write("export PATH=$ROOTSYS/bin:$PATH\n")
        pbs.write("export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH\n")
        pbs.write("export LD_LIBRARY_PATH=/scratch/$PBS_JOBID/AnomCoup_setup/lib:$LD_LIBRARY_PATH\n")
        pbs.write("echo \"Root Version: $(root-config --version)\"\n")

        pbs.write("\nmv -f "+self.workingDir+" /scratch/$PBS_JOBID/AnomCoup_setup\n")
        pbs.write("cd /scratch/$PBS_JOBID/AnomCoup_setup\n")
        pbs.write("ls -l PersonalClasses/Calibrations/JECFiles/\n")
        # pbs.write("rm -rf Monsters/*\n")                                 --> Why do you want to remove all the previous ones??
        # pbs.write("rm -rf LikelihoodResults_ASCII/*\n")
        pbs.write("echo \"Downloading files from "+self.inputPNFSDir+"\"\n")

        # if self.monsterFile == "":
        cmd = "ls "+self.inputPNFSDir+" | grep TopTree_Skimmed_ | grep root"
        inputFiles = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().split("\n")
        del inputFiles[-1]
        pbs.write("echo \"inputFiles:  "+str(inputFiles)+"\"\n")

        if self.inputFileNr == -1:
            for i in inputFiles:
                pbs.write("echo \"Downloading file:  "+i+"\"\n")
                pbs.write("dccp dcap://maite.iihe.ac.be"+self.inputPNFSDir+i+" /scratch/$PBS_JOBID/AnomCoup_setup/"+str(i)+"\n")
        else:
            pbs.write("echo \"Downloading file nr:  "+str(self.inputFileNr-1)+"\"\n")
            pbs.write("echo \"Downloading file:  "+inputFiles[self.inputFileNr-1]+"\"\n")
            pbs.write("dccp dcap://maite.iihe.ac.be"+self.inputPNFSDir+inputFiles[self.inputFileNr-1]+" /scratch/$PBS_JOBID/AnomCoup_setup/TopTree_Skimmed_"+str(self.inputFileNr)+".root\n")

        pbs.write("ls -ltr /scratch/$PBS_JOBID/AnomCoup_setup\n\n")

        # if self.monsterFile == "":
        # pbs.write("./AnomalousCouplingsTreeCreator "+self.systematicOption+" myJESconfig_"+str(self.nJob)+".xml\n")
        pbs.write("mv myAnomCouplconfig_"+str(self.nJob)+".xml myAnomCouplconfigUsed.xml\n")
        pbs.write("rm myAnomCouplconfig_*.xml\n")
        pbs.write("./AnomalousCouplingsTreeCreator myAnomCouplconfigUsed.xml\n")    # _"+str(self.nJob)+".xml\n")
        # pbs.write("rm myAnomCouplconfig_"+str(self.nJob)+".xml\n")
        pbs.write("rm -f TopTree_Skimmed_*.root\n")
        pbs.write("ls -l LightTree/\n")
        #else:
        #    pbs.write("./MTopDiff_Analysis "+self.monsterFileName+"\n")
        #    pbs.write("rm -rf *.root Monsters "+self.monsterFileName+"\n")
        #    pbs.write("ls -lah LikelihoodResults_ASCII/\n")

        #pbs.write("rm -rf resolutions *.so FitResults_ASCII PileUpReweighting Plots* weights JECFiles\n")
        pbs.write("rm -rf AnomalousCouplingsTreeCreator \n")
        pbs.write("mv -fv /scratch/$PBS_JOBID/AnomCoup_setup "+self.workingDir+"\n")
        pbs.write("echo \"THIS IS THE END\"\n")
        pbs.close()

    def submitPBSJob(self):

        global timestamp
        global userName

        cmd = "cd /localgrid/"+userName+"/; qsub -q localgrid@cream02.iihe.ac.be "+self.pbsFile

        #print cmd
        self.pbsID=Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        self.log.output("PBS Job ID: "+self.pbsID)

    def checkPBSJob (self):

        global timestamp
        global userName
        global SleepTime

        status = Popen("qstat "+self.pbsID, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

        while status.find("Unknown Job Id") == -1:

            self.log.output("     -----> It seems that the job is still running, sleeping "+str(SleepTime)+"s")
            time.sleep(SleepTime)
            status = Popen("qstat "+self.pbsID, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

        self.log.output("-> JOB Finished, dumping output")
        self.pbsLog = (Popen("echo /localgrid/"+userName+"/"+self.taskName+".o*", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()).strip()

        if os.path.exists(self.pbsLog):
            for line in open(self.pbsLog):
                self.log.output(line.strip())
        else:
            self.log.output("-> No output file found???")

    def process(self):

        global filesToMerge, plotsToMerge

        if not options.local:
            #if self.monsterFile == "":
            self.createXMLCFG()

            self.setupWorkingDir()
            self.createPBSCFG()
#            sys.exit(1)
            self.submitPBSJob()
            self.checkPBSJob()

            # if self.monsterFile == "":
            self.log.output("self.inputFileNr = "+str(self.inputFileNr))
            if not self.inputFileNr == -1:
                outRootFile = Popen("ls "+self.workingDir+"/LightTree/*.root", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().split("/LightTree/")[1].split("\n")[0]
                outPlotsFile = Popen("ls "+self.workingDir+"/PlotsMacro/*.root", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().split("/PlotsMacro/")[1].split("\n")[0]
                if self.inputFileNr == 1:
                    filesToMerge.append(outRootFile.split(".root")[0])
                    plotsToMerge.append(outPlotsFile.split(".root")[0])
                newOutRootFile = outRootFile.split(".root")[0] + "_" + str(self.inputFileNr) + ".root"
                newOutPlotsFile = outPlotsFile.split(".root")[0] + "_" + str(self.inputFileNr) + ".root"
                self.log.output(Popen("cp -f "+self.workingDir+"/LightTree/"+outRootFile+" "+self.resultsDir+newOutRootFile, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
                self.log.output(Popen("cp -f "+self.workingDir+"/PlotsMacro/"+outPlotsFile+" "+self.plotsMacroDir+newOutPlotsFile, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
# self.log.output(Popen("cp -f "+self.workingDir+"/LightTree/"+outRootFile+" "+self.resultsDir+newOutRootFile+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
            else:
                self.log.output(Popen("cp -f "+self.workingDir+"/LightTree/*.root "+self.resultsDir+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
                self.log.output(Popen("cp -f "+self.workingDir+"/PlotsMacro/*.root "+self.plotsMacroDir+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
            # else:
            #    self.log.output(Popen("cp -f "+self.workingDir+"/LikelihoodResults_ASCII/*.txt* "+self.resultsDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
#           #     self.log.output(Popen("cp -f "+self.workingDir+"/LikelihoodResults_ASCII/*.txt* "+self.resultsDir+"; rm -rf "+self.workingDir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

            # clean up
            os.remove(self.pbsFile)
            os.remove(self.pbsLog)

        else:
            print "NOT YET IMPLEMENTED LOCAL RUNNING"

        return True

##############
## WorkFlow ##
##############


class WorkFlow (threading.Thread ):

    def __init__(self, nThread, *args, **kwds):

        global options
        global log
        global timestamp
        global nJobs

        self.nThread = nThread
        threading.Thread.__init__(self, *args, **kwds)
        self.keepAlive = bool(True)

        if int(nJobs) == 1:
            self.log = log
        else:
            self.log = logHandler("logs/"+options.TaskName+"_"+timestamp+"_thread"+str(self.nThread)+".txt")

    def stop (self):
        self.keepAlive = bool(False)

    def run (self):
        global nJobs
        # our thread runs forever
        while not jobsPool.empty():
            #print self.nThread
            job = jobsPool.get()
            job.setlog(self.log)
            if int(nJobs) == 1:
                log.output("-> Thread "+str(self.nThread)+" Processing Job: "+str(job.nJob+1)+"/"+str(nJobs))
            else:
                log.output("-> Thread "+str(self.nThread)+" Processing Job: "+str(job.nJob+1)+"/"+str(nJobs)+" (LogFile: "+self.log.logFile+")")

            if (job.process()):
                log.output("-> Thread "+str(self.nThread)+" Finished Job: "+str(job.nJob+1)+"/"+str(nJobs))

###############
### OPTIONS ###
###############

optParser = OptionParser()
optParser.add_option("-t","--taskName", dest="TaskName",default="AnomCoup",
                     help="TaskName that will be used as prefix for workingdir and logs", metavar="")
optParser.add_option("-d","--workingDir", dest="WorkingDir",default="",
                     help="Directory containing the needed TopTreeAnalysis setup and the file containing input datasets (inputSamples.txt)", metavar="")
optParser.add_option("-m","--mail", dest="Mail", default="",
                     help="E-mail adress to inform when the script finished", metavar="")
optParser.add_option("-o","--log-stdout", action="store_true", dest="stdout",default=bool(False),
                     help="Write the main log file to the stdout", metavar="")
optParser.add_option("-l","--run-local", action="store_true", dest="local",default=bool(False),
                     help="Use local CPUs", metavar="")
# optParser.add_option("-M","--Monsters", dest="monsterDir", default="",
#                     help="Directory containing the monsters (if this is set, the Likelihood ASCII files will be calculated)", metavar="")
(options, args) = optParser.parse_args()

if options.WorkingDir == "":
    print "No working dir provided! For help use python PBS_KinFit.py -h"
    sys.exit(1)

# -------- #
#   MAIN   #
# -------- #

# -- SETTINGS

RootInstallation = "/user/aolbrech/LocalSoftWare/root_v5.34.05_James/"  # /Software/LocalSoft/root_5.32_stijn/"

SleepTime = int(60)  # time to sleep between checking of job status
# SleepTime = int(10)  # time to sleep between checking of job status

# -- END SETTINGS

# special dirs
if not os.path.exists("logs"):
    os.mkdir("logs")
if not os.path.exists("Results"):
    os.mkdir("Results")

# timestamp
timestamp = strftime("%d%m%Y_%H%M%S")  # need a timestamp for dirs and logfiles

# logging
if not options.stdout:
    # log = logHandler("logs/log_"+timestamp+".txt")
    log = logHandler("logs/"+options.TaskName+"_"+timestamp+".txt")
else:
    log = logHandler("")

# start the loop
log.output("*** TopTreeAnalysis KinFit Batch Job system ***")
log.output("Making a dump of the configuration used:")
log.output("  --> TaskName:\t\t"+options.TaskName)
log.output("  --> WorkingDir:\t\t"+options.WorkingDir)
log.output("  --> log-stdout:\t\t"+str(options.stdout))
log.output("  --> run-local:\t\t"+str(options.local))
log.output("  --> Mail:\t\t"+options.Mail)
# log.output("  --> Directory containing the monsters (to produce Likelihood ASCII files from):\t\t"+options.monsterDir)

# get username
userName = Popen('echo $USER', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().strip()

# copy root to localgrid
if not options.local:
    if not os.path.exists("/localgrid/"+userName+"/root"):
        log.output("-> Copying root to localgrid")
        cmd = "cp -vfr "+RootInstallation+" /localgrid/"+userName+"/root"
        Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        ver = open("/localgrid/"+userName+"/root/ver", "w")
        ver.write(RootInstallation)
        ver.close()
    else:
        rootver = ""
        for ver in open("/localgrid/"+userName+"/root/ver", "r"):
            rootver = ver
            if not rootver == RootInstallation:
                log.output("-> RE-Copying root to localgrid") 

                cmd = "mv /localgrid/"+userName+"/root /localgrid/"+userName+"/root_old_"+timestamp
                cmd += " ;cp -vfr "+RootInstallation+" /localgrid/"+userName+"/root"
                
                Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                ver = open("/localgrid/"+userName+"/root/ver", "w")
                ver.write(RootInstallation)
                ver.close()
        
# sys.exit(1);

# Read in the inputSamples.txt file and create our Queue to store jobs:
jobsPool = Queue.Queue(0)
nJobs = int(0)
filesToMerge = []
plotsToMerge = []

# if options.monsterDir == "":
log.output("Creating the light trees...")
for line in open(options.WorkingDir+"/inputSamples.txt"):

    if len(line.split("#")) < 2:
        splitted = line.split(":")
        if len(splitted) > 1:
            log.output("Process: " + splitted[0] + "  " + splitted[1])
            pnfsPath = splitted[4].split("\n")[0]
            if (len(splitted[0].split("TT")) > 1 and not len(splitted[0].split("Hadronic")) > 1) or (len(splitted[1].split("InvertedIso")) > 1 and len(splitted[0].split("Data")) > 1) : # or len(splitted[0].split("Data_SingleElectron_Run2012D_Prompt_1")) > 1 :
                # process only one inputFile per job
                # cmd = "ls -l "+pnfsPath+" | grep TopTree_Skimmed_1.root | wc -l"        # --> Run only 1 file!
                cmd = "ls -l "+pnfsPath+" | grep TopTree_Skimmed_ | grep root | wc -l"
                nInFiles = int(Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
                print "nInFiles which will be considered is : ", nInFiles
                
                for i in range(1,nInFiles+1):
                    log.output("  --> with inputFile: "+str(i))
                    job = KinFitHandler(nJobs, splitted[0], splitted[1], splitted[2], splitted[3], pnfsPath, i)  # , options.monsterDir)
                    jobsPool.put(job)
                    nJobs += 1
            else:
                job = KinFitHandler(nJobs, splitted[0], splitted[1], splitted[2], splitted[3], pnfsPath, -1)   #, options.monsterDir)
                jobsPool.put(job)
                nJobs += 1
# else:
#    log.output("Calculating Ideogram likelihoods...")
#    inFiles = Popen("ls "+options.monsterDir+"*.root", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().split("\n")
#
#    for file in inFiles:
#        if not file == "":
#            log.output("Proceesing input monster file:  |"+file+"|")
#            job = KinFitHandler(nJobs, "", "", "", "", "", -1, file)
#            jobsPool.put(job)
#            nJobs += 1
        
# sys.exit(1)

# start our threads
workers = []
for x in xrange ( int(nJobs) ):
   workers.append(WorkFlow(x))
   workers[x].start()
   time.sleep(20)
   # if not options.monsterDir == "":
   #    time.sleep(20)

# check if there is a worker that are not done
notDone=bool(True)

while notDone:
    notDone = False
    for worker in workers:
        if worker.isAlive():  # If there is one worker alive, we are still not finished
            notDone = bool(True)
    if not notDone:
        log.output("-> All jobs are DONE")
    time.sleep(10)
		
# Merge (hadd) all the output files from the same sample
# if options.monsterDir == "":
resultsDir = "./Results/RESULTS_"+options.TaskName+"_"+timestamp+"/"
for file in filesToMerge:
    log.output("Merging files: "+file+"_*.root")
    nInFiles = int(Popen("ls -l "+resultsDir+file+"_*.root | wc -l", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
    if nInFiles > 1:
        command = "export ROOTSYS="+RootInstallation+"; export PATH=$ROOTSYS/bin:$PATH; export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH; cd "+resultsDir+"; hadd "+file+".root"
        for i in range(1,nInFiles+1):
            command += " "+file+"_"+str(i)+".root"
        command += "; mkdir backup; mv "+file+"_*.root ./backup/"
#           command += "; rm -rfv "+file+"_*.root"
        log.output("Executing command : |"+command+"|")
        log.output(Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

plotsMacroDir = "./Results/PLOTSMACRO_"+options.TaskName+"_"+timestamp+"/"
for plot in plotsToMerge:
    log.output("Merging files:  "+plot+"_*.root")
    nInPlots = int(Popen("ls -l "+plotsMacroDir+plot+"_*.root | wc -l", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
    if nInPlots > 1:
        command = "export ROOTSYS="+RootInstallation+"; export PATH=$ROOTSYS/bin:$PATH; export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH; cd "+plotsMacroDir+"; hadd "+plot+".root"
        for i in range(1,nInPlots+1):
            command += " "+plot+"_"+str(i)+".root"
        command += "; mkdir backup; mv "+plot+"_*.root ./backup/"
        log.output("Executing command : |"+command+"|")
        log.output(Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
    
# Also merge the data stuff (nominal and invIso)
# for dataSet in ["Mu","Electron"]:
#    for isoType in ["InvertedIso","Nominal"]:
#        inFiles = Popen("cd "+resultsDir+"; ls *Data*"+dataSet+"*"+isoType+"*.root", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().split("\n")
#        nInFiles = len(inFiles)
#        if nInFiles > 1:
#            log.output("Merging "+dataSet+" "+isoType+" files")
#            command = "export ROOTSYS="+RootInstallation+"; export PATH=$ROOTSYS/bin:$PATH; export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH; cd "+resultsDir+"; hadd KinFit_LightTree_TopMassDiff_Data_"+dataSet+"_"+isoType+"_MERGED.root"
#            for i in range(0,nInFiles):
#                command += " "+inFiles[i]
#            command += "; mkdir backup"
#                command += "; rm -rfv "+resultsDir+"*Data*"+dataSet+"*"+isoType+"*.root"
#            log.output("Executing command : |"+command+"|")
#            log.output(Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
#            for i in range(0,nInFiles):
#                log.output(Popen("cd "+resultsDir+"; mv "+inFiles[i]+" ./backup/", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
    
# remove .xml files
log.output("Removing *.xml files from directory:  "+options.WorkingDir)
log.output(Popen("rm "+options.WorkingDir+"/*.xml", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())

# Clean up /localgrid/aolbrech directory !
log.output("Cleaning up /localgrid/aolbrech/ directory by moving all directories to /localgrid/"+userName+"/LightTree_PBSScript/"+options.TaskName+"_"+timestamp)
log.output(Popen("mv /localgrid/"+userName+"/"+options.TaskName+"_"+timestamp+"_job_* /localgrid/"+userName+"/LightTree_PBSScript/"+options.TaskName+"_"+timestamp+"/" , shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read())
# --> Maybe possible to add a hadd of the AnomCouplings.root file which is stored in each of these directories!!

# send mail when finished
if not options.Mail == "":
    log.output("Sending announcement to:  "+options.Mail)
    mailHandler = MailHandler(options.Mail)
    message = "The script finished ;-)\n\n"+Popen("cat "+options.WorkingDir+"/inputSamples.txt", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    mailHandler.sendMail("[PBS_KinFit.py], script finished", message)