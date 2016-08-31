"""
script to merge multiple json file into a single dictionnary writen in an other json file.
The script works as follows:
1) make a list of all the .json files in the current directory
2) create a list of dict where each element corresponds to the content of one json file
3) merge everything in a single dict and save it to a new json file
4) use the merged dictionary to produce the systematics uncertainty table with the total error.

qpython 30.08.2016
"""

import json
from math import sqrt
import os  
from tabulate import tabulate

# json file list
os.chdir("jsonFiles")
jsonFiles=[]
for fn in os.listdir('.'):
     print fn
     if ".json" in fn:
          jsonFiles.append(fn)

print "the list of json files to be use is \n", jsonFiles


# list of dictionaries. (one element per json file/ dictrionary)
dic_list=[]


# loop over json file list previously created
for jsonFile in jsonFiles:
    with open(jsonFile, 'r') as f:
        try:
            dic_list.append(json.load(f))
        except ValueError:
            print "error!!!"
    
# declare new dictionary
Yield_dict={}

# merge all the elements in the merged_dict
for i in dic_list:
    Yield_dict.update(i)
 
   
print "the merged dict is \n"
print Yield_dict


dict_cstSystType={"lumi":0.05, "trigger_sf_electron":0.20, "trigger_sf_muon":0.10, "trk_sf_electron":0.24, "trk_sf_muon":0.12 }

# 
#MuonSystTypes=["Cross-Section","evt_Pusf","Sf_Iso_Muon","sf_id_muon","lumi","trigger_sf_muon", "trk_sf_muon"]
MuonSystTypes=["evt_puSF","sf_iso_muon","sf_id_muon","lumi","trigger_sf_muon", "trk_sf_muon"]
electronSystTypes=["Cross-section","evt_puSF","sf_reco_electron","sf_id_electron","trigger_sf_electron", "trk_sf_electron"]
sampleNames=["ZG","WWTo2l2Nu"]

# double array for table writting 
doubleArray = []
singleArray = []
headers=["Dataset"]


systTypes=list(MuonSystTypes)

# loop over the samples 
for sampleName in sampleNames:

        singleArray = [sampleName]

        SumUnc=0
        CurrentUnc=0
        # loop over the systematics types 
        for systType in systTypes :
            headers.append(systType)


            if systType in dict_cstSystType:
                 print "faco!!!", systType
                 uncertaintyMax = 0
            else :
                # one variable per yield
                YieldDown= Yield_dict[sampleName+systType+"down"]
                YieldUp= Yield_dict[sampleName+systType+"up"]
                YieldCentral =  Yield_dict[sampleName+"Centralup"]
                
                # two diff wrt to central
                diffDown=abs((YieldCentral-YieldDown)/YieldCentral)
                diffUp=abs((YieldCentral-YieldUp)/YieldCentral)
    
                # max of the two effs
                uncertaintyMax=  max(diffDown,diffUp)

            uncertaintyMaxPercent = str(100*uncertaintyMax)+" %"
    
            # Sum the square of the uncertainties
            CurrentUnc=uncertaintyMax
            SumUnc=SumUnc+CurrentUnc*CurrentUnc

            # 
            singleArray.append(uncertaintyMaxPercent)


        # eo systType loop 
        
        # add the tot uncertainty for each samples
        headers.append("total")
        singleArray.append(str(100*sqrt(SumUnc))+" %")
        

        doubleArray.append(singleArray)

# eo loop over the samples 







print tabulate(doubleArray, headers, tablefmt="latex")

# writing results in a tex file 
outputFile = "tables/SystematicsTable"+chan+".tex"
fout = open (outputFile, "w")
fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
fout.write("\\begin{table}"+newLine)
fout.write("\\caption{ " + "Systematic Uncertainty"+chan.replace("_"," ")+ "}"+newLine)

# the actual tabular
fout.write(tabulate(doubleArray, headers, tablefmt="latex"))


# end of table 
fout.write("\\end{table}"+newLine)
fout.write("\\end{document}"+newLine)
fout.close()





# save the dict in a file
with open('merged.json', 'w') as f:
    json.dump(merged_dict, f)



# simple test
"""
t = [{'ComSMS': 'true'}, {'ComMail': 'true'}, {'PName': 'riyaas'}, {'phone': '1', 'phaco':'awesome'}]

fat_big_dict={}

for i in t:
    fat_big_dict.update(i)

print fat_big_dict
"""
