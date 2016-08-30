"""
script to merge multiple json file into a single dictionnary writen in an other json file.
The script works as follows:
1) make a list of all the .json files in the current directory
2) create a list of dict where each element corresponds to the content of one json file
3) merge everything in a single dict and save it to a new json file



qpython 30.08.2016
"""

import json
import os  

# json file list
jsonFiles=[]
for fn in os.listdir('.'):
     if os.path.isfile(fn) and ".json" in fn:
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
merged_dict={}

# merge all the elements in the merged_dict
for i in dic_list:
    merged_dict.update(i)
 
   
print "the merged dict is \n"
print merged_dict


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
