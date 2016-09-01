"""
This file contains various functions that are frequently used.
They are sorted by alphabetic order.
in oder to used them put "import <path>/facoLib as fl" in the desired python file.
Then you can use fl.<anyfucntion> to get the desired function.
All these functions are tested in testFacoLib.py 

qpython 31.08.2016 
"""


import json
from math import sqrt
import os
from tabulate import tabulate


# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"


# function that calculate the combined uncertainty of uncorelated
# uncertainties given in an array
def combinedUncertainty(array, debug=False):
    SumSquare = 0
    for element in array:
        SumSquare = SumSquare + element**2
        if debug:
            print element
            print element**2
            print SumSquare

    return sqrt(SumSquare)


# function that calculate the combined uncertainty (relative to a central
# value) of uncorelated uncertainties given in an array
def combinedRelUncertainty(central, array, debug=False):
    SumSquare = 0
    for element in array:
        SumSquare = SumSquare + element**2
        if debug:
            print element
            print element**2
            print SumSquare

    return sqrt(SumSquare / central**2)


# function that convert a double array into a latex table with possibility
# to save it as pdf
def makeTable(fileName, doubleArray, header, savePDF=False, caption="nice caption bro!!", debug=False):

    # open the output file
    outputfile = "tables/" + fileName + ".tex"
    fout = open(outputfile, "w")

    # begin document and table
    fout.write("\\documentclass{article}" +
               newLine + "\\begin{document}" + newLine)
    fout.write("\\renewcommand{\\arraystretch}{1.2}" + newLine)
    fout.write("\\begin{table}" + newLine)
    fout.write("\\caption{ " + caption + "}" + newLine)

    # the actual tabular
    fout.write(tabulate(doubleArray, header, tablefmt="latex"))

    # end table and doc
    fout.write("\\end{table}" + newLine)
    fout.write("\\end{document}" + newLine)
    fout.close()

    if savePDF:
        # compile tex into pdf
        cmd = "pdflatex " + outputfile
        os.system(cmd)

        # mv table
        cmd = "mv " + fileName + ".pdf" + " tables/"
        os.system(cmd)

        # clean the mess
        cmd = "rm " + fileName + ".aux"
        os.system(cmd)
        cmd = "rm " + fileName + ".log"
        os.system(cmd)


# function that returns a dictionary gotten from the merging of multiple
# json files
def getDictFromJson(matchingPattern, vetoPattern="", debug=False):

    # json file list
    os.chdir("jsonFiles")
    jsonFiles = []
    for fn in os.listdir('.'):
        print fn
        if ".json" in fn and matchingPattern in fn:
            if vetoPattern == "":
                jsonFiles.append(fn)
            elif vetoPattern not in fn:
                jsonFiles.append(fn)

    print "the list of json files to be use is \n", jsonFiles

    # list of dictionaries. (one element per json file/ dictrionary)
    dic_list = []

    # loop over json file list previously created
    for jsonFile in jsonFiles:
        with open(jsonFile, 'r') as f:
            try:
                dic_list.append(json.load(f))
            except ValueError:
                print "error!!!"

    # declare new dictionary
    Yield_dict = {}

    # merge all the elements in the merged_dict
    for i in dic_list:
        Yield_dict.update(i)

    # return the dict
    return Yield_dict
