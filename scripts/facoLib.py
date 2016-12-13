"""
This file contains various functions that are frequently used.
They are sorted by alphabetic order.
In oder to used them put "import <path>/facoLib as fl"
in the desired python file.
Then you can use fl.<anyfucntion> to get the desired function.
All these functions are tested in testFacoLib.py

qpython 31.08.2016
"""


import json
from math import sqrt
import os
import ROOT as rt
from tabulate import tabulate


# usefull variables for writing a tex file
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"


def combinedUncertainty(array, debug=False):
    # type: (list[float], bool) -> float
    """ Function that calculates the combined uncertainty of uncorelated
    uncertainties given in an array.
    """
    SumSquare = 0
    for element in array:
        SumSquare = SumSquare + element**2
        if debug:
            print element
            print element**2
            print SumSquare

    return sqrt(SumSquare)


def combinedRelUncertainty(central, array, debug=False):
    # type: (float, list[float], bool) -> float
    """Function that calculates the combined uncertainty (relative to a central
    value) of uncorelated uncertainties given in an array.
    """
    SumSquare = 0
    for element in array:
        SumSquare = SumSquare + element**2
        if debug:
            print element
            print element**2
            print SumSquare

    return sqrt(SumSquare / central**2)


def makeTable(fileName, doubleArray, header, savePDF=False,
              caption="nice caption bro!!", debug=False):
    # type: (str, list[list[float]], str, bool,
    #       str, bool) -> None
    """Function that convert a array of array into a latex table with the
    possibility to compile it it directly to get the pdf output
    """

    # create new repository
    if not os.path.exists("tables"):
        os.makedirs("tables")

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


def getDictFromJson(matchingPattern, vetoPattern="", debug=False):
    # type: (str, str, bool) -> dict
    """Function that returns a dictionary gotten from the
    merging of multiple json files.
    """

    # json file list
    jsonFiles = []
    pathToJsonFile = '/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/scripts/jsonFiles/' 
    for fn in os.listdir(pathToJsonFile):
        if debug:
            print fn
        if ".json" in fn and matchingPattern in fn:
            if vetoPattern == "":
                jsonFiles.append(fn)
            elif vetoPattern not in fn:
                jsonFiles.append(fn)

    if debug :
        print "the list of json files to be use is \n", jsonFiles

    # list of dictionaries. (one element per json file/ dictrionary)
    dic_list = []

    # loop over json file list previously created
    for jsonFile in jsonFiles:
        with open(pathToJsonFile+jsonFile, 'r') as f:
            try:
                dic_list.append(json.load(f))
            except ValueError:
                print "error!!!"

    # declare new dictionary
    my_dict = {}

    # merge all the elements in the merged_dict
    for i in dic_list:
        my_dict.update(i)

    # return the dict
    return my_dict


def floatToPercent(value):
    # type (float) -> str
    """ Converts a value into a string in percentage."""
    if value > 1:
        raise ValueError("Value is bigger than 1!")
    else:
        return str(100.0 * value) + " %"


def makeEffienciency(hist1, hist2, ymin=False, ymax=False, norm=False):
    # type: (rt.TH1D, rt.TH1D, bool, bool, bool) -> None
    """ Makes a good looking efficiency plots with
    the possibility to show the ratio plots.
    """
    if norm:
        print 'scaling!'
        try:
            print 'scale 1: ', 1 / hist1.Integral()
            print 'scale 2: ', 1 / hist2.Integral()
            hist1.Scale(1 / hist1.Integral())
            hist2.Scale(1 / hist2.Integral())
        except(ZeroDivisionError):
            pass
    retH = hist1.Clone()
    try:
        retH.Divide(hist2)
    except(TypeError):
        # this is the error you get if hist2 is a stack
        hList = hist2.GetHists()
        sumHist = hist1.Clone("sumHist")
        sumHist.Reset()
        for h in hList:
            sumHist.Add(h)
        retH.Divide(sumHist)
    except(AttributeError):
        # this is the error you get if hist1 is a stack
        print "Did you use a stack as argument 1? please use stack as argument 2!"
        raise AttributeError
    if ymax or ymin:
        retH.GetYaxis().SetRangeUser(0.5, 1.5)
        retH.SetLineColor(hist2.GetLineColor())

    retH.GetAxisX().SetLabelSize(hist1.GetHistogram().GetAxisX().GetLabelSize())
    retH.GetAxisX().SetLabelOffset(hist1.GetHistogram().GetAxisX().GetLabelOffset())
    retH.GetAxisX().SetTitleSize(hist1.GetHistogram().GetAxisX().GetTitleSize())
    retH.GetAxisX().GetTitleOffset(hist1.GetHistogram().GetAxisX().GetTitleOffset())
    ROOT.SetOwnership(retH, 0)
    return retH
