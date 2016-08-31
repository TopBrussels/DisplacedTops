"""
This file contains various functions that are frequently used.
They are sorted by alphabetic order.
in oder to used them put "import <path>/facoLib as fl" in the desired python file.
Then you can use fl.<anyfucntion> to get the desired function.
All these functions are tested in testFacoLib.py 

qpython 30.08.2016 
"""




from math import sqrt
import os
from tabulate import tabulate


# usefull variables for writing a tex file                                                                                          
hLine = "\\hline\n"
endLine = " \\\\ "
newLine = " \n"




# function that calculate the combined uncertainty of uncorelated uncertainties given in an array
def combinedUncertainty (array, debug=False):
    SumSquare=0
    for element in array:
        SumSquare=SumSquare+element**2
        if debug:
            print element
            print element**2
            print SumSquare
        
    return sqrt(SumSquare)


# function that calculate the combined uncertainty (relative to a central value) of uncorelated uncertainties given in an array 
def combinedRelUncertainty (central, array, debug=False):
    SumSquare=0
    for element in array:
        SumSquare=SumSquare+element**2
        if debug:
            print element
            print element**2
            print SumSquare
        
    return sqrt(SumSquare/central**2)


# function that convert a double array into a latex table with possibility to save it as pdf
def makeTable (fileName, doubleArray, header, savePDF=False, caption="nice caption bro!!", debug=False):

    # open the output file
    outputfile="tables/"+fileName+".tex"
    fout = open (outputfile, "w")

    # begin document and table
    fout.write("\\documentclass{article}"+newLine+"\\begin{document}"+newLine)
    fout.write ("\\renewcommand{\\arraystretch}{1.2}"+newLine)
    fout.write("\\begin{table}"+newLine)
    fout.write("\\caption{ "+caption+  "}"+newLine)

    # the actual tabular                                                                                                            
    fout.write(tabulate(doubleArray, header, tablefmt="latex"))
    
    # end table and doc
    fout.write("\\end{table}"+newLine)
    fout.write("\\end{document}"+newLine)
    fout.close()


    if savePDF:
        # compile tex into pdf                                                                                                          
        cmd="pdflatex "+outputfile
        os.system(cmd)

        # mv table                                                                                                                      
        cmd="mv "+fileName+".pdf"+" tables/"
        os.system(cmd)

        # clean the mess                                                                                                                
        cmd="rm "+fileName+".aux"
        os.system(cmd)
        cmd="rm "+fileName+".log"
        os.system(cmd)

