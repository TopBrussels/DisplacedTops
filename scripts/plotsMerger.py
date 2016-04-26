"""
script that creates a full presenation from a list of plots
it takes one agrument and it is the path on which to search for the list of plots
qpython April 2016
"""

from datetime import datetime
import os, sys


# path to the list of file
pathToListOfFile = sys.argv[1]

newLine = " \n"

# Define time variables
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
date = dd+"_"+mm+"_"+yyyy


# make a list of files of the current dir
files = os.listdir(pathToListOfFile)

#begining of the tex file
outputFileTrunc = pathToListOfFile+"merged"+date
outputFile = outputFileTrunc+".tex"
fout = open (outputFile, "w")
fout.write("\\documentclass{beamer}")
fout.write(newLine)
fout.write("\\title{Merged Presentation}")
fout.write(newLine)
fout.write("\\author{Quentin}")
fout.write(newLine)
fout.write("\\date{\\today}")
fout.write(newLine)
fout.write("\\begin{document}")
fout.write(newLine)
fout.write("\\frame{\\titlepage}")
fout.write(newLine)
fout.write("\\section[Outline]{}")
fout.write(newLine)
fout.write("\\frame{\\tableofcontents}")
fout.write(newLine)
fout.write("")





# loop over all the plots
for file in files:
    # apply a certain filter (to be tuned)
    if "StackLogY.png" in file and "nvtx" not in file:
        print(file)
        
        # for each file create a new frame and put the figure
        fout.write(newLine)
        fout.write(newLine)
        fout.write("\\frame{")
        fout.write(newLine)
        fout.write("\\frametitle{}")
        fout.write(newLine)
        fout.write("\\begin{figure}[htbp]")
        fout.write(newLine)
        fout.write("\\begin{center}")
        fout.write(newLine)
        fout.write("\\includegraphics[width=0.90\\textwidth,height=0.58\\textwidth ]{"+file+"}")
        fout.write(newLine)
#        fout.write("\\caption{default}")
        fout.write(newLine)
        fout.write("\\label{default}")
        fout.write(newLine)
        fout.write("\\end{center}")
        fout.write(newLine)
        fout.write("\\end{figure}")
        fout.write(newLine)
        fout.write("}")
        fout.write(newLine)
        fout.write(newLine)



# write the end of the document and close
fout.write("\end{document}")
fout.close()


# mv the tex file to the desired directory
#cmd="mv "+outputFile+" "+pathToListOfFile
#os.system(cmd)

#mv to the desired directory
#cmd="cd "+pathToListOfFile
#os.system(cmd)


# compile
#cmd="pdflatex "+outputFile


# clean add waiting time
#cmd="rm "+outputFileTrunc+".aux"
#os.system(cmd)
#cmd="rm "+outputFileTrunc+".log"
#os.system(cmd)
