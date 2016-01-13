#!/bin/bash

# if there is one arg compile only this .cc
if [[ -n $1 ]]
then
    ccfile=$1

    ofile=`echo $ccfile |sed 's/\.cc$//g'`
    echo "compiling : " $ccfile ", executible name: " $ofile
    g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent74 -l TopTreeAna74 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
    cp ~/lib/libTopTreeAnaContent74.so /localgrid/qpython/lib/
    cp ~/lib/libTopTreeAna74.so /localgrid/qpython/lib/
    

# if there is no arg compile all .cc
else
    for ccfile in ./*.cc
    do
	ofile=`echo $ccfile |sed 's/\.cc$//g'`
	echo "compiling : " $ccfile ", executible name: " $ofile
	g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent74 -l TopTreeAna74 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile

    done
fi