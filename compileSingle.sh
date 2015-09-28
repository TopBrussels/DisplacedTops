#!/bin/bash

ccfile=$1

ofile=`echo $ccfile |sed 's/\.cc$//g'`
echo "compiling : " $ccfile ", executible name: " $ofile
g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent76 -l TopTreeAna76 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
