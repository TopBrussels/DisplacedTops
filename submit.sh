#!/bin/bash

# line starting with # P B S will be interpreted


# setting numbers of nodes and cores per nodes
#PBS -l nodes=1:ppn=4

# setting wall time (<HH:MM:SS>)
#PBS -l walltime=05:00:00

# settting the qeue
#PBS -q localgrid


# Setting up
cd /localgrid/qpython/TopBrussels7X/CMSSW_7_4_14/src
eval `scramv1 runtime -sh` #->cmsenv
cd /localgrid/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops

# running the macro
python /localgrid/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops/LaunchParallelMacro.py

# print hostname
echo $HOSTNAME

# copy files over /user/
# to be added..