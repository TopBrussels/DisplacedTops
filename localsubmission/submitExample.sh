#!/bin/bash

###
# This a basic example script that will submit jobs using the Belgian T2 submission queue
# note that here after every "# PBS" will be interpreted even though they look like they are comments!

# you can find more explanation on the official twiki
# http://mon.iihe.ac.be/trac/t2b/wiki/localSubmission

# you can view some information about your jobs here:
# http://mon.iihe.ac.be/jobview/overview.html

# usage: qsub submit.sh
###



# choose the queue to be used
#PBS -q localgrid

# fix the wall time. Make sure that this value is slightly large than the time your job would need
#PBS -l walltime=02:00:00

# setting up your code and your env
cd /user/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops
cmsenv

# want you really want to do!!
python /user/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops/LaunchParallelMacro.py

# to be added...
# merge your file and copy them under the desired location