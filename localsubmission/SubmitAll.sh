#!/bin/bash 


for f in SubmitScripts/*.sh
#for f in SubmitScripts/submit_Data1.sh SubmitScripts/submit_Data2.sh
do
#    echo $f
#    qsub $f -q express
    qsub $f
#    qsub $f -q express -l walltime=00:05:00 
done
