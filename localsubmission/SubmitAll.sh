#!/bin/bash 

for f in SubmitScripts/*.sh
do
    qsub $f
done
