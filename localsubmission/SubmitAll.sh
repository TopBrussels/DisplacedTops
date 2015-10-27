#!/bin/bash 

for f in SumbitScripts/*.sh
do
    qsub $f
done
