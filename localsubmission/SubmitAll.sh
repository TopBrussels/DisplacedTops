#!/bin/bash 
if [[ -n $1 ]] #check if variable is not empty
then
    if [[ $1 == "test" ]]
    then
	cd test
	for f in ./submit*.sh
	do
	    qsub $f
	done
	cd -
    fi

else
    cd output
    for f in ../submit*.sh
    do
	qsub $f
    done
    cd -

fi
