#!/bin/bash 

iterator=0

if [[ -n $1 ]] #check if variable is not empty
then
    if [[ $1 == "test" ]]
    then
	echo "Using test option! Only submiting one file per sample"
	cd test
	for f in ./submit*.sh
	do
	    qsub $f
	    let "iterator++"
	done
	cd -
    else
	echo "A filter is being used! The filter is" $1
	cd output
        for f in ../submit_$1*.sh
        do
            qsub $f
	    let "iterator++"
        done
        cd -
    fi

else
    cd output
    for f in ../submit*.sh
    do
	qsub $f
	let "iterator++"
    done
    cd -

fi

echo "the total number of submited job is " $iterator
