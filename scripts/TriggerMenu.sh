#!/bin/bash

if [ -e mytest.txt ]; then
	rm -f mytest.txt
fi

if  ! [ -f LumiRunSeeker ]
then
	echo "Will now compile LumiRunSeeker.cc"
	g++ -I `root-config --incdir` `root-config --libs` LumiRunSeeker.cc -o LumiRunSeeker
fi

if  ! [ -x  LumiRunSeeker ]
then
	echo "Please make sure that LumiRunSeeker is executable :"
	echo "chmod u+x LumiRunSeeker"
	exit 1
fi

#DATASET[0]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/B1B1ToT1LNuT1LNu_UDD_185_260_250_Tune4C_8TeV-madgraph_v0_Summer12/jdhondt-B1B1ToT1LNuT1LNu_UDD_185_260_250_Tune4C_8TeV-madgraph_v0_Summer12-5e21f1ce3d56eca5bfd8d6eb27e183c7/18112013_140722/TOPTREE_jdhondt-B1B1ToT1LNuT1LNu_UDD_185_260_250_Tune4C_8TeV-madgraph_v0_Summer12-5e21f1ce3d56eca5bfd8d6eb27e183c7_18112013_140722/res/"
DATASET[0]="/home/dhondt/ProductionReleases/V5_0_5/CMSSW_5_3_12_patch2/src/ConfigurationFiles/MuEG/Run2012A-22Jan2013-v1/30012014_232925/TOPTREE_Run2012A-22Jan2013-v1_30012014_232925/res/"
#DATASET[1]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/DoubleMuParked/Run2012B-22Jan2013-v1/23062013_112422/TOPTREE_Run2012B-22Jan2013-v1_23062013_112422/res/"
#DATASET[2]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/DoubleMuParked/Run2012C-22Jan2013-v1/23062013_145427/TOPTREE_Run2012C-22Jan2013-v1_23062013_145427/res/"
#DATASET[3]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/DoubleMuParked/Run2012D-22Jan2013-v1/21062013_125249/TOPTREE_Run2012D-22Jan2013-v1_21062013_125249/res/"
#DATASET[4]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/MuEG/Run2012A-22Jan2013-v1/26062013_092536/TOPTREE_Run2012A-22Jan2013-v1_26062013_092536/res/"
#DATASET[5]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/MuEG/Run2012B-22Jan2013-v1/26062013_095541/TOPTREE_Run2012B-22Jan2013-v1_26062013_095541/res/"
#DATASET[6]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/MuEG/Run2012C-22Jan2013-v1/26062013_095546/TOPTREE_Run2012C-22Jan2013-v1_26062013_095546/res/"
#DATASET[7]="/home/dhondt/AutoMaticTopTreeProducer/CMSSW_5_3_7_patch6_TopTreeProd_53X_v5/src/ConfigurationFiles/MuEG/Run2012D-22Jan2013-v1/26062013_095551/TOPTREE_Run2012D-22Jan2013-v1_26062013_095551/res/"
for DIR in "${DATASET[@]}"
do
	echo "--------------------------------------------------------------------"
	echo -n "Dataset : "
	echo ${DIR} | awk -F/ '{ print $8"/"$9}'
	echo "--------------------------------------------------------------------"
	for FILE in `ls ${DIR} | grep "stdout" | sort -n -k2 -t_`
	do
		grep "HLTAnalyzer-Summary-RUN" "${DIR}${FILE}" >> mytest.txt
	done

	./LumiRunSeeker mytest.txt "$1"

	if [ -e mytest.txt ]; then
		rm -f mytest.txt
	fi
done
