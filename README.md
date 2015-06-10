simple Ntuple maker

To set up the code follow the following recipe


## Firstly, set up CMSSW
export SCRAM_ARCH=slc6_amd64_gcc491

cmsrel CMSSW_7_4_2

cd CMSSW_7_4_2/src

cmsenv

## Get TopTreeProducer from git
git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_74X

cd src 

make

cd ../../..

## Get TopTreeAnalysisBase from git
git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

git checkout CMSSW_74X

make

cd ../../

## Get DisplacedTop directory from git

git clone git@github.com:TopBrussels/DisplacedTops.git TopBrussels/DisplacedTops

cd TopBrussels/DisplacedTops

git checkout master

## Compile
source compile.sh

# Run
./Ntupler testconfig.xml

# make plots
cd scripts

python maked0plot.py
