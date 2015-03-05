#simple Ntuple maker

#To set up the code follow the following recipe

# Firstly, set up CMSSW
cmsrel CMSSW_7_2_1_patch1
cd CMSSW_7_2_1_patch1/src
cmsenv

# Get TopTreeProducer from git
git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer
cd TopBrussels/TopTreeProducer/
git checkout CMSSW_70X
cd src 
make
cd ../../..

# Get TopTreeAnalysisBase from git
git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/
git checkout master
make
cd ../../

# Get DisplacedTop directory from git
cd TopBrussels
git clone git@github.com:TopBrussels/DisplacedTops.git DisplacedTops
cd DisplacedTops
git fetch origin
git checkout CMSSW_7_2_X


# Compile
source compile.sh

# Run
./Ntupler testconfig.xml

# make plots
cd scripts
python maked0plot.py






