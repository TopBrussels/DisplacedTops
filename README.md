simple Ntuple maker

To set up the code follow the following recipe


## Firstly, set up CMSSW
export SCRAM_ARCH=sslc6_amd64_gcc493

cmsrel CMSSW_7_6_3

cd CMSSW_7_6_3/src

cmsenv

git cms-init

## Get TopTreeProducer from git
git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_76X

cd src

make

cd ../../..

## Get TopTreeAnalysisBase from git
git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

git checkout CMSSW_76X

make

cd ../../

## Get DisplacedTop directory from git

git clone git@github.com:TopBrussels/DisplacedTops.git TopBrussels/DisplacedTops

cd TopBrussels/DisplacedTops

git checkout master

## Compile and create executable from local macro

cd ../../

scram b -j8

cd -

# to compile all
source compile.sh

# to compile a single Macro
source compile.sh <Macro.cc>

# Make ntuple
./TreeMaker WWToLNuQQ Diboson 1 390 1 2 1 41297.4650368 47.693 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/WWToLNuQQ_13TeV-powheg/crab_WWToLNuQQ_13TeV-powheg-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151029_124630/0000/TOPTREE_7.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/WWToLNuQQ_13TeV-powheg/crab_WWToLNuQQ_13TeV-powheg-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151029_124630/0000/TOPTREE_6.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/WWToLNuQQ_13TeV-powheg/crab_WWToLNuQQ_13TeV-powheg-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151029_124630/0000/TOPTREE_5.root   ElEl   3  0  2000000

# make MSPlots
./TreeProcesser ElEl


# make plots --obselete--
cd scripts

python maked0plot.py
