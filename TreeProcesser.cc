#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <sstream>
//#include "TDataType.h"
//#include "TTree.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"



using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;

// ---------------------
// bo functions prototype

// make MSPlots for each datasets 
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, bool plotAbs = false);

// draw the MSPlots
void MSPCreator ();

// make TH2F histo for multiple dataset
void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX, int nBinsY, float lowY, float highY, string sVarofinterestY, string xmlNom, string TreePath );

// convert int to string
std::string intToStr (int number);

// make cutflow table for each dataset
void CutFlowMaker (string xmlNom, string TreePath);

// eo functions prototype
// ---------------------

// List of bools that are set to fault as default
// These bools are overwritten depending on the arguments of the main
Bool_t debug = false;
Bool_t debug_plot = false;
bool DileptonElMu = false;
bool DileptonMuMu = false;
bool DileptonElEl = false;
bool bbMu = false;
bool bbEl = false;

// string prefix/sufix
string channelpostfix = "";
string regSuf = "";

//applying all appropriate scale factors for individual objects if the bool is set to true
Bool_t applyElectronSF = false; 
Bool_t applyMuonSF = false; 
Bool_t applyPUSF = false; 
Bool_t applyGlobalSF = false; 


// main fucntion
// get arguments from the cmd line and get appropriate Trees
// make a list of MSplots
// make a list of TH2F
// draw all the MSplots and the TH2F

int main(int argc, char* argv[])
{
  //  if (argc < 7 ) cout << "Not enough arguments! " << endl;
  if (debug){
    cout << "argc = " << argc << endl; 
    for(int i = 0; i < argc; i++)
      {
	cout << "argv[" << i << "] = " << argv[i] << endl; 
      }
  }


  //Placing arguments in properly typed variables

  const string channel = argv[1];
  debug = strtol(argv[2],NULL,10); 
  debug_plot = strtol(argv[3],NULL,10); 
  applyElectronSF = strtol(argv[4],NULL,10);
  applyMuonSF = strtol(argv[5],NULL,10);
  applyPUSF = strtol(argv[6],NULL,10);
  applyGlobalSF = strtol(argv[7],NULL,10);

  cout << "---------------------------------------------------" << endl;
  if (debug) cout << "debug is true --> adding some cout..." << endl;
  if (debug_plot) cout << "debug_plot is true --> using a reduced set of plots to go faster..." << endl;
  if (applyElectronSF) cout << "applying ElectronSF..." << endl;
  if (applyMuonSF) cout << "applying MuonSF..." << endl;
  if (applyPUSF) cout << "applying PUSF..." << endl;
  if (applyGlobalSF) cout << "applying GlobalSF..." << endl;
  cout << "---------------------------------------------------" << endl << endl;
			 
  
  

  string xmlFileName;

  if(channel=="ElMu" || channel== "MuEl")
    {
      cout << " --> Using the Muon-Electron channel..." << endl;
      channelpostfix = "_MuEl";
      xmlFileName = "config/TreeProc_ElMuV4.xml";
      DileptonElMu=true;
    }
  else if(channel=="MuMu")
    {
      cout << " --> Using the Muon-Muon channel..." << endl;
      channelpostfix = "_MuMu";
      xmlFileName = "config/TreeProc_MuMuV4.xml";
      DileptonMuMu=true;
    }
  else if(channel=="ElEl")
    {
      cout << " --> Using the Electron-Electron channel..." << endl;
      channelpostfix = "_ElEl";
      //      channelpostfix += "/unskimmed";
      
      xmlFileName = "config/TreeProc_ElElV4.xml";
      DileptonElEl=true;
    }
  else if(channel=="bbMu")
    {
      cout << " --> Using the bbar+Muon control region..." << endl;
      channelpostfix = "_bbMu";
      xmlFileName = "config/TreeProc_bbMuV4.xml";
      bbMu=true;
    }
  else if(channel=="bbEl")
    {
      cout << " --> Using the bbar+Electron control region..." << endl;
      channelpostfix = "_bbEl";
      xmlFileName = "config/TreeProc_bbElV4.xml";
      bbEl=true;
    }
  else
    {
      cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
      exit(1);
    }
  
  string CraneenPath;

  
  if (debug_plot){
    //    xmlFileName = "config/FullSamplesMuMuV0TreeProc.xml";


    // only few plots!
    if (DileptonMuMu) {
      //      xmlFileName = "config/TreeProc_MuMuV0.xml";
      //      DatasetPlotter(5, -0.5, 4.5, "nMuons", xmlFileName,CraneenPath);
      DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
      //            DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(40, 0, 800, "evt_met", xmlFileName,CraneenPath);
    }
    if (DileptonElEl){
      //      xmlFileName = "config/TreeProc_ElElV0.xml";
      DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
      
    }

  }
    
    
  else{
    
    //    xmlFileName = "config/FullSamplesMuMuV9TreeProc.xml";
     

    //    cout << "CraneenPath is " << CraneenPath << endl;
    cout << "xmlFileName is " << xmlFileName << endl;


    // 
    // calling datasetPlotter to create MSPplots
    



    // bo selecting the right plots and/or variables depending on the final state


    // E-MU FINAL STATE
    if (DileptonElMu){
          
      // event plots
      //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
      DatasetPlotter(70, -0.5, 49.5, "nvtx", xmlFileName,CraneenPath);
    

      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath);

 
      // muon plots
      //      DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath);


      // electron-muon plots
      //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);
    }


    // MU-MU FINAL STATE
    if (DileptonMuMu){
          
      // event plots
      //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
      DatasetPlotter(50, -0.5, 49.5, "nvtx", xmlFileName,CraneenPath);
      DatasetPlotter(40, 0, 800, "evt_met", xmlFileName,CraneenPath);
          
     
      /*
      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath);
      */
      
 
      //      /*
      // muon plots
      DatasetPlotter(6, -0.5, 5.5, "nMuons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_muon[nMuons]", xmlFileName,CraneenPath);
      //	DatasetPlotter(20, 10, 50, "pt_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, -0.015, 0.015, "d0_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, -10, 10, "vz_muon[nMuons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0, 0.5, "v0_muon[nMuons]", xmlFileName,CraneenPath);


      // muonPairs plots
      DatasetPlotter(100, 0.0, 1000, "invMass_mumu[nMuonPairs]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.4, "deltaVz[nMuonPairs]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.06, "deltaV0[nMuonPairs]", xmlFileName,CraneenPath);
      //      */

      // electron-muon plots
      //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);
    }

    // El-El FINAL STATE
    //    if (0){
    if (DileptonElEl){
          
      // event plots
      //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
      DatasetPlotter(50, -0.5, 49.5, "nvtx", xmlFileName,CraneenPath);

      //      /*
      DatasetPlotter(40, 0, 800, "evt_met", xmlFileName,CraneenPath);
              

      // electron plots
      DatasetPlotter(6, -0.5, 5.5, "nElectrons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
      //	DatasetPlotter(20, 10, 50, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, -0.015, 0.015, "d0_electron[nElectrons]", xmlFileName,CraneenPath);

      //      DatasetPlotter(100, -10, 10, "vz_electron[nElectrons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0, 0.5, "v0_electron[nElectrons]", xmlFileName,CraneenPath);

      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, -0.009, 0.009, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
      
      // Dielectron plots
      DatasetPlotter(100, 0.0, 1000, "invMass_elel[nElectronPairs]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.4, "deltaVz[nElectronPairs]", xmlFileName,CraneenPath);
      //      DatasetPlotter(100, 0.0, 0.06, "deltaV0[nElectronPairs]", xmlFileName,CraneenPath);
//      */
      
      
      /*
      // muon plots
      //    DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath);
      */

      // electron-muon plots
      //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);

    }
    if (bbMu){
      //      DatasetPlotter(50, -0.05, 0.05, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(25, 0.0, 0.05, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath, true);
      //      DatasetPlotter(50, 0.0, 1.6, "pfIso_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, 0.05, 1.55, "relIso_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, 0, 1, "CSV_jet[nJets]", xmlFileName,CraneenPath);
      DatasetPlotter(6, 0.922, 1, "CSV_bjet[nBjets]", xmlFileName,CraneenPath);
    }

    if (bbEl){
      //      DatasetPlotter(50, -0.05, 0.05, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(25, 0.0, 0.05, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath, true);
      //      DatasetPlotter(50, 0.0, 1.6, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, 0.05, 1.55, "relIso_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, 0, 1, "CSV_jet[nJets]", xmlFileName,CraneenPath);
      DatasetPlotter(6, 0.922, 1, "CSV_bjet[nBjets]", xmlFileName,CraneenPath);
    }

  }


  // Making cut flow table

  /*
  
  // vector of string
  vector<string> CutFlowPresel;
  CutFlowPresel.push_back("initial");
  CutFlowPresel.push_back("trigger");
  CutFlowPresel.push_back("faco");

  
  SelectionTable CutFlowPreselTable(CutFlowPresel, datasets);
  //    CutFlowPreselTable.SetLuminosity(Luminosity);
  CutFlowPreselTable.SetPrecision(1);


  */
  


  // eo selecting the right plots and/or variables depending on the final state  

  // making 2D histograms


  /*
  if (DileptonMuMu){
    TH2FPlotter(20, -10, 10, "v0_muon[nMuons]", 20, -0.015, 0.015, "d0BeamSpot_muon[nMuons]]",  xmlFileName,CraneenPath);
  }
  */

  /*
  if (DileptonElEl){
    TH2FPlotter(20, -10, 10, "v0_electron[nElectrons]", 20, -0.015, 0.015, "d0BeamSpot_electron[nElectrons]",  xmlFileName,CraneenPath);
    cout << "trying TH2F" << endl;
  }
  //  */

  // calling the function that writes all the MSPlots in a root file
  MSPCreator ();

}
// eo main

void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, bool plotAbs)
{
  cout<<""<<endl;
  cout<<"RUNNING NOMINAL DATASETS"<<endl;
  cout<<""<<endl;

  /*
  // vector with size equals to the number of time we call the DatasetPlotter function
  vector<string> vect_of_plot(4);
  vect_of_plot = {"nElectrons", "pt_electron[nElectrons]", "pt_muon[nMuons]", "d0_muon[nMuons]"};  
  */
  
  const char *xmlfile = xmlNom.c_str();
  cout << "used config file: " << xmlfile << endl;
  
  string pathPNG = "myOutput";
  pathPNG += "_MSPlots/"+channelpostfix+"/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  if (debug) cout << "will start loading from xml file ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  if (debug) cout << "finished loading from xml file ..." << endl;
  
  //***************************************************CREATING PLOTS****************************************************
  //  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  //  outfile->cd();
  string plotname = sVarofinterest;  
  // make for loop here!!!
  MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str()); 

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath , slumi;
  
  int nEntries;
  float ScaleFactor, NormFactor, Luminosity;
  float DataEqLumi= 1.;


  // Declare all the variables that allows to store inforamtion from the trees

  // a vector of double
  int n_object = 0;
  double v_varofInterest_double [20];

  // a simple int
  int varofInterest;

  // a simple double
  double varofInterest_double;

  
  vector<string> v;
  // to avoid modifying original string
  // first duplicate the original string and return a char pointer then free the memory
  
  char delim[] = " []";
  char * dup = strdup(sVarofinterest.c_str());
  char * token = strtok(dup, delim);
  while(token != NULL){
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);

  if (debug){
    //    cout << "n_object is " << n_object << endl;
    cout << "sVarofinterest is " << sVarofinterest << endl;
    cout << "v[0] is "  << v[0] << endl;
    if (v.size()==2) cout << "v[1] is "  << v[1]<< endl;
    else cout << "v[1] is not filled " << endl;
  }
  

  // get the desired directory

  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/24_3_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/30_3_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/1_4_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/15_4_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/NoDisplacedTrigger/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/CMSSW76V4/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/CMSSW76V4_NewCutFlow/";
  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/CMSSW76V4_TTLetp_10_8_2016/";


  CraneenPath=CraneenPath+channelpostfix;
  

  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
    {
      dataSetName = datasets[d]->Name();
      cout << " sVarofinterest is "  << sVarofinterest << endl;
      cout<<"Dataset:  :"<<dataSetName<<endl;

      filepath = CraneenPath+"/DisplacedTop_Run2_TopTree_Study_"+dataSetName + channelpostfix + ".root";  

      //filepath = CraneenPath+dataSetName+ ".root";
      if (debug) cout<<"filepath: "<<filepath<<endl;
	
      // get the tree corresponding to the final state of interest
      string stree = "";
      stree = "tree";


      // get different root file and tree name if using the trimmed trees.
      bool useTrimmedTree = false;
      if (useTrimmedTree)
	{
	  filepath = CraneenPath+"/DisplacedTop_Run2_TopTree_Study_"+dataSetName + channelpostfix + "SkimmedHighPt_OnZ_Lowd0.root";  

	  // add the correct suffix to select the corresponding region. (PCR, DCR, SR1, SR2, SR3)
	  //	  string prefix=stree+"PCR/";
	  //      string prefix=stree+"DCR/";
	  //string prefix=stree+"OnZ/";
	  
	  //	  stree = prefix+stree;
	  //	  stree = "treePCR/tree";
	}

		  
      FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
      string TTreename = stree;	

      cout << "TTreename is " << TTreename << endl;

      // change directory

      //      FileObj[dataSetName.c_str()]->cd("doubleElTreePCR");
	//      doubleElTreePCR->cd();
      //      cout << "current dir is "<< gDirectory->pwd() << endl;

      
      ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
      nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
      cout<<"                 nEntries: "<<nEntries<<endl;


      //      /*
      Int_t nMuons_;
      TBranch        *b_nMuons_;
      Int_t nElectrons_;
      TBranch        *b_nElectrons_;
//      */

      // Set the adress of the Branch that will be used in any case
      //      /*

      Int_t           nMuons;
      TBranch        *b_nMuons;

      Int_t           nElectrons;
      TBranch        *b_nElectrons; 


      ttree[dataSetName.c_str()]->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
      ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
      nMuons_=nMuons;
      nElectrons_=nElectrons;
      b_nMuons_= b_nMuons;
      b_nElectrons_=b_nElectrons;


      //      /*
      //      */
      


//      */

      // bo logic to set the right branch address depending on the string given as argument of the datasetplotter
      if (v.size() == 2){
	ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),v_varofInterest_double);
	ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&n_object); // this work and I am not sure why!
	//	ttree[dataSetName.c_str()]->SetBranchAddress("nMuons", &n_object, &b_nMuons); // this work but might depends on the last argument
	//	n_object=nMuons; // this work but this is not flexible
      }

      else if (v.size() == 1){
	if (v[0].compare(0,4,"evt_") == 0){
	  ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest_double);
	}
	else {
	  ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);
	}
	       
	
      }
      else {
	cout << "Vector of string does not have the good size!!!" << endl;
      }
      // eo logic to set the right branch address depending on the string given as argument of the datasetplotter
      
      

      bool isData= false;
      if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
      

      // get the Data lumi from the xml file

      if (isData) {
        DataEqLumi =  datasets[d]->EquivalentLumi() ; // pb-1 
        if (debug) cout << "DataEqLumi is " << DataEqLumi << endl;
      }
      Luminosity=DataEqLumi;

      
      ///////////////////////////////////////////
      /// Event Scale Factor ///////////////////
      ///////////////////////////////////////////

      
      // bo of event SF
      // -----------
      Double_t  puSF, globalScaleFactor, sf_muons, sf_electrons;


      if (applyGlobalSF && !isData){
	//	cout << "Applying Scale factor " << endl;
      }


      // get the SF from the corresponding branch
	

      //      Double_t sf_electron_[10];
      //      Double_t sf_muon_[10];
      //Double_t evt_puSF_;
      //TBranch  *b_evt_puSF_;

      // for muons final state

      //      /*

      Double_t sf_electron[10]; 
      TBranch        *b_sf_electron;

      Double_t sf_muon[10]; 
      TBranch        *b_sf_muon;

      Double_t        evt_puSF;
      TBranch        *b_evt_puSF;

      ttree[dataSetName.c_str()]->SetBranchAddress("sf_electron", sf_electron, &b_sf_electron); 
      ttree[dataSetName.c_str()]->SetBranchAddress("sf_muon", sf_muon, &b_sf_muon); 
      ttree[dataSetName.c_str()]->SetBranchAddress("evt_puSF", &evt_puSF, &b_evt_puSF);





      
      // -----------
      // eo of event SF
      // Declare the SF


      
      
      //      histo1D[dataSetName.c_str()] = new TH1F((dataSetName+"_"+v[0]).c_str(),(dataSetName+"_"+v[0]).c_str(), nBins, plotLow, plotHigh);

      // bo loop over the entries
      for (int j = 0; j<nEntries; j++)
      //      for (int j = 0; j<100; j++)
        {

	  ttree[(dataSetName).c_str()]->GetEntry(j);


	  puSF = globalScaleFactor = sf_muons = sf_electrons = 1.;
	  
	  //trick to avoid overwriting the number of el/mu in the events
	  //	  if(v.size() == 1 && sVarofinterest.find("nElectrons")!=string::npos) {varofInterest = nEl;}
	  //	  if(v.size() == 1 && sVarofinterest.find("nMuons")!=string::npos) {varofInterest = nMu;}

	  // bo of event SF
	  // -----------
	  
	  globalScaleFactor = 1.0;

	  //	  /*
	  // only apply individual SF if applyGlobalSF is true and sample is not data
	  if (applyGlobalSF && !isData){
	    
	    
	    // electron SF
	    if (applyElectronSF){
	      for (int i = 0; i < nElectrons_; i++)
		{
		  sf_electrons *=sf_electron[i];
		}

	      
	      globalScaleFactor *= sf_electrons;
	      if (debug){
		cout << "sf_electron is " << sf_electrons << endl;
		cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	      }
	    }
	    
	    // muon SF
	    if (applyMuonSF){
	      for (int i = 0 ; i < nMuons_ ; i++ )
		{
		  sf_muons *= sf_muon[i];
		}
	      globalScaleFactor *= sf_muons;
	      
	      if (debug){
		cout << "sf_muon is " << sf_muons << endl;
		cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	      }
	    }	
	    
	    // PU SF
	    if (applyPUSF){
	      puSF=evt_puSF;
	      globalScaleFactor *= puSF;

	      if (debug){
		cout << "puSF is " << puSF << endl;
		cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	      }

	    }
	    
	  }
	  //	  */

	  // eo of event SF
	  // -----------
	  

	  // -----------
	  // bo making MS plots

	  // make MS plot for single value
	  if (v.size() == 1){
	    if (v[0].compare(0,4,"evt_") == 0){ // store a double
	      if (isData) {
		MSPlot[plotname.c_str()]->Fill(varofInterest_double, datasets[d], false, 1); 
	      }
	      else { // not data 
		MSPlot[plotname.c_str()]->Fill(varofInterest_double, datasets[d], true, globalScaleFactor*Luminosity); 
	      }
	    }  
	    else { // not met
	      if (isData) {
		MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, 1);
	      }
	      else { // not data
		MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*Luminosity); 
	      }
	    }
	  }
	  
	  
	  // make MS plot for vector
	  if (v.size() == 2){

	    // bo of loop over the number of object per entry
	    for (int i_object =0 ; i_object < n_object ;i_object ++ )
	      {
		if (0)
		  {
		    cout << "i_object is " << i_object << endl;
		    cout << "n_object is " << n_object << endl;
		  }
		if (isData) 
		  {
		    // for data, fill once per event, weighted with the event scale factor 
		    if (plotAbs) MSPlot[plotname.c_str()]->Fill(fabs(v_varofInterest_double[i_object]), datasets[d], false, 1);
		    if (!plotAbs) MSPlot[plotname.c_str()]->Fill(v_varofInterest_double[i_object], datasets[d], false, 1);
		  }
		else
		  {
		    // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
		    if (plotAbs) MSPlot[plotname.c_str()]->Fill(fabs(v_varofInterest_double[i_object]), datasets[d], true, globalScaleFactor*Luminosity);
		    if (!plotAbs) MSPlot[plotname.c_str()]->Fill(v_varofInterest_double[i_object], datasets[d], true, globalScaleFactor*Luminosity);
		  }
	      
	      }

	  }


	  // eo making MS plots
	  // -----------


	}
      // eo loop over the entries

      
 
      TCanvas *canv = new TCanvas(("canv_"+v[0]+dataSetName).c_str(),("canv_"+v[0]+dataSetName).c_str());
      
      
      //      histo1D[dataSetName.c_str()]->Draw();
      string writename = "";
      if(isData)
	{
	  writename = "data_nominal";
	}
      else
	{
	  writename = dataSetName +"__nominal";
	}
      //cout<<"writename  :"<<writename<<endl;
      //      histo1D[dataSetName.c_str()]->Write((writename).c_str());
     
      //      canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
      //      canv->SaveAs((pathPNG+dataSetName+".C").c_str());
    }


  //  treeLoader.UnLoadDataset();
  
  if (debug){
    cout << "before cleaning" << endl;
    if (v.size() == 2){
      cout << " v[0] is " << v[0] << " and v[1] is " << v[1] << endl;
    }
    
    else if (v.size() == 1){
      cout << " v[0] is " << v[0] << endl;
      
    }
  }
  

  // clearing vector
  v.clear();
  if (debug){
    cout << "after cleaning" << endl ;
    cout << "v.size() is " << v.size() << endl;
  }
  

};

// TH2D function

void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX, int nBinsY, float lowY, float highY, string sVarofinterestY, string xmlNom, string TreePath )
{
  const char *xmlfile = xmlNom.c_str();
  cout << "entering TH2FPlotter" << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 
  treeLoader.LoadDatasets (datasets, xmlfile);
  if (debug) cout << "will start loading from xml file ..." << endl; 

  string plotname;
  plotname=sVarofinterestY+"Vs"+sVarofinterestX;
  


  // bo trial example
  histo2D[plotname.c_str()] =  new TH2F (plotname.c_str() ,plotname.c_str() ,nBinsX , lowX , highX , nBinsY, lowY, highY);


  int nEntries;
  string dataSetName, filepath ;
  int n_objectX = 0;
  double v_varofInterest_doubleX [20];
  int n_objectY = 0;
  double v_varofInterest_doubleY [20];
  char delim[] = " []";



  // get the the x axis
  vector<string> vX;
  // to avoid modifying original string                                               
  // first duplicate the original string and return a char pointer then free the memory                                                                                    
  
  char * dupX = strdup(sVarofinterestX.c_str());
  char * tokenX = strtok(dupX, delim);
  while(tokenX != NULL){
    vX.push_back(string(tokenX));
    // the call is treated as a subsequent calls to strtok:                           
    // the function continues from where it left in previous invocation               
    tokenX = strtok(NULL, delim);
  }
  free(dupX);

  if (debug){
    //    cout << "n_object is " << n_object << endl;                                 
    cout << "sVarofinterest is " << sVarofinterestX << endl;
    cout << "vX[0] is "  << vX[0] << endl;
    if (vX.size()==2) cout << "vX[1] is "  << vX[1]<< endl;
    else cout << "vX[1] is not filled " << endl;
  }


  // get the Y axis
  vector<string> vY;
  // to avoid modifying original string                                               
  // first duplicate the original string and return a char pointer then free the memory                                                                                    
  
  char * dupY = strdup(sVarofinterestY.c_str());
  char * tokenY = strtok(dupY, delim);
  while(tokenY != NULL){
    vY.push_back(string(tokenY));
    // the call is treated as a subsequent calls to strtok:                           
    // the function continues from where it left in previous invocation               
    tokenY = strtok(NULL, delim);
  }
  free(dupY);

  if (debug){
    //    cout << "n_object is " << n_object << endl;                                 
    cout << "sVarofinterest is " << sVarofinterestY << endl;
    cout << "vY[0] is "  << vY[0] << endl;
    if (vY.size()==2) cout << "vY[1] is "  << vY[1]<< endl;
    else cout << "vY[1] is not filled " << endl;
  }
  
  // defining x and y avis titles
  histo2D[plotname.c_str()]->GetXaxis()->SetTitle(vX[0].c_str());
  histo2D[plotname.c_str()]->GetYaxis()->SetTitle(vY[0].c_str());



  // get the desired directory


  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/2_3_2016/";
  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/14_3_2016/";


  CraneenPath=CraneenPath+channelpostfix;


  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets       
    {
      dataSetName = datasets[d]->Name();
      cout << " sVarofinterestX is "  << sVarofinterestX << endl;
      cout << " sVarofinterestY is "  << sVarofinterestY << endl;
      cout<<"Dataset:  :"<<dataSetName<<endl;
      //      filepath = CraneenPath+"/DisplacedTop_Run2_TopTree_Study_"+dataSetName + channelpostfix + ".root";
      //      filepath = CraneenPath+"/"+dataSetName + channelpostfix + "MuMuSkimmed.root";
      //filepath = CraneenPath+dataSetName+ ".root";                                                         
      if (debug) cout<<"filepath: "<<filepath<<endl;

      // get the tree corresponding to the final state of interest                                           
      string stree = "";
      if (DileptonMuMu) stree = "doubleMuTree";
      if (DileptonElEl) stree = "doubleElTree";
      if (DileptonElMu) stree = "tree";
      if (bbMu || bbEl) stree = "preCutTree";
      

      FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
      string TTreename = stree;
      cout << "TTreename is " << TTreename << endl;
      ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset              
      nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
      cout<<"                 nEntries: "<<nEntries<<endl;
      

      // set value for x 
      ttree[dataSetName.c_str()]->SetBranchAddress(vX[0].c_str(),v_varofInterest_doubleX);
      ttree[dataSetName.c_str()]->SetBranchAddress(vX[1].c_str(),&n_objectX);
      cout << "n_objectX is " << n_objectX << endl;

      // set value for y 
      ttree[dataSetName.c_str()]->SetBranchAddress(vY[0].c_str(),v_varofInterest_doubleY);
      ttree[dataSetName.c_str()]->SetBranchAddress(vY[1].c_str(),&n_objectY);



      // loop over the entries
      for (int j = 0; j<nEntries; j++)
        {

	  //	  histo2D[plotname.c_str()]->Fill(1.2,5);	  // dummy fill trial
	  ttree[(dataSetName).c_str()]->GetEntry(j);
	  // loop over the number of object 
	  // FIXME in this case the n_objectX is always equal to zero!!!
	  for (int i_object =0 ; i_object < n_objectX ;i_object ++ )// n_objectX = n_objectY
	    {      // Filling histo
	      histo2D[plotname.c_str()]->Fill(v_varofInterest_doubleX[i_object],v_varofInterest_doubleY[i_object]);
	    }
	}
    }
}


// function that writes all the MSPlots created in a root file
void MSPCreator (){
  Bool_t debug = false;

  string pathPNG = "myOutput";
  pathPNG += "_MSPlots/"+channelpostfix+"/"+regSuf+"/";
  mkdir(pathPNG.c_str(),0777);
  
  //  cout << "removing all the file of the directory " << pathPNG << endl;
  //  system("exec rm -r /"+pathPNG+"/");


  //  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  outfile->cd();



  // loop over TH2F
  /*
  for(map<string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
      cout << "entering the TH2F histo loop " << endl;
      string name = it->first;
      TH2F *temp = it->second;
      if (debug){
	cout << "Saving the MSP" << endl;
	cout << " and it->first is " << it->first << endl;
	cout << " and it->second is " << it->second << endl;
      }
      temp->Draw();
      temp->Write();
    }
  */

  
  
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      string name = it->first;
      MultiSamplePlot *temp = it->second;
      if (debug){
	cout << "Saving the MSP" << endl;
	cout << " and it->first is " << it->first << endl;
	cout << " and it->second is " << it->second << endl;
      }
      temp->Draw("MyMSP", 1, false, false, false, 10);
      temp->Write(outfile, "MyMSP"+it->first, true,"myOutput_MSPlots/"+channelpostfix+"/" , "png");
      temp->Write(outfile, "MyMSP"+it->first, true,"myOutput_MSPlots/"+channelpostfix+"/" , "pdf");
      //      vector<string> temp_histo = it->GetTH1FNames();
      //      for (int i_hist=0; i_hist < temp_histo.size();i_hist++  ){
      //	cout << "hist is" << temp_histo[i_hist] << endl;
      //	cout << "integral is " << it->GetTH1F(temp_histo[i_hist].GetSum()) << endl;
      //      }
    }



  
  outfile->Write("kOverwrite");
}


// function that converts an int into a string
std::string intToStr (int number){
  std::ostringstream buff;
  buff<<number;
  return buff.str();
}

