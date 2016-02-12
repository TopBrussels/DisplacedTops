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


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath);
void MSPCreator ();

void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );


// faco TO BE CHANGED
Bool_t debug = false;
Bool_t debug_plot = false;
bool DileptonElMu = false;
bool DileptonMuMu = false;
bool DileptonElEl = false;
string channelpostfix = "";
//applying all appropriate scale factors for individual objects if the bool is set to true
Bool_t applyElectronSF = false; 
Bool_t applyMuonSF = false; 
Bool_t applyPUSF = false; 
Bool_t applyGlobalSF = false; 



int main(int argc, char* argv[])
{
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
  applyElectronSF = strtol(argv[3],NULL,10);
  applyMuonSF = strtol(argv[4],NULL,10);
  applyPUSF = strtol(argv[5],NULL,10);
  applyGlobalSF = strtol(argv[6],NULL,10);


  /*
    Double_t aDouble;
    const char aDoubleType =  aDouble.TDataMember::GetTypeName();
    cout << "aDoubleType is " << aDoubleType << endl; //faco
  */
  
  

  string xmlFileName;

  if(channel=="ElMu" || channel== "MuEl")
    {
      cout << " --> Using the Muon-Electron channel..." << endl;
      channelpostfix = "_MuEl";
      //      xmlFileName = "config/FullSamplesElMuV0TreeProc.xml";
      xmlFileName = "config/TreeProc_FullSamplesElMuV0.xml";
      DileptonElMu=true;
    }
  else if(channel=="MuMu")
    {
      cout << " --> Using the Muon-Muon channel..." << endl;
      channelpostfix = "_MuMu";
      //      xmlFileName = "config/FullSamplesMuMuV0TreeProc.xml";
      xmlFileName = "config/TreeProc_FullSamplesMuMuV0.xml";
      DileptonMuMu=true;
      //      cout << "DileptonMuMu is " << DileptonMuMu << endl;
    }
  else if(channel=="ElEl")
    {
      cout << " --> Using the Electron-Electron channel..." << endl;
      channelpostfix = "_ElEl";
      //      xmlFileName = "config/FullSamplesElElV0TreeProc.xml";
      xmlFileName = "config/TreeProc_FullSamplesElElV0.xml";
      DileptonElEl=true;
    }
  else
    {
      cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
      exit(1);
    }
  
  string CraneenPath;

  
  if (debug_plot){
    //    xmlFileName = "config/FullSamplesMuMuV0TreeProc.xml";
    xmlFileName = "config/TreeProc_FullSamplesMuMuV0.xml";

    // only few plots!
    if (DileptonMuMu) {
      //      DatasetPlotter(5, -0.5, 4.5, "nMuons_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(70, -0.5, 69.5, "nvtx_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(40, 0, 800, "evt_met_mumu", xmlFileName,CraneenPath);
    }
    if (DileptonElEl){
      //      DatasetPlotter(5, -0.5, 4.5, "nElectrons_elel", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
    }

  }
    
    
  else{
    
    //    xmlFileName = "config/FullSamplesMuMuV9TreeProc.xml";
     

    //    cout << "CraneenPath is " << CraneenPath << endl;
    cout << "xmlFileName is " << xmlFileName << endl;



    // calling datasetPlotter to create MSPplots



    // bo selecting the right plots and/or variables depending on the final state


    // E-MU FINAL STATE
    if (DileptonElMu){
          
      // event plots
      //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
      DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath);
    

      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_electron[nElectrons]", xmlFileName,CraneenPath);
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

      DatasetPlotter(70, -0.5, 69.5, "nvtx_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(40, 0, 800, "evt_met_mumu", xmlFileName,CraneenPath);
          
     
      /*
      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      */
      
 
      //      /*
      // muon plots
      DatasetPlotter(11, -0.5, 10.5, "nMuons_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      //	DatasetPlotter(20, 10, 50, "pt_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -10, 10, "vz_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);

      // muonPairs plots
      DatasetPlotter(100, 0.0, 1.5, "deltaVz_mumu[nMuonPairs_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 1000, "invMass_mumu[nMuonPairs_mumu]", xmlFileName,CraneenPath);
      //      */

      // electron-muon plots
      //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);
    }

    // El-El FINAL STATE
    if (DileptonElEl){
          
      // event plots
      //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
      DatasetPlotter(70, -0.5, 69.5, "nvtx_elel", xmlFileName,CraneenPath);
      DatasetPlotter(40, 0, 800, "evt_met_elel", xmlFileName,CraneenPath);
              

      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons_elel", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      //	DatasetPlotter(20, 10, 50, "pt_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -10, 10, "vz_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron_elel[nElectrons_elel]", xmlFileName,CraneenPath);

      
      // Dielectron plots
      DatasetPlotter(100, 0.0, 1.5, "deltaVz_elel[nElectronPairs_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 1000, "invMass_elel[nElectronPairs_elel]", xmlFileName,CraneenPath);
      
      /*
      // muon plots
      //    DatasetPlotter(11, -0.5, 10.5, "nMuons_elel", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_muon_elel[nMuons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_muon_elel[nMuons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_muon_elel[nMuons_elel]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon_elel[nMuons_elel]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_muon_elel[nMuons_elel]", xmlFileName,CraneenPath);
      */

      // electron-muon plots
      //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);

    }
  }

  // eo selecting the right plots and/or variables depending on the final state  


  // calling the function that writes all the MSPlots in a root file
  MSPCreator ();

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath)
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
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/3_2_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/4_2_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/5_2_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/7_2_2016/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/9_2_2016/";
  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_6_3/src/TopBrussels/DisplacedTops/MergedTrees/10_2_2016/";

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
      if (DileptonMuMu) stree = "doubleMuTree";
      if (DileptonElEl) stree = "doubleElTree";
      if (DileptonElMu) stree = "tree";
		  
      FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
      string TTreename = stree;	




      ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
      nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
      cout<<"                 nEntries: "<<nEntries<<endl;
      

      // bo logic to set the right branch address depending on the string given as argument of the datasetplotter
      if (v.size() == 2){
	ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&n_object); 
	ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),v_varofInterest_double);
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
      
      
      Int_t data_index;
      if (isData) {
        Luminosity =  datasets[d]->EquivalentLumi() ; // pb-1 
        cout << "Luminosity is " << Luminosity << endl;
	data_index=d;
	
      }

      //Luminosity =datasets[data_index]->EquivalentLumi() ;
      Luminosity=2610.49;


      
      ///////////////////////////////////////////
      /// Event Scale Factor ///////////////////
      ///////////////////////////////////////////

      
      // bo of event SF
      // -----------
      Double_t  puSF, globalScaleFactor, sf_muon, sf_electron;


      /*
	Bool_t applyElectronSF, applyMuonSF , applyPUSF, applyGlobalSF;
      
	if (applyGlobalSF && !isData){
	cout << "Applying Scale factor " << endl;
	}

	// get the SF from the corresponding branch

	Int_t nEl;
	ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons_elel",&nEl);
	Int_t nMu;
	ttree[dataSetName.c_str()]->SetBranchAddress("nMuons_elel",&nMu);

	Double_t electronSF[10];
	ttree[dataSetName.c_str()]->SetBranchAddress("sf_electron_elel",&electronSF);
	Double_t muonSF[10];
	ttree[dataSetName.c_str()]->SetBranchAddress("sf_muon_elel", &muonSF);

      */
      // -----------
      // eo of event SF
      // Declare the SF


      
      
      //      histo1D[dataSetName.c_str()] = new TH1F((dataSetName+"_"+v[0]).c_str(),(dataSetName+"_"+v[0]).c_str(), nBins, plotLow, plotHigh);

      // bo of loop through entries and fill plots
      for (int j = 0; j<nEntries; j++)
        {

	  ttree[(dataSetName).c_str()]->GetEntry(j);


	  puSF = globalScaleFactor = sf_muon = sf_electron = 1.;
	  
	  //trick to avoid overwriting the number of el/mu in the events
	  //	  if(v.size() == 1 && sVarofinterest.find("nElectrons_elel")!=string::npos) {varofInterest = nEl;}
	  //	  if(v.size() == 1 && sVarofinterest.find("nMuons_elel")!=string::npos) {varofInterest = nMu;}

	  // bo of event SF
	  // -----------
	  
	  globalScaleFactor = 1.0;

	  /*
	  // only apply individual SF if applyGlobalSF is true and sample is not data
	  if (applyGlobalSF && !isData){
	    
	  puSF = globalScaleFactor = sf_muon = sf_electron = 1;
	    
	  // electron SF
	  if (applyElectronSF){
	  for (int i = 0; i < nEl; i++)
	  {
	  sf_electron *=electronSF[i];
	  }

	      
	  globalScaleFactor *= sf_electron;
	  if (debug){
	  cout << "sf_electron is " << sf_electron << endl;
	  cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	  }
	  }
	    
	  // muon SF
	  if (applyMuonSF){
	  sf_muon=0.5;

	  globalScaleFactor *= sf_muon;
	  if (debug){
	  cout << "sf_muon is " << sf_muon << endl;
	  cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	  }
	  }	
	    
	  // PU SF
	  if (applyPUSF){
	  globalScaleFactor *= puSF;
	  if (debug){
	  cout << "puSF is " << puSF << endl;
	  cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	  }

	  }
	    
	  }
	  */


	  // -----------
	  // eo of event SF


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
		    MSPlot[plotname.c_str()]->Fill(v_varofInterest_double[i_object], datasets[d], false, 1); 	
		  }
		else
		  {
		    // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
		    MSPlot[plotname.c_str()]->Fill(v_varofInterest_double[i_object], datasets[d], true, globalScaleFactor*Luminosity);
		  }
	      
	      }

	  }

	  //		  MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, globalScaleFactor); 

	}
      
 
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


// function that writes all the MSPlots created in a root file
void MSPCreator (){
  Bool_t debug = false;

  string pathPNG = "myOutput";
  pathPNG += "_MSPlots/"+channelpostfix+"/";
  mkdir(pathPNG.c_str(),0777);
  //  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  outfile->cd();
  
  
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      string name = it->first;
      MultiSamplePlot *temp = it->second;
      if (debug){
	cout << "Saving the MSP" << endl;
	cout << " and it->first is " << it->first << endl;
      }
      temp->Draw("MyMSP", 1, false, false, false, 10);
      temp->Write(outfile, "MyMSP"+it->first, true,"myOutput_MSPlots/"+channelpostfix+"/" , "png");
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

