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



int main()
{
    int NumberOfBins = 3;	//fixed width nBins


    //------- Set Channel --------//
    bool DileptonMuEl = false;
    bool DileptonMuMu = false;
    bool DileptonElEl = false;
    bool SingleMu = true;
    bool SingleEl = false;
    float lBound = 0;   //-1->0.2 topness
    float uBound = 3;
    float lumiScale = -1;
    //    int lumiScale = 0.10651;  //Amount of luminosity to scale to in fb^-1
    //    int lumiScale = 0.10651;  //Amount of luminosity to scale to in fb^-1
    string xmlFileName;
    string xmlFileNameSys;
    string CraneenPath;
    string splitting;


    if(1)
    {
      xmlFileName = "config/FullSamplesV8TreeProc.xml";
     

      //      CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/";

      //      CraneenPath="/user/qpython/TopBrussels7X/CMSSW_7_4_12_patch1/src/TopBrussels/DisplacedTops/Craneens_MuEl/Craneens29_9_2015/";
    }

    //    cout << "CraneenPath is " << CraneenPath << endl;
    cout << "xmlFileName is " << xmlFileName << endl;



    // calling datasetPlotter to create MSPplots

    // event plots
    //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath);
    DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath);
    DatasetPlotter(100, 0, 5, "puSF", xmlFileName,CraneenPath);


    

    // electron plots
    //    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath);
    DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", xmlFileName,CraneenPath);
    DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath);
    DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath);
    //    DatasetPlotter(100, -0.1, 0.1, "d0_electron[nElectrons]", xmlFileName,CraneenPath);
    DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron[nElectrons]", xmlFileName,CraneenPath);
    DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath);
    
 
    // muon plots
    //    DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath);
    DatasetPlotter(100, 0, 1000, "pt_muon[nMuons]", xmlFileName,CraneenPath);
    DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath);
    DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath);
    //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
    DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath);
    DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath);


    // electron-muon plots
    //    DatasetPlotter(100, 0.0, 0.2, "to_Add_Invmass", xmlFileName,CraneenPath);






    // calling the function that writtes all the MSPlots in a root file
    MSPCreator ();

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath)
{
  Bool_t debug = false;
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
  pathPNG += "_MSPlots/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  
    //***************************************************CREATING PLOTS****************************************************
  //  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  //  outfile->cd();
  string plotname = sVarofinterest;   ///// Non Jet Split plot
  // make for loop here!!!
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str()); 

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;
  float ScaleFactor, NormFactor, Luminosity;
  int varofInterest;
  double varofInterest_double [20];


  
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



  // declaring list of variables that contains the SF info
  double_t sf_muon, sf_electron, PUSF;


  //cout << v[0] << "  " << v[1] << endl;
  

  //  string CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_4_12_patch1/src/TopBrussels/DisplacedTops/Craneens_MuEl/Craneens29_9_2015/";
  //  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops/MergedTrees/19_11_2015/";
  TString CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/DisplacedTops/MergedTrees/25_11_2015/";

  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
    {
      dataSetName = datasets[d]->Name();
      cout<<"Dataset:  :"<<dataSetName<<endl;
      filepath = CraneenPath+"/DisplacedTop_Run2_TopTree_Study_"+dataSetName +"_MuEl"+ ".root";
      //filepath = CraneenPath+dataSetName+ ".root";
      if (debug) cout<<"filepath: "<<filepath<<endl;
	

      FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
      string TTreename = "tree";	
      ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
      nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
      cout<<"                 nEntries: "<<nEntries<<endl;
      
      
      // bo logic to set the right branch address depending on the string given as argument of the datasetplotter
      if (v.size() == 2){
	ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&varofInterest); 
	ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),varofInterest_double);
      }

      else if (v.size() == 1){
	//	if (debug)	cout << "v.size is to 1" << " and v[0] is " << v[0] << endl ;
	ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);//&varofInterest
	
      }
      else {
	cout << "Vector of string does not have the good size!!!" << endl;
      }
      // eo logic to set the right branch address depending on the string given as argument of the datasetplotter
      
      

      bool isData= false;
      if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
      



      
      ///////////////////////////////////////////                                                                    
      //  Event Scale Factor                                                                                         
      ///////////////////////////////////////////                                                                    
      
      // bo of event SF
      // -----------

      //applying all appropriate scale factors for individual objects if the bool is set to true
      Bool_t applyElectronSF, applyMuonSF , applyPUSF, applyGlobalSF;
      applyElectronSF = true;
      applyMuonSF = true;
      applyPUSF = true;
      applyGlobalSF = true;
      
      if (applyGlobalSF && !isData){
	cout << "Applying Scale factor " << endl;
      }

      // get the SF from the corresponding branch
      Double_t sf_electron, sf_muon, puSF, globalScaleFactor;
      ttree[dataSetName.c_str()]->SetBranchAddress("sf_muon",&sf_muon);
      ttree[dataSetName.c_str()]->SetBranchAddress("sf_electron",&sf_electron);
      ttree[dataSetName.c_str()]->SetBranchAddress("puSF",&puSF); // change to sf_pu

      // -----------
      // eo of event SF
      // Declare the SF



      Luminosity = 552.672; // pb-1
      TString slumi = "552.672";

      
      //      histo1D[dataSetName.c_str()] = new TH1F((dataSetName+"_"+v[0]).c_str(),(dataSetName+"_"+v[0]).c_str(), nBins, plotLow, plotHigh);

      // bo of loop through entries and fill plots
      for (int j = 0; j<nEntries; j++)
        {
	  ttree[(dataSetName).c_str()]->GetEntry(j);

	  // bo of event SF
	  // -----------
	  
	  globalScaleFactor = 1.0;

	  // only apply individual SF if applyGlobalSF is true and sample is not data
	  if (applyGlobalSF && !isData){
	  
	    // electron SF
	    if (applyElectronSF){
	      globalScaleFactor *= sf_electron;
	      if (debug){
		cout << "sf_electron is " << sf_electron << endl;
		cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	      }
	    }
	    
	    // muon SF
	    if (applyMuonSF){

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

	  // -----------
	  // eo of event SF



	  //	    cout << "Var of interest is " << varofInterest << endl;
	  
	  // bo of loop over the number of object per entry
	  for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
	    {
	      if (debug) cout << "varofInterest is " << varofInterest << endl;
	      if (isData) 
		{
		  // for data, fill once per event, weighted with the event scale factor only
		  MSPlot[plotname.c_str()]->Fill(varofInterest/*_double[i_object]*/, datasets[d], false, globalScaleFactor); 
		}
	      else
		{
		  // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
		  MSPlot[plotname.c_str()]->Fill(varofInterest/*_double[i_object]*/, datasets[d], true, globalScaleFactor*Luminosity);
		}
	      
	    }

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
  pathPNG += "_MSPlots/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
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
      temp->Draw("MyMSP", 1, false, false, false, 100);
      temp->Write(outfile, "MyMSP"+it->first, true,"myOutput_MSPlots" , "png");
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

