#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"

#include <sstream>

using namespace std;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;

std::string intToStr (int number);

void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string CraneenPath);




int main()
{
    int NumberOfBins = 3;	//fixed width nBins

    //------- Set Channel --------//
    bool DileptonMuEl = false;
    bool DileptonMuMu = false;
    bool DileptonElEl = false;
    bool SingleMu = true;
    bool SingleEl = false;
    bool jetSplit = false;
    bool jetTagsplit = true;

    
    //    string VoI = "BDT"; //variable of interest for plotting
    //    string VoI = "nElectrons"; //variable of interest for plotting
    float lBound = 0;   //-1->0.2 topness
    float uBound = 3;
    int lumiScale = 50;  //Amount of luminosity to scale to in fb^-1

    /*
    vector<string> vars;
    vars.push_back("HT");
    vars.push_back("LeadingMuonPt");
    vars.push_back("LeadingElectronPt");
    vars.push_back("LeadingBJetPt");
    vars.push_back("HT2M");
    vars.push_back("MVAvals1");
    */

    string xmlFileName;
    string xmlFileNameSys;
    string CraneenPath;
    string splitVar, splitVar1, splitVar2;
    string splitting;

    float bSplit, tSplit, wSplit, bSplit1, tSplit1, wSplit1, bSplit2, tSplit2, wSplit2;  //the bottom, top, and width of the splitting

    if(SingleMu)
    {

        xmlFileName = "config/Run2SingleLepton_samples.xml";
        xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
	//	CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens24_3_2015_merge/Craneen_";// faco here we set the craneen path
	CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/";
    }

    std::string slumiScale = intToStr(lumiScale);

    if(jetSplit == true) splitting = "JS";
    else if (jetTagsplit == true)  splitting = "JTS";
    else{splitting = "inc";}



    cout << "HERE!!!!" << endl;
    //        SystematicsAnalyser(NumberOfBins, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, xmlFileNameSys, CraneenPath);
    cout << "CraneenPath is " << CraneenPath << endl;
    cout << "xmlFileName is " << xmlFileName << endl;
    //    DatasetPlotter(NumberOfBins, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath);
    //DatasetPlotter(10, 0, 10, "nElectrons", "config/Run2SingleLepton_samples.xml","/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/");//xmlFileName, CraneenPath);
    DatasetPlotter(100, 0, 1000, "pt_electron[nElectrons]", "config/Run2SingleLepton_samples.xml","/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/");//xmlFileName, CraneenPath);
    DatasetPlotter(100, 0, 1000, "pt_muon[nMuons]", "config/Run2SingleLepton_samples.xml","/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/");//xmlFileName, CraneenPath);
    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", "config/Run2SingleLepton_samples.xml","/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/");//xmlFileName, CraneenPath);


}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;


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
    TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
    outfile->cd();
    string plotname = sVarofinterest;   ///// Non Jet Split plot
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity;
    int varofInterest;
    double varofInterest_double [20];

    string varofinterest_iterator;
    string varofinterest_fixed;



    //    string sdelim = " []";
    //    char* delim =  " []";
    
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

    cout << v[0] << "  " << v[1] << endl;



    //    cout << varofinterest_fixed << " " << varofinterest_iterator << endl;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

	//        filepath = CraneenPath+dataSetName + "_Run2_TopTree_Study.root";
        filepath = CraneenPath+dataSetName + ".root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset

	//        string nTuplename = "Craneen__"+ leptoAbbr;
	//        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
	//        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();

	string TTreename = "tree";	
        ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
	nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

	

	//        ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
	//        ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
	//        ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
	//        ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);

	//        ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&varofInterest); 
        ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),varofInterest_double);
	//        ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest_double.c_str(),&varofInterest_double);
	//	ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons",&nElectrons);

	bool isData= false;
	if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
	  

	ScaleFactor = 1;
	Luminosity =1;
	
        //for fixed bin width
        histo1D[dataSetName.c_str()] = new TH1F((dataSetName+"_"+v[0]).c_str(),(dataSetName+"_"+v[0]).c_str(), nBins, plotLow, plotHigh);
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            ttree[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
	    //	    cout << "Var of interest is " << varofInterest << endl;
	    
	    for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
	      {
	      if (isData)
		{
		  MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, NormFactor*ScaleFactor*Luminosity);
		  histo1D[dataSetName.c_str()]->Fill(varofInterest_double[i_object],NormFactor*ScaleFactor*Luminosity);
		}
	      else
		{
		  
		  MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, ScaleFactor*Luminosity);
		  histo1D[dataSetName.c_str()]->Fill(varofInterest_double[i_object],NormFactor*ScaleFactor*Luminosity);
		}
	      }

        }

 
        TCanvas *canv = new TCanvas(("canv"+v[0]).c_str(),("canv"+v[0]).c_str());


        histo1D[dataSetName.c_str()]->Draw();
        string writename = "";
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            writename = "data_nominal";
        }
        else
        {
            writename = dataSetName +"__nominal";
        }
        //cout<<"writename  :"<<writename<<endl;
        histo1D[dataSetName.c_str()]->Write((writename).c_str());

        canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
        canv->SaveAs((pathPNG+dataSetName+".C").c_str());
    }


    treeLoader.UnLoadDataset();


    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 100);
	temp->Write(outfile, "MyMSP"+v[0], false,"myOutput_MSPlots" , "png");
    }
    outfile->Write("kOverwrite");
};



std::string intToStr (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

