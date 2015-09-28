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

void SystematicsAnalyser(int nBins, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath);

void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath);
void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath);

void Split2DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, string CraneenPath);
void Split2SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlSys, string CraneenPath);
 
void Split2_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, bool jetTagsplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom);

int main()
{
    int NumberOfBins = 30;	//fixed width nBins

    //------- Set Channel --------//
    bool DileptonMuEl = false;
    bool DileptonMuMu = false;
    bool DileptonElEl = false;
    bool SingleMu = true;
    bool SingleEl = false;
    bool jetSplit = false;
    bool jetTagsplit = true;

    string VoI = "BDT"; //variable of interest for plotting
    float lBound = -0.5;   //-1->0.2 topness
    float uBound = 1.0;
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
    string leptoAbbr;
    string channel;
    string xmlFileName;
    string xmlFileNameSys;
    string CraneenPath;
    string splitVar, splitVar1, splitVar2;
    string splitting;

    float bSplit, tSplit, wSplit, bSplit1, tSplit1, wSplit1, bSplit2, tSplit2, wSplit2;  //the bottom, top, and width of the splitting

    if(SingleMu)
    {
        leptoAbbr = "Mu";
        channel = "ttttmu";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
        xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
	//	CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens24_3_2015_merge/Craneen_";// faco here we set the craneen path
	CraneenPath = "/user/qpython/TopBrussels7X/CMSSW_7_5_3/src/TopBrussels/DisplacedTops/MACRO_Output_MuEl/";
    }

    std::string slumiScale = intToStr(lumiScale);

    if(jetSplit == true) splitting = "JS";
    else if (jetTagsplit == true)  splitting = "JTS";
    else{splitting = "inc";}

    string shapefileName = "";
    shapefileName = "shapefile"+leptoAbbr+"_"+slumiScale+"_"+VoI+"_"+splitting+".root";
    TFile *shapefile = new TFile((shapefileName).c_str(), "RECREATE");
    TFile *errorfile = new TFile(("ScaleFiles"+leptoAbbr+"_light/Error.root").c_str(),"RECREATE");


    cout << "HERE!!!!" << endl;
    //        SystematicsAnalyser(NumberOfBins, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, xmlFileNameSys, CraneenPath);
    cout << "CraneenPath is " << CraneenPath << endl;
    cout << "xmlFileName is " << xmlFileName << endl;
    DatasetPlotter(NumberOfBins, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath);


    errorfile->Close();
    shapefile->Close();
    delete shapefile;
    delete errorfile;
}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname = sVarofinterest;   ///// Non Jet Split plot
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

	//        filepath = CraneenPath+dataSetName + "_Run2_TopTree_Study.root";
        filepath = CraneenPath+dataSetName + ".root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset

	//        string nTuplename = "Craneen__"+ leptoAbbr;
	//        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
	//        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();

	string TTreename = "tree";	
        ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
	nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);



        //for fixed bin width
        histo1D[dataSetName.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            ttree[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
            Luminosity = 15000;

            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
                histo1D[dataSetName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
            else
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
                histo1D[dataSetName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }

        }
        if(dataSetName == "TTJets")  //to put nominal histo into error file
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[dataSetName.c_str()]->Write("Nominal");
            errorfile->Write();
        }

        shapefile->cd();
        TCanvas *canv = new TCanvas();

        histo1D[dataSetName.c_str()]->Draw();
        string writename = "";
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            writename = channel + "__data_obs__nominal";
        }
        else
        {
            writename = channel + "__" + dataSetName +"__nominal";
        }
        //cout<<"writename  :"<<writename<<endl;
        histo1D[dataSetName.c_str()]->Write((writename).c_str());

        canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
    }


    treeLoader.UnLoadDataset();

    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error.root";
    MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 100);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
    }
};


void SystematicsAnalyser(int nBins, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING SYS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"	Dataset:  :"<<dataSetName<<endl;
        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

        cout<<""<<endl;
        cout<<"file in use: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
        cout<<"initialised file"<<endl;
        string nTupleName = "Craneen__"+leptoAbbr;

         ttree[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
        nEntries = (int) ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        //sVarofinterest = "HT";
         ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
         ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);



        ////****************************************************Define plots***********************************************
        string plotname = sVarofinterest+"_"+dataSetName;
        histo1D[plotname.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);

        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
             ttree[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi
            Luminosity = 15000;
            histo1D[plotname.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
        }

        ////****************************************************Fill plots***********************************************
        shapefile->cd();

        TCanvas *canv2 = new TCanvas();
        if(Normalise)
        {
            double dIntegral = histo1D[plotname.c_str()]->Integral();
            histo1D[plotname.c_str()]->Scale(1./dIntegral);
        }
        histo1D[plotname.c_str()]->Draw();
        string writename = "";
        writename = channel + "__TTJets__" + dataSetName;
        //cout<<"writename  :"<<writename<<endl;

        histo1D[plotname.c_str()]->Write((writename).c_str());
        canv2->SaveAs(("Sys_"+plotname+".pdf").c_str());

        if(dataSetName == "TTScaleDown")
        {
            errorfile->cd();
            errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[plotname.c_str()]->Write("Minus");
            //errorfile->Write();
        }

        if(dataSetName == "TTScaleUp")
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[plotname.c_str()]->Write("Plus");
            //errorfile->Write();
        }
    }
};

void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname;   ///// Jet Split plot
    string numStr;
    string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
    string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
    string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
    for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
    {
        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        plotname = sVarofinterest + numStr + sSplitVar;
        MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
    }
    plotname = "";

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath, histoName;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
         ttree[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int) ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

         ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
         ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);



        //for fixed bin width
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
        }
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
             ttree[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi

            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }

            
            if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
            {
                plotname = sVarofinterest + stSplit + sSplitVar;
                histoName = dataSetName + stSplit + sSplitVar;
                //outsiderange = false;
            }
            else //If it doesn't belong in the last bin, find out which it belongs in
            {
                bool outsiderange = true;
                for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
                {
                    if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
                    {
                        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                        plotname = sVarofinterest + numStr + sSplitVar;
                        histoName = dataSetName + numStr + sSplitVar;
                        outsiderange = false;
                        break;
                    }
                }
                if (outsiderange == true){
                    continue;
                }
            }


            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
            else
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }

        }
        if(dataSetName == "TTJets")  //to put nominal histo into error file
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[histoName.c_str()]->Write("Nominal");
            errorfile->Write();
        }

        shapefile->cd();
        TCanvas *canv = new TCanvas();
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            plotname = sVarofinterest + numStr + sSplitVar;
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()]->Draw();
            string writename = "";
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                writename = channel + numStr + sSplitVar + "__data_obs__nominal";
            }
            else
            {
                writename = channel + numStr + sSplitVar + "__" + dataSetName +"__nominal";
            }
            //cout<<"writename  :"<<writename<<endl;
            histo1D[histoName.c_str()]->Write((writename).c_str());

            canv->SaveAs((pathPNG+histoName+".pdf").c_str());
        }
    }


    treeLoader.UnLoadDataset();

    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error.root";
    //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
        temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 100);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
    }
};


void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING SYS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath, numStr, histoName;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar;
    string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
    string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
    string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
    errorfile->cd();
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"	Dataset:  :"<<dataSetName<<endl;
        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

        cout<<""<<endl;
        cout<<"file in use: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
        cout<<"initialised file"<<endl;
        string nTupleName = "Craneen__"+leptoAbbr;

         ttree[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
        nEntries = (int) ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        //sVarofinterest = "HT";
         ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
         ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);



        ////****************************************************Define plots***********************************************
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
        }

        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
             ttree[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi

            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }
            
            //bool outsiderange = true;
            if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
            {
                histoName = dataSetName + stSplit + sSplitVar;
                //outsiderange = false;
            }
            else //If it doesn't belong in the last bin, find out which it belongs in
            {
                bool outsiderange = true;
                for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
                {
                    if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
                    {
                        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                        histoName = dataSetName + numStr + sSplitVar;
                        outsiderange = false;
                        break;
                    }
                }
                if (outsiderange == true){
                    continue;
                }
            }


            histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
        }

        ////****************************************************Fill plots***********************************************


        TCanvas *canv2 = new TCanvas();

        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            shapefile->cd();

            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            if(Normalise)
            {
                double dIntegral = histo1D[histoName.c_str()]->Integral();
                histo1D[histoName.c_str()]->Scale(1./dIntegral);
            }
            histo1D[histoName.c_str()]->Draw();
            string writename = "";

//            writename = channel + numStr + sSplitVar + "__" + dataSetName +"__nominal";
//
//            cout<<"writename  :"<<writename<<endl;
//            histo1D[histoName.c_str()]->Write((writename).c_str());
//
//            canv->SaveAs(("Sys_"+histoName+".pdf").c_str());

            if(dataSetName == "TTScaleDown")
            {
                writename = channel + numStr + sSplitVar + "__TTJets__scaleDown";
                //cout<<"writename  :"<<writename<<endl;
                histo1D[histoName.c_str()]->Write((writename).c_str());
                canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                errorfile->cd();
                //errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
                errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                histo1D[histoName.c_str()]->Write("Minus");
                //errorfile->Write();
            }

            if(dataSetName == "TTScaleUp")
            {
                writename = channel + numStr + sSplitVar + "__TTJets__scaleUp";
                //cout<<"writename  :"<<writename<<endl;
                histo1D[histoName.c_str()]->Write((writename).c_str());
                canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                errorfile->cd();
                errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                histo1D[histoName.c_str()]->Write("Plus");
                //errorfile->Write();
            }
        }


//        histo1D[plotname.c_str()]->Draw();
//        string writename = "";
//        writename = channel + "__TTJets__" + dataSetName;
//        cout<<"writename  :"<<writename<<endl;
//
//        histo1D[plotname.c_str()]->Write((writename).c_str());
//        canv2->SaveAs(("Sys_"+plotname+".pdf").c_str());
//
//        if(dataSetName == "TTScaleDown")
//        {
//            errorfile->cd();
//            errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//            histo1D[plotname.c_str()]->Write("Minus");
//            //errorfile->Write();
//        }
//
//        if(dataSetName == "TTScaleUp")
//        {
//            errorfile->cd();
//            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//            histo1D[plotname.c_str()]->Write("Plus");
//            //errorfile->Write();
//        }
    }
};

void Split2DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname;   ///// Jet Split plot
    string numStr1;
    string numStr2;
    string sbSplit1 = static_cast<ostringstream*>( &(ostringstream() << fbSplit1) )->str();
    string stSplit1 = static_cast<ostringstream*>( &(ostringstream() << ftSplit1) )->str();
    string swSplit1 = static_cast<ostringstream*>( &(ostringstream() << fwSplit1) )->str();
    string sbSplit2 = static_cast<ostringstream*>( &(ostringstream() << fbSplit2) )->str();
    string stSplit2 = static_cast<ostringstream*>( &(ostringstream() << ftSplit2) )->str();
    string swSplit2 = static_cast<ostringstream*>( &(ostringstream() << fwSplit2) )->str();
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
            plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
        }  
    }
    plotname = "";

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath, histoName;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar1, splitVar2;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
         ttree[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int) ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

         ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
         ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);

        //for fixed bin width
        for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
        {
            numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
                numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();   
                histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
            }
        }
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
             ttree[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }
            bool outsiderange = true;
            if(splitVar1 >= ftSplit1) //Check if this entry belongs in the last bin in var1.  Done here to optimize number of checks
            {   

                if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
                    plotname = sVarofinterest + stSplit1 + sSplitVar1 +stSplit2 + sSplitVar2;
                    histoName = dataSetName + stSplit1 + sSplitVar1 +stSplit2 + sSplitVar2;
                    outsiderange = false;
                    //                                cout<<"splitvar2: "<<splitVar2<<endl;
                }    
                else //If it doesn't belong in the last bin in var2, find out which it belongs in
                {

                    for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
                    {
                        if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
                        {
                            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                            plotname = sVarofinterest + stSplit1 + sSplitVar1 + numStr2 + sSplitVar2;
                            histoName = dataSetName + stSplit1 + sSplitVar1 + numStr2 + sSplitVar2;
                            outsiderange = false;
                            break;
                        }
                    }
                }            
            }
            else //If it doesn't belong in the last bin, find out which it belongs in
            {
                for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
                {
                    if(splitVar1>=s && splitVar1<(s+fwSplit1)) //splitVar falls inside one of the bins
                    {
                        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                        if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
                            plotname = sVarofinterest + numStr1 + sSplitVar1 +stSplit2 + sSplitVar2;
                            histoName = dataSetName + numStr1 + sSplitVar1 +stSplit2 + sSplitVar2;
                            outsiderange = false;
                        }    
                        else //If it doesn't belong in the last bin, find out which it belongs in
                        {
                            for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
                            {
                                if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
                                {
                                    numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                                    plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                                    histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                                    outsiderange = false;
                                    break;
                                }
                            }
                        } 
                        //plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                        //histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;;
                        //break;
                    }
                }
            }
            if (outsiderange == true){
                continue;
            }
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
            else
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
        }
        if(dataSetName == "TTJets")  //to put nominal histo into error file
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[histoName.c_str()]->Write("Nominal");
            errorfile->Write();
        }

        shapefile->cd();
        TCanvas *canv = new TCanvas();
        for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
        {
            numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
                numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                histo1D[histoName.c_str()]->Draw();
                string writename = "";
                if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
                {
                    writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__data_obs__nominal";
                }
                else
                {
                    writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__" + dataSetName +"__nominal";
                }
                //cout<<"writename  :"<<writename<<endl;
                histo1D[histoName.c_str()]->Write((writename).c_str());

                canv->SaveAs((pathPNG+histoName+".pdf").c_str());
            }
        }
    }


    treeLoader.UnLoadDataset();

    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error.root";
    //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
        temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 100);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
    }
};


void Split2SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlSys, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING SYS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr;   //add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath, numStr1, numStr2, histoName;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar1, splitVar2;
    string sbSplit1 = static_cast<ostringstream*>( &(ostringstream() << fbSplit1) )->str();
    string stSplit1 = static_cast<ostringstream*>( &(ostringstream() << ftSplit1) )->str();
    string swSplit1 = static_cast<ostringstream*>( &(ostringstream() << fwSplit1) )->str();
    string sbSplit2 = static_cast<ostringstream*>( &(ostringstream() << fbSplit2) )->str();
    string stSplit2 = static_cast<ostringstream*>( &(ostringstream() << ftSplit2) )->str();
    string swSplit2 = static_cast<ostringstream*>( &(ostringstream() << fwSplit2) )->str();

    errorfile->cd();
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<" Dataset:  :"<<dataSetName<<endl;
        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

        cout<<""<<endl;
        cout<<"file in use: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
        cout<<"initialised file"<<endl;
        string nTupleName = "Craneen__"+leptoAbbr;

         ttree[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
        nEntries = (int) ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        //sVarofinterest = "HT";
         ttree[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
         ttree[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
         ttree[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
         ttree[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);

        ////****************************************************Define plots***********************************************
        for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
        {
            numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
                numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();            
                histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                //cout<<"histoName: "<< histoName<<endl;
                histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
            }
        }
        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
             ttree[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi
            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }
            bool outsiderange = true;
            if(splitVar1 >= ftSplit1) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
            {
                if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
                    histoName = dataSetName + stSplit1 + sSplitVar1 + stSplit2 + sSplitVar2;
                    outsiderange = false;
                } 
                else //If it doesn't belong in the last bin, find out which it belongs in
                {
                    for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
                    {
                        if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
                        {
                            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                            histoName = dataSetName + stSplit1+ sSplitVar1 + numStr2 + sSplitVar2;
                            outsiderange = false;
                            break;
                        }
                    }
                } 
            }
            else //If it doesn't belong in the last bin for var1, find out which it belongs in
            {
                for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
                {
                    numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                    if(splitVar1>=s && splitVar1<(s+fwSplit1)) //splitVar falls inside one of the bins
                    {

                        if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
                            histoName = dataSetName + numStr1 + sSplitVar1 + stSplit2 + sSplitVar2;
                            outsiderange = false;
                        }    
                        else //If it doesn't belong in the last bin, find out which it belongs in
                        {
                            for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
                            {
                                if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
                                {
                                    numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                                    histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                                    outsiderange =false;
                                    break;
                                }
                            }
                        } 
                        //break;
                    }
                }
            }
            if (outsiderange == true){
                continue;
            }
            histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
        }

        ////****************************************************Fill plots***********************************************

        TCanvas *canv2 = new TCanvas();

        for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
        {
            numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
                shapefile->cd();
                numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
                histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
                if(Normalise)
                {
                    double dIntegral = histo1D[histoName.c_str()]->Integral();
                    histo1D[histoName.c_str()]->Scale(1./dIntegral);
                }
                histo1D[histoName.c_str()]->Draw();
                string writename = "";

                if(dataSetName == "TTScaleDown")
                {
                    writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__TTJets__scaleDown";
                    //cout<<"writename  :"<<writename<<endl;
                    histo1D[histoName.c_str()]->Write((writename).c_str());
                    canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                    errorfile->cd();
                    //errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
                    errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                    histo1D[histoName.c_str()]->Write("Minus");
                    //errorfile->Write();
                }

                if(dataSetName == "TTScaleUp")
                {
                    writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__TTJets__scaleUp";
                    //cout<<"writename  :"<<writename<<endl;
                    histo1D[histoName.c_str()]->Write((writename).c_str());
                    canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                    errorfile->cd();
                    errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                    histo1D[histoName.c_str()]->Write("Plus");
                    //errorfile->Write();
                }
            }
        }
    }
};

void Split2_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, bool jetTagsplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom){
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    const char *xmlfile = xmlNom.c_str();
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;
    int nDatasets = datasets.size();
    TH1F *tempHisto;
    float tempEntries;
    int nChannels = 0;

    string numStr1, numStr2, binname, histoName, dataSetName;
    ofstream card;
    card.open ("datacard.txt");

    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            //for(int i=0; i<datasets.size()-1; i++){
                nChannels += 1;
            //}
        }  
    } 

    card << "imax " + static_cast<ostringstream*>( &(ostringstream() << nChannels) )->str() + "\n";
    card << "jmax " + static_cast<ostringstream*>( &(ostringstream() << nDatasets-2) )->str() + "\n";
    card << "kmax *\n";
    card << "---------------\n";
    card << "shapes * * "+shapefileName+"  $CHANNEL__$PROCESS__nominal  $CHANNEL__$PROCESS__$SYSTEMATIC\n";
    card << "---------------\n";
    card << "bin                               ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
            binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            card << binname + " ";
        }  
    }    
    card << "\n";
    card << "observation                               ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
            binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            tempEntries = 0;
            for (int j = 0; j<nDatasets; j++){
                dataSetName=datasets[j]->Name();
                if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
                {
                    histoName = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__data_obs__nominal";
                    tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
                    tempEntries = tempHisto->GetSumOfWeights();  cout<<"tempEntries: "<<tempEntries<<"  "<<tempHisto->GetTitle()<<endl;
                    card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";                }
                else{
                    continue;
                }
            }
        }  
    } 
    card << "\n";

    card << "---------------------------\n";
    card << "bin                               ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
            binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            for(int i = 0; i<nDatasets-1; i++)    card << binname + " ";
        }  
    }
    card << "\n";
    card << "process                             ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            for (int j = 0; j<nDatasets; j++){
                dataSetName=datasets[j]->Name();
                if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
                {
                    continue;
                }
                else{
                    card<<dataSetName+"           ";
                }
            }
        }  
    }    

    card << "\n";
    card << "process                                ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            for(int i=0; i<datasets.size()-1; i++){
                card << static_cast<ostringstream*>( &(ostringstream() << i) )->str() + "                    ";
            }
        }  
    }  
    card << "\n";
    card << "rate                                ";
    for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
            numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
            binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
            tempEntries = 0;
            for (int j = 0; j<nDatasets; j++){
                dataSetName=datasets[j]->Name();
                if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
                {
                    continue;
                }
                else{
                    histoName = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__" + dataSetName +"__nominal";
                    tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
                    tempEntries = tempHisto->GetSumOfWeights();  cout<<"tempEntries: "<<tempEntries<<endl;
                    card << static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str() + "               ";
                }
            }
        }  
    }     
    card << "\n";
    card << "---------------------------\n";
    card << "lumi                  lnN           ";
    for (int i=0; i<nChannels*nDatasets; i++){
        card << "1.026                ";
    }
    card << "\n";
    for (int d = 0; d<nDatasets; d++){
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            continue;
        }
        else {
            if(dataSetName.find("tttt")!=string::npos){
                card << "tttt_norm              lnN          ";
                cout<<"chan "<<nChannels<<endl;
                for (int j = 0; j<nChannels; j++){
                    cout<<"J: "<<j<<endl;
                    card << "1.1                 -                    -                    -                  -                    ";
                }
                card << "\n";
            }
            else{
                card << dataSetName + "_norm      lnN           ";
                for(int k = 0; k<nChannels; k++){
                    for (int dash1 = 0; dash1<d-1; dash1++){
                        card << "-                  ";
                    }
                    card << "1.04                  ";                   
                    for (int dash2 = 5; dash2>d; dash2-- ){
                        card << "-                  ";
                    }
                }
                card<<"\n";
            }   
        }        
    }

    card << "scale                shape           ";
    for (int d = 0; d<nDatasets; d++){
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("TTJets")!=string::npos)
        {
            for(int k = 0; k<nChannels; k++){
                for (int dash1 = 0; dash1<d-1; dash1++){
                    card << "-                  ";
                }
                card << "1                      ";                   
                for (int dash2 = 5; dash2>d; dash2-- ){
                    card << "-                  ";
                }
            }
            card<<"\n";
        }
        else {
            continue;
        }        
    }    

    card.close();
};


std::string intToStr (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

