//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:
//all-had ->679 * .46 = 312.34
//semi-lep ->679 *.45 = 305.55
//di-lep-> 679* .09 = 61.11

//ttbar @ NNLO 8 TeV:
//all-had -> 245.8 * .46 = 113.068
//semi-lep-> 245.8 * .45 = 110.61
//di-lep ->  245.8 * .09 = 22.122

#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <iostream>
#include <map>
#include <cstdlib>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
//#include "TopTreeAnalysisBase/Selection/interface/FourTopSelectionTable.h"

#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = false;
bool debug = false;
float topness;



int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TProfile*> histoProfile;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
float ElectronRelIso(TRootElectron* el, float rho);

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 14)
    {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
	for (int n_arg=1; n_arg<argc; n_arg++)
	  {
	    std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl; 
	  }
        return 1;
    }

    //Placing arguments in properly typed variables for Dataset creation

    const string dName              = argv[1];
    const string dTitle             = argv[2];
    const int color                 = strtol(argv[4], NULL, 10);
    const int ls                    = strtol(argv[5], NULL, 10);
    const int lw                    = strtol(argv[6], NULL, 10);
    const float normf               = strtod(argv[7], NULL);
    const float EqLumi              = strtod(argv[8], NULL);
    const float xSect               = strtod(argv[9], NULL);
    const float PreselEff           = strtod(argv[10], NULL);
    string fileName                 = argv[11];
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);
    vector<string> vecfileNames;
    for(int args = 11; args < argc-2; args++)
    {
        vecfileNames.push_back(argv[args]);
    }




    cout << "---Dataset accepted from command line---" << endl;
    cout << "Dataset Name: " << dName << endl;
    cout << "Dataset Title: " << dTitle << endl;
    cout << "Dataset color: " << color << endl;
    cout << "Dataset ls: " << ls << endl;
    cout << "Dataset lw: " << lw << endl;
    cout << "Dataset normf: " << normf << endl;
    cout << "Dataset EqLumi: " << EqLumi << endl;
    cout << "Dataset xSect: " << xSect << endl;
    cout << "Dataset File Name: " << vecfileNames[0] << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "----------------------------------------" << endl;
//    cin.get();



    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");

    int passed = 0;
    int ndefs =0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;


    bool bx25 = false;

    clock_t start = clock();


    int doJESShift = 0; // 0: off 1: minus 2: plus
    cout << "doJESShift: " << doJESShift << endl;

    int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
    cout << "doJERShift: " << doJERShift << endl;


    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the Displaced Top search ! "           << endl;
    cout << "*************************************************************" << endl;


    string postfix = "_Run2_TopTree_Study_" + dName; // to relabel the names of the output file

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel when dilepton is also true)
    bool dilepton = true;
    bool Muon = true;
    bool Electron = true;

    if(Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(Muon && !Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuMu";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(!Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_ElEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else
    {
        cerr<<"Correct Di-lepton Channel not selected."<<endl;
        exit(1);
    }



    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
//    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = true;
    anaEnv.loadGenJetCollection = false;
    anaEnv.loadGenEventCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.loadTrackMETCollection = false;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;
    int verbose = 2;//anaEnv.Verbose;



    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi*normf);
    datasets.push_back(theDataset);
    float Luminosity = 5.59; //pb^-1??


    string dataSetName;

    // declare Scale factor related objects                                                                               
    LeptonTools ElectronSF = new LeptonTools(true); // true-> will enable some cout
    LeptonTools MuonsSF = new LeptonTools(true);  // true-> will enable some cout

    // load Sf
    ElectronSF.readElectronSF();
    double el_SfValue;
    //    string pathToCaliDir="/user/qpython/TopBrussels7X/CMSSW_7_4_12_patch4/src/TopBrussels/TopTreeAnalysisBase/Calibrations/LeptonSF/";

    //    double el_SfValue = ElectronSF.getElectronSF(1.2,30.0,"Nominal");



    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    dataSetName = theDataset->Name();
    if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
    {
        Luminosity = theDataset->EquivalentLumi();
        cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    //Output ROOT file
    string outputDirectory("MACRO_Output"+channelpostfix);
    mkdir(outputDirectory.c_str(),0777);
    string rootFileName (outputDirectory+"/FourTop"+postfix+channelpostfix+".root");
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_fatjets;
    vector < TRootMET* >      mets;

    //Global variable
    TRootEvent* event = 0;

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////

    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    //Muons
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]                                       = new MultiSamplePlot(datasets, "MuonEta", 40,-4, 4, "Muon #eta");
    MSPlot["MuonRelIsolation"]                              = new MultiSamplePlot(datasets, "MuonRelIsolation", 10, 0, .25, "RelIso");
    //Electrons
    MSPlot["ElectronRelIsolation"]                          = new MultiSamplePlot(datasets, "ElectronRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["ElectronPt"]                                    = new MultiSamplePlot(datasets, "ElectronPt", 30, 0, 300, "PT_{e}");
    MSPlot["ElectronEta"]                                   = new MultiSamplePlot(datasets, "ElectronEta", 40,-4, 4, "Jet #eta");
    MSPlot["NbOfElectronsPreSel"]                           = new MultiSamplePlot(datasets, "NbOfElectronsPreSel", 10, 0, 10, "Nb. of electrons");
    

    //Init Electron Plots
    MSPlot["InitElectronPt"]                                = new MultiSamplePlot(datasets, "InitElectronPt", 30, 0, 300, "PT_{e}");
    /*
    MSPlot["InitElectronEta"]                               = new MultiSamplePlot(datasets, "InitElectronEta", 40, -4, 4, "#eta");
    MSPlot["NbOfElectronsInit"]                             = new MultiSamplePlot(datasets, "NbOfElectronsInit", 10, 0, 10, "Nb. of electrons");
    MSPlot["InitElectronRelIsolation"]                      = new MultiSamplePlot(datasets, "InitElectronRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["InitElectronSuperClusterEta"]                   = new MultiSamplePlot(datasets, "InitElectronSuperClusterEta", 10, 0, 2.5, "#eta");
    MSPlot["InitElectrondEtaI"]                             = new MultiSamplePlot(datasets, "InitElectrondEtaI", 20, 0, .05, "#eta");
    MSPlot["InitElectrondPhiI"]                             = new MultiSamplePlot(datasets, "InitElectrondPhiI", 20, 0, .2, "#phi");
    MSPlot["InitElectronHoverE"]                            = new MultiSamplePlot(datasets, "InitElectronHoverE", 10, 0, .15, "H/E");
    MSPlot["InitElectrond0"]                                = new MultiSamplePlot(datasets, "InitElectrond0", 20, 0, .1, "d0");
    MSPlot["InitElectrondZ"]                                = new MultiSamplePlot(datasets, "InitElectrondZ", 10, 0, .25, "dZ");
    MSPlot["InitElectronEminusP"]                           = new MultiSamplePlot(datasets, "InitElectronEminusP", 10, 0, .25, "1/GeV");
    MSPlot["InitElectronConversion"]                        = new MultiSamplePlot(datasets, "InitElectronConversion", 2, 0, 2, "Conversion Pass");
    MSPlot["InitElectronMissingHits"]                       = new MultiSamplePlot(datasets, "InitElectronMissingHits", 10, 0, 10, "MissingHits");
    MSPlot["InitElectronCutFlow"]                           = new MultiSamplePlot(datasets, "InitElectronCutFlow", 12, 0, 12, "CutNumber");
    MSPlot["InitElectronDiagRelIso"]                        = new MultiSamplePlot(datasets, "InitElectronDiagRelIso", 100, 0, 1, "RelIso");
    MSPlot["InitElectronDiagChIso"]                         = new MultiSamplePlot(datasets, "InitElectronDiagChIso", 100, 0, 1, "RelIso");
    MSPlot["InitElectronDiagNIso"]                          = new MultiSamplePlot(datasets, "InitElectronDiagNIso", 100, 0, 1, "RelIso");
    MSPlot["InitElectronDiagPhIso"]                         = new MultiSamplePlot(datasets, "InitElectronDiagPhIso", 100, 0, 1, "RelIso");
    MSPlot["InitElectronDiagPUChIso"]                       = new MultiSamplePlot(datasets, "InitElectronDiagPUChIso", 100, 0, 1, "RelIso");

    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]                             = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 20, 0, 1, "CSV b-disc.");
    //Jets
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["3rdJetPt"]                                      = new MultiSamplePlot(datasets, "3rdJetPt", 30, 0, 300, "PT_{jet3}");
    MSPlot["4thJetPt"]                                      = new MultiSamplePlot(datasets, "4thJetPt", 30, 0, 300, "PT_{jet4}");
    MSPlot["5thJetPt"]                                      = new MultiSamplePlot(datasets, "5thJetPt", 30, 0, 300, "PT_{jet5}");
    MSPlot["6thJetPt"]                                      = new MultiSamplePlot(datasets, "6thJetPt", 30, 0, 300, "PT_{jet6}");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    MSPlot["HTExcess2M"]                                    = new MultiSamplePlot(datasets, "HTExcess2M", 30, 0, 1500, "HT_{Excess 2 M b-tags}");
    MSPlot["HTH"]                                           = new MultiSamplePlot(datasets, "HTH", 20, 0, 1, "HTH");
    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");
    MSPlot["METCutAcc"]                                     = new MultiSamplePlot(datasets, "METCutAcc", 30, 0, 150, "MET cut pre-jet selection");
    MSPlot["METCutRej"]                                     = new MultiSamplePlot(datasets, "METCutRej", 30, 0, 150, "MET cut pre-jet selection");
    //Topology Plots
    MSPlot["TotalSphericity"]                               = new MultiSamplePlot(datasets, "TotalSphericity", 20, 0, 1, "Sphericity_{all}");
    MSPlot["TotalCentrality"]                               = new MultiSamplePlot(datasets, "TotalCentrality", 20, 0, 1, "Centrality_{all}");
    MSPlot["DiLepSphericity"]                               = new MultiSamplePlot(datasets, "DiLepSphericity", 20, 0, 1, "Sphericity_{ll}");
    MSPlot["DiLepCentrality"]                               = new MultiSamplePlot(datasets, "DiLepCentrality", 20, 0, 1, "Centrality_{ll}");
    MSPlot["TopDiLepSphericity"]                            = new MultiSamplePlot(datasets, "TopDiLepSphericity", 20, 0, 1, "Sphericity_{tll}");
    MSPlot["TopDiLepCentrality"]                            = new MultiSamplePlot(datasets, "TopDiLepCentrality", 20, 0, 1, "Centrality_{tll}");

    //ZMass window plots
    MSPlot["ZMassWindowWidthAcc"]                           = new MultiSamplePlot(datasets, "ZMassWindowWidthAcc", 20, 0, 100, "Z mass window width");
    MSPlot["ZMassWindowWidthRej"]                           = new MultiSamplePlot(datasets, "ZMassWindowWidthRej", 20, 0, 100, "Z mass window width");
    MSPlot["DiLepMass"]                                     = new MultiSamplePlot(datasets, "DiLepMass", 30, 0, 150, "m_{ll}");

    //n-1 Cut Plots
    MSPlot["NMinusOne"]                                     = new MultiSamplePlot(datasets, "NMinusOne", 8, 0, 8, "CutNumber");
    */
    
    ///////////////////
    // 1D histograms //
    ///////////////////

    ///////////////////
    // 2D histograms //
    ///////////////////
    histo2D["HTLepSep"] = new TH2F("HTLepSep","dR_{ll}:HT",50,0,1000, 20, 0,4);

    //Plots
    string pathPNG = "MSPlots_FourTop"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);

    cout <<"Making directory :"<< pathPNG  <<endl;
    vector<string> CutsselecTable;
    if(Muon && Electron)
      {
	CutsselecTable.push_back(string("initial"));
	CutsselecTable.push_back(string("Event cleaning and Trigger"));
	CutsselecTable.push_back(string("Exactly 2 Loose Electron"));
	CutsselecTable.push_back(string("Z Mass Veto"));
	CutsselecTable.push_back(string("At least 4 Jets"));
	CutsselecTable.push_back(string("At least 1 CSVM Jet"));
	CutsselecTable.push_back(string("At least 2 CSVM Jets"));
	CutsselecTable.push_back(string("At Least 500 GeV HT"));
      }
    


    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;
        bool nlo = false;
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("bx50") != std::string::npos) bx25 = false;
        else bx25 = true;

        if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;

        if(bx25) cout << "Dataset with 25ns Bunch Spacing!" <<endl;
        else cout << "Dataset with 50ns Bunch Spacing!" <<endl;
        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;


        //////////////////////////////////////////////
        // Setup Date string and nTuple for output  //
        //////////////////////////////////////////////

        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );

        int year = now->tm_year + 1900;
        int month =  now->tm_mon + 1;
        int day = now->tm_mday;
        int hour = now->tm_hour;
        int min = now->tm_min;
        int sec = now->tm_sec;

        string year_str;
        string month_str;
        string day_str;
        string hour_str;
        string min_str;
        string sec_str;

        ostringstream convert;   // stream used for the conversion
        convert << year;      // insert the textual representation of 'Number' in the characters in the stream
        year_str = convert.str();
        convert.str("");
        convert.clear();
        convert << month;      // insert the textual representation of 'Number' in the characters in the stream
        month_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        day_str = convert.str();
        convert.str("");
        convert.clear();
        convert << hour;      // insert the textual representation of 'Number' in the characters in the stream
        hour_str = convert.str();
        convert.str("");
        convert.clear();
        convert << min;      // insert the textual representation of 'Number' in the characters in the stream
        min_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        sec_str = convert.str();
        convert.str("");
        convert.clear();


        string date_str = day_str + "_" + month_str + "_" + year_str;

        cout <<"DATE STRING   "<<date_str << endl;

        //string dataSetName = datasets[d]->Name();
        string channel_dir = "Craneens"+channelpostfix;
        string date_dir = channel_dir+"/Craneens" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);



        //     string Ntupname = "Craneens/Craneen_" + dataSetName +postfix + "_" + date_str+  ".root";

        string Ntupname = "Craneens"+channelpostfix+"/Craneens"+ date_str  +"/Craneen_" + dataSetName +postfix + ".root";
        string Ntuptitle = "Craneen_" + channelpostfix;

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");

        // TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"nJets:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2M:HTb:HTH:HTRat:topness:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");

	//        TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"BDT:nJets:nFatJets:nWTags:nTopTags:nLtags:nMtags:nTtags:3rdJetPt:4thJetPt:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2L:HTb:HTH:HTRat:topness:EventSph:EventCen:DiLepSph:DiLepCen:TopDiLepSph:TopDiLepCen:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");


	Int_t nElectronsPostCut;
        Int_t nElectrons;
        Double_t pt_electron[10];
        Double_t phi_electron[10];
        Double_t eta_electron[10];
        Double_t E_electron[10];
        Double_t d0_electron[10];
        Double_t chargedHadronIso_electron[10];
        Double_t neutralHadronIso_electron[10];
        Double_t photonIso_electron[10];
        Double_t pfIso_electron[10];
	Double_t sf_electron[10];
        Int_t charge_electron[10];

        Int_t nMuonsPostCut;
        Int_t nMuons;
        Double_t pt_muon[10];
        Double_t phi_muon[10];
        Double_t eta_muon[10];
        Double_t E_muon[10];
        Double_t d0_muon[10];
        Double_t chargedHadronIso_muon[10];
        Double_t neutralHadronIso_muon[10];
        Double_t photonIso_muon[10];
        Double_t pfIso_muon[10];
        Int_t charge_muon[10];



	// define the output tree                                                                                         
	fout->cd();
        TTree* myTree = new TTree("tree","tree");

	// electrons
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//                                                            
        myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("chargedHadronIsoelectron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        myTree->Branch("neutralHadronIsoelectron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        myTree->Branch("photonIsoelectron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        myTree->Branch("pfIsoelectron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
	myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
	

	// muons
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("chargedHadronIsomuon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        myTree->Branch("neutralHadronIsomuon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        myTree->Branch("photonIsomuon",photonIso_muon,"photonIso_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");


	/*
	// define a second tree that will be filled after the selection
        TTree* postCutTree = new TTree("tree","tree");

	// electrons
        postCutTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        postCutTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        postCutTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        postCutTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        postCutTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        postCutTree->Branch("chargedHadronIsoelectron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        postCutTree->Branch("neutralHadronIsoelectron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        postCutTree->Branch("photonIsoelectron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        postCutTree->Branch("pfIsoelectron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        postCutTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        postCutTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
	

	// muons
        postCutTree->Branch("nMuons",&nMuons, "nMuons/I");
        postCutTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        postCutTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        postCutTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        postCutTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        postCutTree->Branch("chargedHadronIsomuon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        postCutTree->Branch("neutralHadronIsomuon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        postCutTree->Branch("photonIsomuon",photonIso_muon,"photonIso_muon[nMuons]/D");
        postCutTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        postCutTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        postCutTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");

	*/




	
        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

	int itrigger = -1, previousRun = -1;

        int start = 0;
        cout << "teh bugz!" << endl;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events in total dataset = "<<  ending  <<endl;

        int event_start = startEvent;
        if (verbose > 1) cout << " - Loop over events " << endl;


        double currentfrac =0.;
        double end_d;
        if(endEvent > ending)
            end_d = ending;
        else
            end_d = endEvent;

	end_d = 10000; //artifical ending for debug
        int nEvents = end_d - event_start;
        cout <<"Will run over "<<  (end_d - event_start) << " events..."<<endl;
        cout <<"Starting event = = = = "<< event_start  << endl;
        if(end_d < startEvent)
        {
            cout << "Starting event larger than number of events.  Exiting." << endl;
            return 1;
        }

        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootSubstructureJet*>    selectedFatJets;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);



        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////

        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
	  // Set default value for evertything that goes in the Tree
	  //	  for (imuons = 0){ 
	  //	  }
	  nMuons = nElectrons = 0;
	  
	  
            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, init_fatjets,  mets, debug);  //load event

            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;

            if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
	    //            MSPlot["NbOfElectronsInit"]->Fill(init_electrons.size(), datasets[d], true, Luminosity*scaleFactor );
            for (Int_t initel =0; initel < init_electrons.size(); initel++ )
            {
                float initreliso = ElectronRelIso(init_electrons[initel], rho);
		
                MSPlot["InitElectronPt"]->Fill(init_electrons[initel]->Pt(), datasets[d], true, Luminosity*scaleFactor);

            }

            float weight_0 = event->weight0();
            if (debug)cout <<"Weight0: " << weight_0 <<endl;
            if(nlo)
            {
                if(weight_0 < 0.0)
                {
                    scaleFactor = -1.0;  //Taking into account negative weights in NLO Monte Carlo
                    negWeights++;
                }
            }


            string graphName;

            //////////////////
            //Loading Gen jets
            //////////////////

	    // Synch cut 1: pt, eta, veto, DR                                                                                                                           
            // Declare one bool per cut                                                                                                                                 
            ///----                                                                                                                                                     

            // cut on electrons                                                                                                                                         
            Bool_t passedPtEl = false;
            Bool_t passedEtaEl = false;
            Bool_t passedIdEl = false;
            Bool_t passedIsoEl = false;

            //cut on muons                                                                                                                                              
            Bool_t passedPtMu = false;
            Bool_t passedEtaMu = false;
            Bool_t passedIdMu = false;
            Bool_t passedIsoMu = false;

            // vetos                                                                                                                                                    
            Bool_t passedExtraElVeto =false;
            Bool_t passedExtraMuVeto =false;

            // opposite signed leptons                                                                                                                                  
            Bool_t passedElMuOS = false;

            // non overlapping leptons                                                                                                                                  
            Bool_t passedElMuNotOverlaping = false;

	    
	    

            std::string filterName = "";
            int currentRun = event->runId();
            if(previousRun != currentRun)
            {
                cout <<"What run? "<< currentRun<<endl;
                previousRun = currentRun;


                if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
                {
                    if (debug)cout <<"event loop 6a"<<endl;

                    // cout << " RUN " << event->runId() << endl;
                }

            } //end previousRun != currentRun

            ///////////////////////
            // JER smearing
            //////////////////////

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////


            // Declare selection instance
            Run2Selection selection(init_jets, init_fatjets, init_muons, init_electrons, mets);

            // Define object selection cuts

	    if (debug)cout<<"Getting Jets"<<endl;
	    selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()
	    if (debug)cout<<"Getting Loose Muons"<<endl;
	    //	    selectedMuons                                       = selection.GetSelectedDisplacedMuons();
	    selectedMuons = init_muons;
	    if (debug)cout<<"Getting Loose Electrons"<<endl;
	    //	    selectedElectrons                                   = selection.GetSelectedDisplacedElectrons(); // VBTF ID
	    selectedElectrons = init_electrons;


            vector<TRootJet*>      selectedLightJets;

            int JetCut =0;
            int nMu, nEl, nLooseIsoMu;



            //////////////////////////////////
            // Preselection Lepton Operations //
            //////////////////////////////////
	    





            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
	    //            selecTable.Fill(d,0,scaleFactor);



	    /*
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
              
            }
	    */

            weightCount += scaleFactor;
            eventCount++;

            //////////////
            // N-1 Plot //
            //////////////
	    // missing...

            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection

            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
////            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event
//
//          cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  "<< nMu <<"  "  <<nEl<<" "<< selectedJets.size()   <<"  " <<  nLtags   <<"  "<<JetCut  <<"  "<<Muon<<" "<<Electron<<endl;


            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections






            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////

	    //            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);




            //////////////////////
            // Muon Based Plots //
            //////////////////////
	    
	    nMuons=0;
            for (Int_t selmu =0; selmu < selectedMuons.size() && selmu < 10; selmu++ )
            {
	      pt_muon[nMuons]=selectedMuons[selmu]->Pt();
	      phi_muon[nMuons]=selectedMuons[selmu]->Phi();
	      eta_muon[nMuons]=selectedMuons[selmu]->Eta();
	      E_muon[nMuons]=selectedMuons[selmu]->E();
	      d0_muon[nMuons]=selectedMuons[selmu]->d0();
	      chargedHadronIso_muon[nMuons]=selectedMuons[selmu]->chargedHadronIso(4);
	      neutralHadronIso_muon[nMuons]=selectedMuons[selmu]->neutralHadronIso(4);
	      photonIso_muon[nMuons]=selectedMuons[selmu]->photonIso(4);
	      pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
	      charge_muon[nMuons]=selectedMuons[selmu]->charge();
	      nMuons++;

	    }

            //////////////////////////
            // Electron Based Plots //
            //////////////////////////
	    
	    nElectrons=0;
            for (Int_t selel =0; selel < selectedElectrons.size() && selel < 10; selel++ )
	    {
	      pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
	      phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
	      eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
	      E_electron[nElectrons]=selectedElectrons[selel]->E();
	      d0_electron[nElectrons]=selectedElectrons[selel]->d0();
	      chargedHadronIso_electron[nElectrons]=selectedElectrons[selel]->chargedHadronIso(3);
	      neutralHadronIso_electron[nElectrons]=selectedElectrons[selel]->neutralHadronIso(3);
	      photonIso_electron[nElectrons]=selectedElectrons[selel]->photonIso(3);
	      pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
	      charge_electron[nElectrons]=selectedElectrons[selel]->charge();
	      if (selectedElectrons[selel]->Pt() >= 35 ) {
                sf_electron[nElectrons]=ElectronSF.getElectronSF(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),"Nominal");
              }
	      nElectrons++;
            }


            //////////////////////
            // Jets Based Plots //
            //////////////////////

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
            }



            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();

            ///////////////////
            //MET Based Plots//
            ///////////////////

	    //            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);

	    /*

            if((dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0) && Muon && Electron)
            {
                cout <<"Data Event Passed Selection.  Run: "<< event->runId() << " LumiSection: " << event->lumiBlockId() << " Event: "<< event->eventId() <<endl;
                cout <<"Muon Eta: " << selectedMuons[0]->Eta() << " Muon Pt: " << selectedMuons[0]->Pt() << " Electron Eta: " << selectedElectrons[0]->Eta() << " Electron Pt: " << selectedElectrons[0]->Pt() << endl;
            }

            if((dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0) && Muon && !Electron)
            {
                cout <<"Data Event Passed Selection.  Run: "<< event->runId() << " LumiSection: " << event->lumiBlockId() << " Event: "<< event->eventId() << " HT: " << HT << " nMTags: " << nMtags <<endl;
                cout <<"Muon1 Eta: " << selectedMuons[0]->Eta() << " Muon1 Pt: " << selectedMuons[0]->Pt() << " Muon2 Eta: " << selectedMuons[1]->Eta() << " Muon2 Pt: " << selectedMuons[1]->Pt() << endl;
            }
            if((dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0) && !Muon && Electron)
            {
                cout <<"Data Event Passed Selection.  Run: "<< event->runId() << " LumiSection: " << event->lumiBlockId() << " Event: "<< event->eventId() << " HT: " << HT << " nMTags: " << nMtags <<endl;
                cout <<"Electron1 Eta: " << selectedElectrons[0]->Eta() << " Electron1 Pt: " << selectedElectrons[0]->Pt() << " Electron2 Eta: " << selectedElectrons[1]->Eta() << " Electron2 Pt: " << selectedElectrons[1]->Pt() << endl;
            }
	    */

            //////////////////
            //Filling Tree//
            //////////////////
	    //	    /*
	    if (nElectrons + nMuons >= 2){
	      passed++;
	      myTree->Fill();
	    }
	    //	    */
	    
	    /*
	    for(int iele=0; iele<selectedElectrons.size(); iele++){
              if (abs(selectedElectrons[iele]->Pt()) > 25){
                passedPtEl = true;
                if (abs(selectedElectrons[iele]->Eta()) < 2.5){
                  passedEtaEl = true;
                  for(int imuo=0; imuo<selectedMuons.size(); imuo++){
                    if (abs(selectedMuons[imuo]->Pt()) > 25){
                      passedPtMu=true;
                      if (abs(selectedMuons[imuo]->Eta()) < 2.5){
                        passedEtaMu=true;
			myTree->Fill();
                        if(postCut_electronsTLV.size() == 1) {
                          passedExtraElVeto=true;
                          if(postCut_muonsTLV.size() == 1){
                            passedExtraMuVeto=true;
                            if(selected_electrons[iele]->charge() * selected_muons[imuo]->charge() == -1){ // to do! make a new el/muon collection                      
                              passedElMuOS=true;
                              if(postCut_electronsTLV[iele].DeltaR(postCut_muonsTLV[imuo]) > 0.5){
                                passedElMuNotOverlaping=true;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
*/



        } //End Loop on Events

        myTree->Write();
        tupfile->Close();
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    eventlist.close();

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;

    //////////////////////
    // Selection tables //
    //////////////////////

    //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    selecTable.TableCalculator(  true, true, true, true, true);

    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTable.Write(  outputDirectory+"/FourTop"+postfix+"_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);

    fout->cd();
    cout <<" after cd .."<<endl;

    string pathPNGJetCombi = pathPNG + "JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);

//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin();
            it != MSPlot.end();
            it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Draw(name.c_str(), 0, false, false, false, 1);
        temp->Write(fout, name, false, pathPNG, "pdf");
    }


    TDirectory* th2dir = fout->mkdir("Histos2D");
    th2dir->cd();


    for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {

        TH2F *temp = it->second;
        temp->Write();
        //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
        //tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    }
    delete fout;
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}

float ElectronRelIso(TRootElectron* el, float rho)
{
    double EffectiveArea = 0.;
    // Updated to 2015 EA from https://indico.cern.ch/event/370494/contribution/2/attachments/736984/1011061/Rami_update_on_CB_ELE_ID_PHYS14PU20bx25.pdf
    if (fabs(el->superClusterEta()) >= 0.0   && fabs(el->superClusterEta()) < 0.8   ) EffectiveArea = 0.1013;
    if (fabs(el->superClusterEta()) >= 0.8   && fabs(el->superClusterEta()) < 1.3 ) EffectiveArea = 0.0988;
    if (fabs(el->superClusterEta()) >= 1.3 && fabs(el->superClusterEta()) < 2.0   ) EffectiveArea = 0.0572;
    if (fabs(el->superClusterEta()) >= 2.0   && fabs(el->superClusterEta()) < 2.2   ) EffectiveArea = 0.0842;
    if (fabs(el->superClusterEta()) >= 2.2   && fabs(el->superClusterEta()) < 2.5   ) EffectiveArea = 0.1530;
    if (fabs(el->superClusterEta()) >= 2.5) EffectiveArea = -9999;

    double isoCorr = 0;
    isoCorr = el->neutralHadronIso(3) + el->photonIso(3) - rho*EffectiveArea;
    float isolation = (el->chargedHadronIso(3) + (isoCorr > 0.0 ? isoCorr : 0.0))/(el->Pt());

    return isolation;
}
