//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:
//all-had ->679 * .46 = 312.34
//semi-lep ->679 *.45 = 305.55
//di-lep-> 679* .09 = 61.11

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


Bool_t debug = false;
Bool_t testTree = false;


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
//float ElectronRelIso(TRootElectron* el, float rho);

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO

  if (debug)
    {
      cout << "list of arguments are ..." << endl;
      for (int n_arg=1; n_arg<argc; n_arg++)
	{
	  std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl;
	}
    }
    


  if(argc < 14 )
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
    // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list) 
    const int JobNum                = strtol(argv[argc-3], NULL, 10);
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);


    // all the files are stored from arg 11 to argc-2
    vector<string> vecfileNames;
    for(int args = 11; args < argc-2; args++)
      {
	vecfileNames.push_back(argv[args]);
      }
    
    if (debug){
      cout << "The list of file to run over will be printed..." << endl;
      for ( int nfiles = 0; nfiles < vecfileNames.size(); nfiles++)
	{
	  cout << "file number " << nfiles << " is " << vecfileNames[nfiles] << endl;
	}
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
    cout << "JobNum: " << JobNum << endl;
    bool isData= false;
    if(dName.find("Data")!=string::npos || dName.find("data")!=string::npos || dName.find("DATA")!=string::npos){
      isData = true;
      cout << "running on data !!!!" << endl;
    }
     
    cout << "----------------------------------------" << endl;
   
//    cin.get();



    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");


    int passed = 0;
    int passed_pc = 0;
    int ndefs =0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;


    bool bx25 = false; // faco

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


    // The luminosity has to be adapted in function of the data!!
    //    float Luminosity =16.344; //pb^-1 
    float Luminosity =1000.0; //pb^-1 


    string dataSetName;




    ////////////////////////////////
    //  Event Scale Factor
    ////////////////////////////////

    string pathToCaliDir="/user/qpython/TopBrussels7X/CMSSW_7_4_14/src/TopBrussels/TopTreeAnalysisBase/Calibrations/";


    // Muon SF
    string muonFile= "Muon_SF_TopEA.root";
    MuonSFWeight *muonSFWeight_ = new MuonSFWeight (pathToCaliDir+"LeptonSF/"+muonFile,"SF_totErr", false, false); // (... , ... , debug, print warning)

    // Electron SF
    string electronFile= "Elec_SF_TopEA.root";
    ElectronSFWeight *electronSFWeight_ = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electronFile,"GlobalSF", false, false); // (... , ... , debug, print warning)

    // PU SF
    LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root","pileup","pileup");

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

    // add jobs number at the end of file if 
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();

    string rootFileName (outputDirectory+"/DisplacedTop"+postfix+channelpostfix+".root");
    if (strJobNum != "0")
      {
	cout << "strJobNum is " << strJobNum << endl;
        rootFileName = outputDirectory+"/DisplacedTop"+postfix+channelpostfix+"_"+strJobNum+".root";
      }
    
    
    

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
    /*
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
    */

    ///////////////////
    // 1D histograms //
    ///////////////////

    ///////////////////
    // 2D histograms //
    ///////////////////
    //    histo2D["HTLepSep"] = new TH2F("HTLepSep","dR_{ll}:HT",50,0,1000, 20, 0,4);

    //Plots
    string pathPNG = "MSPlots_FourTop"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);


    // define cuts here

    // electron
    float el_pt_cut =30.; // 42
    float el_eta_cut = 2.4;


    // muon
    float mu_pt_cut = 30.;
    float mu_eta_cut = 2.4;
    float mu_iso_cut = 0.12;

    // convert into string

    std::ostringstream el_pt_cut_strs, el_eta_cut_strs, mu_pt_cut_strs, mu_eta_cut_strs, mu_iso_cut_strs;
    std::string el_pt_cut_str, el_eta_cut_str, mu_pt_cut_str, mu_eta_cut_str, mu_iso_cut_str;
    el_pt_cut_strs << el_pt_cut;
    el_eta_cut_strs << el_eta_cut;
    mu_pt_cut_strs << mu_pt_cut;
    mu_eta_cut_strs << mu_eta_cut;
    mu_iso_cut_strs << mu_iso_cut;

    el_pt_cut_str = el_pt_cut_strs.str();
    el_eta_cut_str = el_eta_cut_strs.str();
    mu_pt_cut_str = mu_pt_cut_strs.str();
    mu_eta_cut_str = mu_eta_cut_strs.str();
    mu_iso_cut_str = mu_iso_cut_strs.str();

    
    // check
    cout <<"Making directory :"<< pathPNG  <<endl;
    vector<string> CutFlowPresel;
    
    CutFlowPresel.push_back(string("initial"));
    CutFlowPresel.push_back(string("ecal crack veto"));
    CutFlowPresel.push_back(string("at least one good electron: pt $<$ "+el_pt_cut_str+", eta $<$ "+el_eta_cut_str));
    CutFlowPresel.push_back(string("at least one good muon: pt $<$ "+mu_pt_cut_str+", eta $<$ "+mu_eta_cut_str+", iso $<$ "+mu_iso_cut_str));
    CutFlowPresel.push_back(string("extra electron veto"));
    CutFlowPresel.push_back(string("extra muon veto"));
    CutFlowPresel.push_back(string("electron blinding d0"));
    CutFlowPresel.push_back(string("muon blinding d0"));
    //    CutFlowPresel.push_back(string("OS leptons"));
    CutFlowPresel.push_back(string(""));
    //    CutFlowPresel.push_back(string("Non overlaping leptons"));
    CutFlowPresel.push_back(string(" "));
    


    SelectionTable CutFlowPreselTable(CutFlowPresel, datasets);
    CutFlowPreselTable.SetLuminosity(Luminosity);
    CutFlowPreselTable.SetPrecision(1);


    // cutflow for one good electron
    
    vector<string> CutFlow_oneEl;
    
    CutFlow_oneEl.push_back(string("initial"));
    CutFlow_oneEl.push_back(string("electron crack veto"));
    CutFlow_oneEl.push_back(string("at least one good electron: pt $<$ "+el_pt_cut_str+", eta $<$ "+el_eta_cut_str));

    SelectionTable CutFlow_oneElTable(CutFlow_oneEl, datasets);
    CutFlow_oneElTable.SetLuminosity(Luminosity);
    CutFlow_oneElTable.SetPrecision(1);
    

    // cutflow for one good muon
    
    vector<string> CutFlow_oneMu;
    
    CutFlow_oneMu.push_back(string("initial"));
    CutFlow_oneMu.push_back(string("at least one good muon: pt $<$ "+mu_pt_cut_str+", eta $<$ "+mu_eta_cut_str+", iso $<$ "+mu_iso_cut_str));

    SelectionTable CutFlow_oneMuTable(CutFlow_oneMu, datasets);
    CutFlow_oneMuTable.SetLuminosity(Luminosity);
    CutFlow_oneMuTable.SetPrecision(1);
    



    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {

        cout << "Load Dataset" << endl;
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

	cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
	

	// Lumi scale
	Bool_t applyLumiScale = false;
	double lumiScale = -99.;

	if (isData || !applyLumiScale) { 
	  cout << "Lumi scale is not applied!! " << endl << endl;
	  if (debug){
	    cout << "isData is " << isData << endl;
	    cout << "applyLumiScale is " << !applyLumiScale << endl;
	  }
	  lumiScale = 1. ;
	}
	else {
	  lumiScale = Luminosity*xSect/datasets[d]->NofEvtsToRunOver();
	  cout << "dataset has equilumi = " << datasets[d]->EquivalentLumi() << endl;
	  cout << "the weight to apply for each event of this data set is " << "Lumi * (xs/NSample) -->  " << Luminosity << " * (" << xSect << "/" << datasets[d]->NofEvtsToRunOver() << ") = " << Luminosity*xSect/datasets[d]->NofEvtsToRunOver()  <<  endl;
	}



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


	// variables for electrons
	Int_t nElectronsPostCut;
        Int_t nElectrons;
        Double_t pt_electron[10];
        Double_t phi_electron[10];
        Double_t eta_electron[10];
        Double_t eta_superCluster_electron[10];
        Double_t E_electron[10];
        Double_t d0_electron[10];
        Double_t d0BeamSpot_electron[10];
        Double_t chargedHadronIso_electron[10];
        Double_t neutralHadronIso_electron[10];
        Double_t photonIso_electron[10];
        Double_t pfIso_electron[10];
        Int_t charge_electron[10];


	//bo id related variables
	Double_t sigmaIEtaIEta_electron[10];
	Double_t deltaEtaIn_electron[10];
	Double_t deltaPhiIn_electron[10];
	Double_t hadronicOverEm_electron[10];
	Int_t missingHits_electron[10];
	Bool_t passConversion_electron[10];
	// eo id realted variables

	Bool_t isEBEEGap;
	Double_t sf_electron[10];	

	// variables for muons
        Int_t nMuonsPostCut;
        Int_t nMuons;
        Double_t pt_muon[10];
        Double_t phi_muon[10];
        Double_t eta_muon[10];
        Double_t E_muon[10];
        Double_t d0_muon[10];
        Double_t d0BeamSpot_muon[10];
        Double_t chargedHadronIso_muon[10];
        Double_t neutralHadronIso_muon[10];
        Double_t photonIso_muon[10];
        Double_t pfIso_muon[10];
	Double_t sf_muon[10];
        Int_t charge_muon[10];


	// event related variables
	Int_t run_num;
	Int_t evt_num;
	Int_t lumi_num;
	Int_t nvtx;
	Int_t npu;
	Double_t puSF;


	// bo MytreePreCut

	// variables for electrons
        Int_t nElectrons_pc;
        Double_t pt_electron_pc[10];
        Double_t phi_electron_pc[10];
        Double_t eta_electron_pc[10];
        Double_t eta_superCluster_electron_pc[10];
        Double_t E_electron_pc[10];
        Double_t d0_electron_pc[10];
        Double_t d0BeamSpot_electron_pc[10];
        Double_t chargedHadronIso_electron_pc[10];
        Double_t neutralHadronIso_electron_pc[10];
        Double_t photonIso_electron_pc[10];
        Double_t pfIso_electron_pc[10];
        Int_t charge_electron_pc[10];


	//bo id related variables
	Double_t sigmaIEtaIEta_electron_pc[10];
	Double_t deltaEtaIn_electron_pc[10];
	Double_t deltaPhiIn_electron_pc[10];
	Double_t hadronicOverEm_electron_pc[10];
	Int_t missingHits_electron_pc[10];
	Bool_t passConversion_electron_pc[10];
	// eo id realted variables

	Bool_t isEBEEGap_pc;
	Double_t sf_electron_pc[10];	

	// variables for muons
        Int_t nMuons_pc;
        Double_t pt_muon_pc[10];
        Double_t phi_muon_pc[10];
        Double_t eta_muon_pc[10];
        Double_t E_muon_pc[10];
        Double_t d0_muon_pc[10];
        Double_t d0BeamSpot_muon_pc[10];
        Double_t chargedHadronIso_muon_pc[10];
        Double_t neutralHadronIso_muon_pc[10];
        Double_t photonIso_muon_pc[10];
        Double_t pfIso_muon_pc[10];
	Double_t sf_muon_pc[10];
        Int_t charge_muon_pc[10];


	// eo MytreePreCut




	// event related variables
	/*
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
	*/

	Int_t NEvent = datasets[d]->NofEvtsToRunOver();
	//	Double_t xs [1];
       



	// Define the bookkeeping tree
	/*
	TTree *bookkeeping = new TTree("startevents","startevents");
	//	bookkeeping->Branch("xs",xs,"xs/D");
	bookkeeping->Branch("NEvent",NEvent,"NEvent/I");
	*/

	// define the output tree                                                                                         
	fout->cd();
        TTree* myTree = new TTree("tree","tree");



	// electrons
	//	myTree->Branch("templatevar",templatevar,"templatevar[nElectrons]/D");
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//                                                            
        myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        myTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        myTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        myTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
	myTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
	myTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
	myTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
	myTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
	myTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
	myTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
	myTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
	myTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
	myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
	

	// muons
	//        myTree->Branch("templatevar",templatevar,"templatevar[nMuons]/D");
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
	myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
	myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");


	// event related variables
	myTree->Branch("run_num",&run_num,"run_num/I");
	myTree->Branch("evt_num",&evt_num,"evt_num/I");
	myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
	myTree->Branch("nvtx",&nvtx,"nvtx/I");
	myTree->Branch("npu",&npu,"npu/I");
	myTree->Branch("puSF",&puSF,"puSF/D");	


	
	// Define a secondary tree that is filled before the whole list of cut
	TTree* myPreCutTree = new TTree("preCutTree","preCutTree");


	// electrons
	//	myPreCutTree->Branch("templatevar",templatevar,"templatevar[nElectrons_pc]/D");
        myPreCutTree->Branch("nElectrons_pc",&nElectrons_pc, "nElectrons_pc/I");//                                                            
        myPreCutTree->Branch("pt_electron_pc",pt_electron_pc,"pt_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("phi_electron_pc",phi_electron_pc,"phi_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("eta_electron_pc",eta_electron_pc,"eta_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("eta_superCluster_electron_pc",eta_superCluster_electron_pc,"eta_superCluster_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("E_electron_pc",E_electron_pc,"E_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("chargedHadronIso_electron_pc",chargedHadronIso_electron_pc,"chargedHadronIso_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("neutralHadronIso_electron_pc",neutralHadronIso_electron_pc,"neutralHadronIso_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("photonIso_electron_pc",photonIso_electron_pc,"photonIso_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("pfIso_electron_pc",pfIso_electron_pc,"pfIso_electron_pc[nElectrons_pc]/D");
        myPreCutTree->Branch("charge_electron_pc",charge_electron_pc,"charge_electron_pc[nElectrons_pc]/I");
        myPreCutTree->Branch("d0_electron_pc",d0_electron_pc,"d0_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("d0BeamSpot_electron_pc",d0BeamSpot_electron_pc,"d0BeamSpot_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("sigmaIEtaIEta_electron_pc",sigmaIEtaIEta_electron_pc,"sigmaIEtaIEta_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("deltaEtaIn_electron_pc",deltaEtaIn_electron_pc,"deltaEtaIn_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("deltaPhiIn_electron_pc",deltaPhiIn_electron_pc,"deltaPhiIn_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("hadronicOverEm_electron_pc",hadronicOverEm_electron_pc,"hadronicOverEm_electron_pc[nElectrons_pc]/D");
	myPreCutTree->Branch("missingHits_electron_pc",missingHits_electron_pc,"missingHits_electron_pc[nElectrons_pc]/I");
	myPreCutTree->Branch("passConversion_electron_pc",passConversion_electron_pc,"passConversion_electron_pc[nElectrons_pc]/O)");
	myPreCutTree->Branch("isEBEEGap_pc",isEBEEGap_pc,"isEBEEGap_pc[nElectrons]/O)");
	myPreCutTree->Branch("sf_electron_pc",sf_electron_pc,"sf_electron_pc[nElectrons_pc]/D");
	

	// muons
	//        myPreCutTree->Branch("templatevar",templatevar,"templatevar[nMuons_pc]/D");
        myPreCutTree->Branch("nMuons_pc",&nMuons_pc, "nMuons_pc/I");
        myPreCutTree->Branch("pt_muon_pc",pt_muon_pc,"pt_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("phi_muon_pc",phi_muon_pc,"phi_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("eta_muon_pc",eta_muon_pc,"eta_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("E_muon_pc",E_muon_pc,"E_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("chargedHadronIso_muon_pc",chargedHadronIso_muon_pc,"chargedHadronIso_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("neutralHadronIso_muon_pc",neutralHadronIso_muon_pc,"neutralHadronIso_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("photonIso_muon_pc",photonIso_muon_pc,"photonIso_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("pfIso_muon_pc",pfIso_muon_pc,"pfIso_muon_pc[nMuons_pc]/D");
        myPreCutTree->Branch("charge_muon_pc",charge_muon_pc,"charge_muon_pc[nMuons_pc]/I");
        myPreCutTree->Branch("d0_muon_pc",d0_muon_pc,"d0_muon_pc[nMuons_pc]/D");
	myPreCutTree->Branch("d0BeamSpot_muon_pc",d0BeamSpot_muon_pc,"d0BeamSpot_muon_pc[nMuons_pc]/D");
	myPreCutTree->Branch("sf_muon_pc",sf_muon_pc,"sf_muon_pc[nMuons_pc]/D");






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

	//	end_d = 10000; //artifical ending for debug
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
	  if (debug) cout << "just entered the event loop!" << endl;
	  // Set default value for evertything that goes in the Tree
	  //	  for (imuons = 0){ 
	  //	  }
	  nMuons = nElectrons = 0;
	  
	  
            double ievt_d = ievt;
	    currentfrac = ievt_d/end_d;
            if (debug)cout << endl << endl << "Starting a new event loop!"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, init_fatjets,  mets, debug);  //load event

            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;

            if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
	    //            MSPlot["NbOfElectronsInit"]->Fill(init_electrons.size(), datasets[d], true, Luminosity*globalScaleFactor );


	    /*
	    for (Int_t initel =0; initel < init_electrons.size(); initel++ )
            {
	      float initreliso = ElectronRelIso(init_electrons[initel], rho);
	      MSPlot["InitElectronPt"]->Fill(init_electrons[initel]->Pt(), datasets[d], true, Luminosity*globalScaleFactor);

            }
	    */



            string graphName;

            //////////////////
            //Loading Gen jets
            //////////////////

	    // Synch cut 1: pt, eta, veto, DR                                                                                                                           
            // Declare one bool per cut
	    /// ---

            // cut on electrons
	    Bool_t passedEcalCrackVeto = false;
	    Bool_t passedGoodEl =false;
	    Bool_t passedBlindingEl =false;
	    /*
            Bool_t passedPtEl = false;
            Bool_t passedEtaEl = false;
            Bool_t passedIdEl = false;
            Bool_t passedIsoEl = false;
	    */

            //cut on muons
	    Bool_t passedGoodMu = false;
	    Bool_t passedBlindingMu =false;
	    /*
            Bool_t passedPtMu = false;
            Bool_t passedEtaMu = false;
            Bool_t passedIdMu = false;
            Bool_t passedIsoMu = false;
	    */

            // vetos
            Bool_t passedExtraElVeto =false;
            Bool_t passedExtraMuVeto =false;

            // opposite signed leptons
            Bool_t passedElMuOS = false;

            // non overlapping leptons
            Bool_t passedElMuNotOverlaping = false;

	    /// ---
	    
	    
	    string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	    if(previousFilename != currentFilename)
	      {
                previousFilename = currentFilename;
                iFile++;
                cout<<"File changed!!! => "<<currentFilename<<endl;
	      }	    

	    run_num=event->runId();
	    evt_num=event->eventId();
	    lumi_num=event->lumiBlockId();
	    nvtx=vertex.size();
	    npu=(int)event->nTruePU();



	    ///////////////////////////////////////////
            //  Trigger
            ///////////////////////////////////////////
	    /*
	    bool trigged = false;
            std::string filterName = "";
            int currentRun = event->runId();
            if(previousRun != currentRun){

	      cout << "Changing run number. Run is  "<< currentRun<<endl;
	      previousRun = currentRun;
	      
	      if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos){
		if (debug)cout <<"applying trigger for data!"<<endl;
		itrigger = treeLoader.iTrigger (string ("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_v2"), currentRun, iFile);
		
		cout << " RUN " << event->runId() << endl;
	      }

            } //end previousRun != currentRun
	    */



            ///////////////////////
            // JER smearing
            //////////////////////

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////




	    // Apply trigger selection
	    /*
            trigged = treeLoader.EventTrigged (itrigger);
            //trigged = true;  // Disabling the HLT requirement
            if (trigged) cout << " trigged is true!!" <<endl;

            if(itrigger == 9999 ) cout << "Lumi Block: " << event->lumiBlockId() << " Event: " << event->eventId() << endl;
	    */

            // Declare selection instance
            Run2Selection selection(init_jets, init_fatjets, init_muons, init_electrons, mets);

            // Define object selection cuts

	    if (debug)cout<<"Getting Jets"<<endl;
	    selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()

	    // make a new collections of muons
	    if (debug)cout<<"Getting Muons"<<endl;
	    selectedMuons = selection.GetSelectedDisplacedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut); // pt, eta, iso
	    //selectedMuons = selection.GetSelectedMuons();
	    //selectedMuons = init_muons;

	    // make a new collections of electrons
	    if (debug)cout<<"Getting Electrons"<<endl;
	    selectedElectrons = selection.GetSelectedDisplacedElectrons(el_pt_cut, el_eta_cut);// pt, eta
	    //selectedElectrons = selection.GetSelectedElectrons();
	    //selectedElectrons = init_electrons;


	    //	    vector<TRootJet*>      selectedLightJets;

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
	    //            CutFlowTable.Fill(d,0,globalScaleFactor);



	    /*
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
              
            }
	    */

	    //            weightCount += globalScaleFactor;
            eventCount++;

            //////////////
            // N-1 Plot //
            //////////////
	    // missing...

            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection

	    //            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
////            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event
//
	    if (debug) cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  " << endl << init_muons.size() <<" "  <<init_electrons.size()<<" "<< selectedJets.size()   <<" "   <<" "<<JetCut  <<" "<<Muon<<" "<<Electron<<endl;


            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections






            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////

	    //            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*globalScaleFactor);




            //////////////////////////
            // Electron Based Plots //
            //////////////////////////
	    if (debug) cout << "before electrons loop" << endl;

	    
	    nElectrons=0;
            for (Int_t selel =0; selel < selectedElectrons.size() && selel < 10; selel++ )
	    {
	      pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
	      phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
	      eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
	      eta_superCluster_electron[nElectrons]=selectedElectrons[selel]->superClusterEta();
	      E_electron[nElectrons]=selectedElectrons[selel]->E();
	      d0_electron[nElectrons]=selectedElectrons[selel]->d0();
	      d0BeamSpot_electron[nElectrons]=selectedElectrons[selel]->d0BeamSpot();
	      chargedHadronIso_electron[nElectrons]=selectedElectrons[selel]->chargedHadronIso(3);
	      neutralHadronIso_electron[nElectrons]=selectedElectrons[selel]->neutralHadronIso(3);
	      photonIso_electron[nElectrons]=selectedElectrons[selel]->photonIso(3);
	      pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
	      charge_electron[nElectrons]=selectedElectrons[selel]->charge();
	      sigmaIEtaIEta_electron[nElectrons]=selectedElectrons[selel]->sigmaIEtaIEta();
	      deltaEtaIn_electron[nElectrons]=selectedElectrons[selel]->deltaEtaIn();
	      deltaPhiIn_electron[nElectrons]=selectedElectrons[selel]->deltaPhiIn();
	      hadronicOverEm_electron[nElectrons]=selectedElectrons[selel]->hadronicOverEm();
	      missingHits_electron[nElectrons]=selectedElectrons[selel]->missingHits();
	      passConversion_electron[nElectrons]=selectedElectrons[selel]->passConversion();
	      // following code found in http://cmslxr.fnal.gov/source/RecoEgamma/PhotonIdentification/src/PhotonIsolationCalculator.cc#0520
	      //	      isEBEEGap_pc = false;
	      isEBEEGap = false;
	      Double_t eta =  eta_superCluster_electron[nElectrons];
	      Double_t feta = fabs(eta);
	      if (fabs(feta-1.479)<.1) isEBEEGap = true ;

	      sf_electron[nElectrons]=electronSFWeight_->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
	      if (debug) cout << "in electrons loops, nelectrons equals to " << nElectrons << " and pt equals to " << pt_electron[nElectrons] << endl;
	      nElectrons++;
            }


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
	      d0BeamSpot_muon[nMuons]=selectedMuons[selmu]->d0BeamSpot();
	      chargedHadronIso_muon[nMuons]=selectedMuons[selmu]->chargedHadronIso(4);
	      neutralHadronIso_muon[nMuons]=selectedMuons[selmu]->neutralHadronIso(4);
	      photonIso_muon[nMuons]=selectedMuons[selmu]->photonIso(4);
	      pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
	      charge_muon[nMuons]=selectedMuons[selmu]->charge();
	      sf_muon[nMuons]=muonSFWeight_->at(selectedMuons[selmu]->Eta(),selectedMuons[selmu]->Pt(),0);
	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons << " and pt equals to " << pt_muon[nMuons] << endl;
	      nMuons++;


	    }



	    // Precut Tree
	    
            //////////////////////////
            // Electron Based Plots //
            //////////////////////////
	    if (debug) cout << "before electrons loop" << endl;

	    
	    nElectrons_pc=0;
            for (Int_t initel =0; initel < init_electrons.size() && initel < 10; initel++ )
	    {
	      pt_electron_pc[nElectrons_pc]=init_electrons[initel]->Pt();
	      phi_electron_pc[nElectrons_pc]=init_electrons[initel]->Phi();
	      eta_electron_pc[nElectrons_pc]=init_electrons[initel]->Eta();
	      eta_superCluster_electron_pc[nElectrons_pc]=init_electrons[initel]->superClusterEta();
	      E_electron_pc[nElectrons_pc]=init_electrons[initel]->E();
	      d0_electron_pc[nElectrons_pc]=init_electrons[initel]->d0();
	      d0BeamSpot_electron_pc[nElectrons_pc]=init_electrons[initel]->d0BeamSpot();
	      chargedHadronIso_electron_pc[nElectrons_pc]=init_electrons[initel]->chargedHadronIso(3);
	      neutralHadronIso_electron_pc[nElectrons_pc]=init_electrons[initel]->neutralHadronIso(3);
	      photonIso_electron_pc[nElectrons_pc]=init_electrons[initel]->photonIso(3);
	      pfIso_electron_pc[nElectrons_pc]=init_electrons[initel]->relPfIso(3,0);
	      charge_electron_pc[nElectrons_pc]=init_electrons[initel]->charge();
	      sigmaIEtaIEta_electron_pc[nElectrons_pc]=init_electrons[initel]->sigmaIEtaIEta();
	      deltaEtaIn_electron_pc[nElectrons_pc]=init_electrons[initel]->deltaEtaIn();
	      deltaPhiIn_electron_pc[nElectrons_pc]=init_electrons[initel]->deltaPhiIn();
	      hadronicOverEm_electron_pc[nElectrons_pc]=init_electrons[initel]->hadronicOverEm();
	      missingHits_electron_pc[nElectrons_pc]=init_electrons[initel]->missingHits();
	      passConversion_electron_pc[nElectrons_pc]=init_electrons[initel]->passConversion();


	      sf_electron_pc[nElectrons_pc]=electronSFWeight_->at(init_electrons[initel]->Eta(),init_electrons[initel]->Pt(),0);
	      if (debug) cout << "in electrons loops, nelectrons equals to " << nElectrons_pc << " and pt equals to " << pt_electron_pc[nElectrons_pc] << endl;
	      nElectrons_pc++;
            }


            //////////////////////
            // Muon Based Plots //
            //////////////////////
	    
	    nMuons_pc=0;
            for (Int_t initmu =0; initmu < init_muons.size() && initmu < 10; initmu++ )
            {
	      pt_muon_pc[nMuons_pc]=init_muons[initmu]->Pt();
	      phi_muon_pc[nMuons_pc]=init_muons[initmu]->Phi();
	      eta_muon_pc[nMuons_pc]=init_muons[initmu]->Eta();
	      E_muon_pc[nMuons_pc]=init_muons[initmu]->E();
	      d0_muon_pc[nMuons_pc]=init_muons[initmu]->d0();
	      d0BeamSpot_muon_pc[nMuons_pc]=init_muons[initmu]->d0BeamSpot();
	      chargedHadronIso_muon_pc[nMuons_pc]=init_muons[initmu]->chargedHadronIso(4);
	      neutralHadronIso_muon_pc[nMuons_pc]=init_muons[initmu]->neutralHadronIso(4);
	      photonIso_muon_pc[nMuons_pc]=init_muons[initmu]->photonIso(4);
	      pfIso_muon_pc[nMuons_pc]=init_muons[initmu]->relPfIso(4,0);
	      charge_muon_pc[nMuons_pc]=init_muons[initmu]->charge();
	      sf_muon_pc[nMuons_pc]=muonSFWeight_->at(init_muons[initmu]->Eta(),init_muons[initmu]->Pt(),0);
	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons_pc << " and pt equals to " << pt_muon_pc[nMuons_pc] << endl;
	      nMuons_pc++;


	    }


            //////////////////////
            // Jets Based Plots //
            //////////////////////
	    if (debug) cout << "before jets loop" << endl;

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
            }



            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();

            ///////////////////
            //MET Based Plots//
            ///////////////////

	    //            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*globalScaleFactor);

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
	    if (debug) cout << "filling the tree, sum of leptons equals to " << nElectrons + nMuons << endl;
	    if (nElectrons_pc >= 1 && nMuons_pc >= 1){
		myPreCutTree->Fill(); 
		passed_pc++;
	    }

	    if (testTree && (nElectrons >= 1 || nMuons >= 1 )){
	      myTree->Fill();
	      passed++;
	    }
	    
	    if (debug) cout << " DONE filling the tree, sum of leptons equals to " <<nElectrons + nMuons << endl;


	    // filling the cutflow

	    // preselection cut
	    CutFlowPreselTable.Fill(d,0,globalScaleFactor*lumiScale);
	    if(!isEBEEGap){
	      passedEcalCrackVeto = true;
	      CutFlowPreselTable.Fill(d,1,globalScaleFactor*lumiScale);
	      if (selectedElectrons.size() >= 1 ){
		passedGoodEl=true;
		CutFlowPreselTable.Fill(d,2,globalScaleFactor*lumiScale);
		if (selectedMuons.size() >= 1 ){
		  passedGoodMu=true;
		  CutFlowPreselTable.Fill(d,3,globalScaleFactor*lumiScale);
		  if (selectedElectrons.size() == 1 ){
		    passedExtraElVeto = true;
		    CutFlowPreselTable.Fill(d,4,globalScaleFactor*lumiScale);
		    if (selectedMuons.size() == 1 ){
		      passedExtraMuVeto = true;
		      CutFlowPreselTable.Fill(d,5,globalScaleFactor*lumiScale);
		      if(abs(selectedElectrons[0]->d0BeamSpot()) < 0.01){
			passedBlindingEl = true;
			CutFlowPreselTable.Fill(d,6,globalScaleFactor*lumiScale);
			if (abs(selectedMuons[0]->d0BeamSpot()) < 0.01){
			  passedBlindingMu=true;
			  CutFlowPreselTable.Fill(d,7,globalScaleFactor*lumiScale);
			  //			  if(selectedElectrons[0]->charge() * selectedMuons[0]->charge() == -1){
			  if (1){
			    passedElMuOS=true;
			    CutFlowPreselTable.Fill(d,8,globalScaleFactor*lumiScale);
			    Double_t DeltaR = sqrt (2);// to be done
			    if (1){
			      //			  if(selectedElectrons[0]->DeltaR(selectedMuons[0]) > 0.5){
			      passedElMuNotOverlaping=true;
			      CutFlowPreselTable.Fill(d,9,globalScaleFactor*lumiScale);
			      passed++;
			      if (debug) cout << "About to fill the tree!! The number of event that have passed all the cuts is " << passed << endl;
			      myTree->Fill(); 
			      //			      cout << "Filling Tree!!!" << endl;
			      if (debug) cout << "puSF is " << puSF << endl;
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }


	    // single el 
	    CutFlow_oneElTable.Fill(d,0,globalScaleFactor*lumiScale);
	    if(!isEBEEGap){
	      passedEcalCrackVeto = true;
	      CutFlow_oneElTable.Fill(d,1,globalScaleFactor*lumiScale);
	      if (selectedElectrons.size() >= 1 ){
		passedGoodEl=true;
		CutFlow_oneElTable.Fill(d,2,globalScaleFactor*lumiScale);
	      }
	    }


	    // single el 
	    CutFlow_oneMuTable.Fill(d,0,globalScaleFactor*lumiScale);
	    if (selectedMuons.size() >= 1 ){
	      passedGoodMu=true;
	      CutFlow_oneMuTable.Fill(d,1,globalScaleFactor*lumiScale);
	    }


	    
	    /*


	      
	    for(int iele=0; iele<selectedElectrons.size(); iele++){
              if (abs(selectedElectrons[iele]->Pt()) > 40){
                passedPtEl = true;
                if (abs(selectedElectrons[iele]->Eta()) < 2.5){
                  passedEtaEl = true;
                  for(int imuo=0; imuo<selectedMuons.size(); imuo++){
                    if (abs(selectedMuons[imuo]->Pt()) > 40){
                      passedPtMu=true;
                      if (abs(selectedMuons[imuo]->Eta()) < 2.5){
                        passedEtaMu=true;
		
			/*
                        if(postCut_electronsTLV.size() == 1) {
                          passedExtraElVeto=true;
                          if(postCut_muonsTLV.size() == 1){
                            passedExtraMuVeto=true;

                            if(selectedElectrons[iele]->charge() * selectedMuons[imuo]->charge() == -1){ // to do! make a new el/muon collection                      
                              passedElMuOS=true;
			      myTree->Fill();
			      passed++;
			      /*
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

	    //cout << "ievt is " << ievt << " and end_d is " << end_d << endl;            
	    //            if( run_num > 10000){//data
	    //	      isdata=1;
	    //            }
            

	    

        } //End Loop on Events

	//	bookkeeping->Fill();	    
	//	bookkeeping->Write();
	

	myTree->Write();
	myPreCutTree->Write();
	if (debug) cout << "Done writing the Tree" << endl;
        tupfile->Close();
        cout <<"n events before all the cuts is  =  "<< passed_pc <<endl;
        cout <<"n events after all the cuts is  =  "<< passed <<endl;
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
    CutFlowPreselTable.TableCalculator(  true, true, true, true, true);
    //    CutFlow_oneElTable.TableCalculator(  true, true, true, true, true);
    //    CutFlow_oneMuTable.TableCalculator(  true, true, true, true, true);

    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    if (strJobNum != "0")
      {
	CutFlowPreselTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Presel_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,true,false,false,true);
	//	CutFlow_oneElTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneEl_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,true,false,false,true);
	//	CutFlow_oneMuTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneMu_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,true,false,false,true);
      }
    else 
      {
	CutFlowPreselTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Presel_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);
	//	CutFlow_oneElTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneEl_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);
	//	CutFlow_oneMuTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneMu_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);
      }

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




    delete fout;
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}


/*
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
*/
