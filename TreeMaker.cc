

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
#include "TopTreeAnalysisBase/Selection/interface/MuonSelection.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronSelection.h"
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
#include "TopTreeAnalysisBase/Tools/interface/Trigger.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;


Bool_t debug = false;
Bool_t testTree = false;
//if (debug) testTree=true;


int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TProfile*> histoProfile;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;


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
  const string channel            = argv[argc-4];
  const int JobNum                = strtol(argv[argc-3], NULL, 10);
  const int startEvent            = strtol(argv[argc-2], NULL, 10);
  const int endEvent              = strtol(argv[argc-1], NULL, 10);


  // all the files are stored from arg 11 to argc-4
  vector<string> vecfileNames;
  for(int args = 11; args < argc-4; args++) 
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
  cout << "Channel is " << channel << endl;
  cout << "Beginning Event: " << startEvent << endl;
  cout << "Ending Event: " << endEvent << endl;
  cout << "JobNum: " << JobNum << endl;
  bool isData= false;
  if(dName.find("Data")!=string::npos || dName.find("data")!=string::npos || dName.find("DATA")!=string::npos){
    isData = true;
    cout << "running on data !!!!" << endl;
  }
     
  cout << "----------------------------------------" << endl;
   

  //  ofstream eventlist;
  //  eventlist.open ("interesting_events_mu.txt");


  int passed = 0;
  int passed_pc = 0;
  int passed_elel =0;
  int passed_mumu =0;
  int ndefs =0;
  int negWeights = 0;
  float weightCount = 0.0;
  int eventCount = 0;
  int skippedEvent =0;


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

  bool printTriggers = false;
  bool applyTriggers = true;
  string channelpostfix = "";
  string xmlFileName = "";

  //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel when dilepton is also true)
  bool elel = false; 
  bool elmu = false; 
  bool mumu = false; 



  if(channel=="ElMu" || channel== "MuEl")
    {
      cout << " --> Using the Muon-Electron channel..." << endl;
      elmu=true;
      channelpostfix = "_MuEl";
    }
  else if(channel=="MuMu")
    {
      cout << " --> Using the Muon-Muon channel..." << endl;
      mumu=true;
      channelpostfix = "_MuMu";
    }
  else if(channel=="ElEl")
    {
      cout << " --> Using the Electron-Electron channel..." << endl;
      elel=true;
      channelpostfix = "_ElEl";
    }
  else
    {
      cerr << "The channel --" << channel << "-- is not in the list of authorised channels !!"<<endl;
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
  //  anaEnv.TrackMETCollection = "";
  //  anaEnv.GenEventCollection = "GenEvent";
  anaEnv.NPGenEventCollection = "NPGenEvent";
  anaEnv.MCParticlesCollection = "MCParticles";
  anaEnv.loadFatJetCollection = true;
  //  anaEnv.loadGenJetCollection = false;
  //  anaEnv.loadGenEventCollection = false;
  anaEnv.loadNPGenEventCollection = false;
  anaEnv.loadMCParticles = true;
  //  anaEnv.loadTrackMETCollection = false;
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
  //    theDataset->SetEquivalentLuminosity(EqLumi*normf);
  datasets.push_back(theDataset);


  string dataSetName;

  ////////////////////////////////
  //  Event Scale Factor
  ////////////////////////////////

  string pathToCaliDir="../TopTreeAnalysisBase/Calibrations/";


  /// Leptons

  // Muon SF
  double muonSFID, muonSFIso;
  //    MuonSFWeight *muonSFWeight_ = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"Muon_SF_TopEA.root","SF_totErr", true, true); // (... , ... , debug, print warning)
  MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, false, false);
  //    MuonSFWeight *muonSFWeightID_M = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
  //  MuonSFWeight *muonSFWeightID_L = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
  
  MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, false, false);  // Tight RelIso, Tight ID
  //   MuonSFWeight *muonSFWeightIso_TM = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Tight RelIso, Medium ID
  //   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Tight ID
  //   MuonSFWeight *muonSFWeightIso_LM = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Medium ID
  //   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Loose ID


  // Electron SF
  string electronFile= "Elec_SF_TopEA.root";
  ElectronSFWeight *electronSFWeight_ = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electronFile,"GlobalSF",true, false, false); // (... , ... , debug, print warning)



  

  // PU SF
  LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert_Silver.root","pileup","pileup");

  /////////////////////////////////
  //  Loop over Datasets
  /////////////////////////////////

  dataSetName = theDataset->Name();
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

  map <string,TH1F*> histo1D;
    
  std::string  titlePlot = ""; 
  titlePlot = "cutFlow"+channelpostfix; 
  histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 13,-0.5,12.5);



  ///////////////////
  // 2D histograms //
  ///////////////////
  //    histo2D["HTLepSep"] = new TH2F("HTLepSep","dR_{ll}:HT",50,0,1000, 20, 0,4);

  //Plots
  string pathPNG = "MSPlots_FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/";
  //pathPNG = pathPNG +"/";
  mkdir(pathPNG.c_str(),0777);



  // -------------------------
  // bo defining cuts value --
  // -------------------------
  // think to do it in a loop as it will be faster and much less error prompt

  // electron
  float el_pt_cut = 42.; // 42
  float el_eta_cut = 2.4; // 2.4
  float el_d0_cut = 0.02; // 0.02


  // muon
  float mu_pt_cut = 40.; // 40
  float mu_eta_cut = 2.4; // 2.4
  float mu_iso_cut = 0.15; // 0.15
  float mu_d0_cut = 0.02; //0.02



  // change value depending on the trigger
  if (mumu){
    mu_pt_cut=35;
  }
  else if (elel){
    el_pt_cut=40;
  }
  else if (elmu){
    mu_pt_cut=40;
    el_pt_cut=40;
  }
  
  
  


  // convert into string

  std::ostringstream el_pt_cut_strs, el_eta_cut_strs, el_d0_cut_strs, mu_pt_cut_strs, mu_eta_cut_strs, mu_iso_cut_strs, mu_d0_cut_strs; 
  std::string el_pt_cut_str, el_eta_cut_str, el_d0_cut_str, mu_pt_cut_str, mu_eta_cut_str, mu_iso_cut_str, mu_d0_cut_str; 
  el_pt_cut_strs << el_pt_cut;
  el_eta_cut_strs << el_eta_cut;
  el_d0_cut_strs << el_d0_cut;
  mu_pt_cut_strs << mu_pt_cut;
  mu_eta_cut_strs << mu_eta_cut;
  mu_iso_cut_strs << mu_iso_cut;
  mu_d0_cut_strs << mu_d0_cut;
    
  el_pt_cut_str = el_pt_cut_strs.str();
  el_eta_cut_str = el_eta_cut_strs.str();
  el_d0_cut_str = el_d0_cut_strs.str();
  mu_pt_cut_str = mu_pt_cut_strs.str();
  mu_eta_cut_str = mu_eta_cut_strs.str();
  mu_iso_cut_str = mu_iso_cut_strs.str();
  mu_d0_cut_str = mu_d0_cut_strs.str();


  // -------------------------
  // eo defining cuts value --
  // -------------------------

    
  // check
  cout <<"Making directory :"<< pathPNG  <<endl;
  vector<string> CutFlowPresel;
    
  CutFlowPresel.push_back(string("initial"));
  CutFlowPresel.push_back(string("ecal crack veto"));
  CutFlowPresel.push_back(string("at least one good electron: pt $>$ "+el_pt_cut_str+", eta $<$ "+el_eta_cut_str));
  CutFlowPresel.push_back(string("at least one good muon: pt $>$ "+mu_pt_cut_str+", eta $<$ "+mu_eta_cut_str+", iso $<$ "+mu_iso_cut_str));
  CutFlowPresel.push_back(string("extra electron veto"));
  CutFlowPresel.push_back(string("extra muon veto"));
  CutFlowPresel.push_back(string("electron d0 $<$ "+el_d0_cut_str));
  CutFlowPresel.push_back(string("muon d0 $<$ "+mu_d0_cut_str));
  //    CutFlowPresel.push_back(string("OS leptons"));
  CutFlowPresel.push_back(string(""));
  //    CutFlowPresel.push_back(string("Non overlaping leptons"));
  CutFlowPresel.push_back(string(" "));
    


  SelectionTable CutFlowPreselTable(CutFlowPresel, datasets);
  //    CutFlowPreselTable.SetLuminosity(Luminosity);
  CutFlowPreselTable.SetPrecision(1);

  // ---------------------
  // bo Synch cut flows --
  // ---------------------



  // start a new table (electrons only)
  vector<string> CutFlow_oneEl;
  CutFlow_oneEl.push_back(string("initial"));
  CutFlow_oneEl.push_back(string("At least one electron with: pt $>$ "+el_pt_cut_str));
  CutFlow_oneEl.push_back(string("ecal crack veto (!isEBEEGAP)"));
  CutFlow_oneEl.push_back(string("electron with abs eta $<$ "+el_eta_cut_str));
  CutFlow_oneEl.push_back(string("electron Id"));
  CutFlow_oneEl.push_back(string("delta rho isolation"));
  CutFlow_oneEl.push_back(string("electron d0 $<$ "+el_d0_cut_str));


  SelectionTable CutFlow_oneElTable(CutFlow_oneEl, datasets);
  //    CutFlow_oneElTable.SetLuminosity(Luminosity);
  CutFlow_oneElTable.SetPrecision(1);
    

  // start a new table (muon only)
  vector<string> CutFlow_oneMu;
  CutFlow_oneMu.push_back(string("initial"));
  CutFlow_oneMu.push_back(string("at least one muon with pt $>$ "+mu_pt_cut_str));
  CutFlow_oneMu.push_back(string("muon with abs eta  $<$ "+mu_eta_cut_str));
  CutFlow_oneMu.push_back(string("muon id")); 
  CutFlow_oneMu.push_back(string("delta beta iso $<$ "+mu_iso_cut_str));
  CutFlow_oneMu.push_back(string("muon d0 $<$ "+mu_d0_cut_str));
    
  SelectionTable CutFlow_oneMuTable(CutFlow_oneMu, datasets);
  //    CutFlow_oneMuTable.SetLuminosity(Luminosity);
  CutFlow_oneMuTable.SetPrecision(1);


  // start a table (full synch)
  vector<string> CutFlow;
  CutFlow.push_back(string("initial"));
  CutFlow.push_back(string("At least one electron with: pt $>$ "+el_pt_cut_str));
  //    CutFlow.push_back(string("ecal crack veto (!isEBEEGAP)"));

  CutFlow.push_back(string("electron with abs eta $<$ "+el_eta_cut_str));
  CutFlow.push_back(string("electron Id"));
  CutFlow.push_back(string("delta rho isolation"));

  CutFlow.push_back(string("at least one muon with pt $>$ "+mu_pt_cut_str));
  CutFlow.push_back(string("muon with abs eta  $<$ "+mu_eta_cut_str));
  CutFlow.push_back(string("muon Id"));

  CutFlow.push_back(string("delta beta iso $<$ "+mu_iso_cut_str));
  CutFlow.push_back(string("extra electron veto"));
  CutFlow.push_back(string("extra muon veto"));

  CutFlow.push_back(string("electron d0 $<$ "+el_d0_cut_str));
  CutFlow.push_back(string("muon d0 $<$ "+mu_d0_cut_str));
  CutFlow.push_back(string("electron and muon with OS"));

  CutFlow.push_back(string("electron-muon pair with deltaR $>$ 0.5"));
    
  SelectionTable CutFlowTable(CutFlow, datasets);
  //    CutFlow_oneMuTable.SetLuminosity(Luminosity);
  CutFlowTable.SetPrecision(1);


  // ---------------------
  // eo Synch cut flows --
  // ---------------------


  ////////////////////////////
  ///  Initialise trigger  ///
  ////////////////////////////

  //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);

  Trigger * trigger;
  
  if (channel=="ElMu"){
    trigger = new Trigger(1, 1, 0, 1);
  }
  else if (channel=="MuMu"){
    trigger = new Trigger(1, 0 , 0, 1);
  }
  else if (channel=="ElEl"){
    trigger = new Trigger(0, 1 , 0, 1);
  }
  else cout << "Wrong chanel name" << endl;
  
  


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
	

      /// book triggers
      if (applyTriggers) {
	trigger->bookTriggers(isData);
      }
      
      // Lumi scale
      Bool_t applyLumiScale = false;
      double lumiScale = -99.;


      /*
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
      */



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
      Bool_t isId_electron[10];
      Bool_t isIso_electron[10];
      // eo id realted variables

      Bool_t isEBEEGap[10]; // change at some point.
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
      //        Double_t relIso_muon[10];
      Bool_t isId_muon[10];
      Bool_t isIso_muon[10];
      Double_t pfIso_muon[10];
      Double_t sf_muon[10];
      Int_t charge_muon[10];


      // event related variables
      Int_t run_num;
      Int_t event_num;
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
      Bool_t isId_electron_pc[10];
      Bool_t isIso_electron_pc[10];
      Bool_t isEBEEGap_pc[10];
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
      Bool_t isId_muon_pc[10];
      Bool_t isIso_muon_pc[10];

      // event related variables
      Int_t run_num_pc;
      Int_t event_num_pc;
      Int_t lumi_num_pc;
      Int_t nvtx_pc;
      Int_t npu_pc;
      Double_t puSF_pc;


      // eo MytreePreCut


      // bo MyDoubleElTree

      // variables for electrons
      Int_t nElectrons_elel;
      Double_t pt_electron_elel[10];
      Double_t phi_electron_elel[10];
      Double_t eta_electron_elel[10];
      Double_t eta_superCluster_electron_elel[10];
      Double_t E_electron_elel[10];
      Double_t vz_electron_elel[10]; 
      Double_t v0_electron_elel[10]; 
      Double_t d0_electron_elel[10];
      Double_t d0BeamSpot_electron_elel[10];
      Double_t chargedHadronIso_electron_elel[10];
      Double_t neutralHadronIso_electron_elel[10];
      Double_t photonIso_electron_elel[10];
      Double_t pfIso_electron_elel[10];
      Int_t charge_electron_elel[10];

      //bo id related variables
      Double_t sigmaIEtaIEta_electron_elel[10];
      Double_t deltaEtaIn_electron_elel[10];
      Double_t deltaPhiIn_electron_elel[10];
      Double_t hadronicOverEm_electron_elel[10];
      Int_t missingHits_electron_elel[10];
      Bool_t passConversion_electron_elel[10];
      // eo id realted variables
      Bool_t isId_electron_elel[10];
      Bool_t isIso_electron_elel[10];
      Bool_t isEBEEGap_elel[10];
      Double_t sf_electron_elel[10];	


      // variables for electronPairs
      Int_t nElectronPairs_elel; // if there is n electrons there is (n*n - n)/2 distinct pairs
      Double_t deltaVz_elel[10]; // max 5 electrons -> max (25 -5)/2 = 10 electronPairs
      Double_t deltaV0_elel[10];
      Double_t invMass_elel[10];

      // variables for muons
      Int_t nMuons_elel;
      Double_t pt_muon_elel[10];
      Double_t phi_muon_elel[10];
      Double_t eta_muon_elel[10];
      Double_t E_muon_elel[10];
      Double_t d0_muon_elel[10];
      Double_t d0BeamSpot_muon_elel[10];
      Double_t chargedHadronIso_muon_elel[10];
      Double_t neutralHadronIso_muon_elel[10];
      Double_t photonIso_muon_elel[10];
      Double_t pfIso_muon_elel[10];
      Double_t sf_muon_elel[10];
      Int_t charge_muon_elel[10];
      Bool_t isId_muon_elel[10];
      Bool_t isIso_muon_elel[10];


      // event related variables
      Int_t run_num_elel;
      Int_t event_num_elel;
      Int_t lumi_num_elel;
      Int_t nvtx_elel;
      Int_t npu_elel;
      Double_t evt_puSF_elel;
      Double_t evt_met_elel; //faco here


      // eo MyDoubleElTree


      // bo MyDoubleMuTree

      // variables for electrons
      Int_t nElectrons_mumu;
      Double_t pt_electron_mumu[10];
      Double_t phi_electron_mumu[10];
      Double_t eta_electron_mumu[10];
      Double_t eta_superCluster_electron_mumu[10];
      Double_t E_electron_mumu[10];
      Double_t d0_electron_mumu[10];
      Double_t d0BeamSpot_electron_mumu[10];
      Double_t chargedHadronIso_electron_mumu[10];
      Double_t neutralHadronIso_electron_mumu[10];
      Double_t photonIso_electron_mumu[10];
      Double_t pfIso_electron_mumu[10];
      Int_t charge_electron_mumu[10];


      //bo id related variables
      Double_t sigmaIEtaIEta_electron_mumu[10];
      Double_t deltaEtaIn_electron_mumu[10];
      Double_t deltaPhiIn_electron_mumu[10];
      Double_t hadronicOverEm_electron_mumu[10];
      Int_t missingHits_electron_mumu[10];
      Bool_t passConversion_electron_mumu[10];
      // eo id realted variables
      Bool_t isId_electron_mumu[10];
      Bool_t isIso_electron_mumu[10];
      Bool_t isEBEEGap_mumu[10];
      Double_t sf_electron_mumu[10];	

      // variables for muons
      Int_t nMuons_mumu;
      Double_t pt_muon_mumu[10];
      Double_t phi_muon_mumu[10];
      Double_t eta_muon_mumu[10];
      Double_t E_muon_mumu[10];
      Double_t vz_muon_mumu[10]; 
      Double_t v0_muon_mumu[10]; 
      Double_t d0_muon_mumu[10];
      Double_t d0BeamSpot_muon_mumu[10];
      Double_t chargedHadronIso_muon_mumu[10];
      Double_t neutralHadronIso_muon_mumu[10];
      Double_t photonIso_muon_mumu[10];
      Double_t pfIso_muon_mumu[10];
      Double_t sf_muon_mumu[10];
      Int_t charge_muon_mumu[10];
      Bool_t isId_muon_mumu[10];
      Bool_t isIso_muon_mumu[10];

      // variables for muonPairs 
      Int_t nMuonPairs_mumu; // if there is n muons there is (n*n - n)/2 distinct pairs
      Double_t deltaVz_mumu[10]; // max 5 muons -> max (25 -5)/2 = 10 muonPairs
      Double_t deltaV0_mumu[10];
      Double_t invMass_mumu[10];

      // event related variables
      Int_t run_num_mumu;
      Int_t event_num_mumu;
      Int_t lumi_num_mumu;
      Int_t nvtx_mumu;
      Int_t npu_mumu;
      Double_t evt_puSF_mumu;
      Double_t evt_met_mumu;

      // eo MyDoubleMuTree



      Int_t NEvent = datasets[d]->NofEvtsToRunOver();
      //	Double_t xs [1];
       



      // Define the bookkeeping tree
      /*
	TTree *bookkeeping = new TTree("startevents","startevents");
	//	bookkeeping->Branch("xs",xs,"xs/D");
	bookkeeping->Branch("NEvent",NEvent,"NEvent/I");
      */

	
      fout->cd();



      // bo of the tree with all the cuts applied for e-mu final state
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
      myTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
      myTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
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
      myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
      myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
      myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
      myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
      myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
      myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
      myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");


      // event related variables
      myTree->Branch("run_num",&run_num,"run_num/I");
      myTree->Branch("event_num",&event_num,"event_num/I");
      myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
      myTree->Branch("nvtx",&nvtx,"nvtx/I");
      myTree->Branch("npu",&npu,"npu/I");
      myTree->Branch("puSF",&puSF,"puSF/D");	
	
      // eo the tree with all the cuts applied for e-mu final state

	
      // bo a secondary tree that is filled before the whole list of cut
      TTree* myPreCutTree = new TTree("preCutTree","preCutTree");


      // electrons
      //	myPreCutTree->Branch("templatevar",templatevar,"templatevar[nElectrons_elel_pc]/D");
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
      myPreCutTree->Branch("isId_electron_pc",isId_electron_pc,"isId_electron_pc[nElectrons_pc]/O)");
      myPreCutTree->Branch("isIso_electron_pc",isIso_electron_pc,"isIso_electron_pc[nElectrons_pc]/O)");
      myPreCutTree->Branch("isEBEEGap_pc",isEBEEGap_pc,"isEBEEGap_pc[nElectrons_pc]/O)");
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
      myPreCutTree->Branch("isId_muon_pc",isId_muon_pc,"isId_muon_pc[nMuons_pc]/D");
      myPreCutTree->Branch("isIso_muon_pc",isIso_muon_pc,"isIso_muon_pc[nMuons_pc]/D");
      myPreCutTree->Branch("d0_muon_pc",d0_muon_pc,"d0_muon_pc[nMuons_pc]/D");
      myPreCutTree->Branch("d0BeamSpot_muon_pc",d0BeamSpot_muon_pc,"d0BeamSpot_muon_pc[nMuons_pc]/D");
      myPreCutTree->Branch("sf_muon_pc",sf_muon_pc,"sf_muon_pc[nMuons_pc]/D");


      // eo a secondary tree that is filled before the whole list of cut


      // bo a third tree that is filled if there is at least two electrons (el-el)

      TTree* myDoubleElTree = new TTree("doubleElTree","doubleElTree");


      // electrons
      //	myDoubleElTree->Branch("templatevar",templatevar,"templatevar[nElectrons_elel]/D");
      myDoubleElTree->Branch("nElectrons_elel",&nElectrons_elel, "nElectrons_elel/I");//                                                            
      myDoubleElTree->Branch("pt_electron_elel",pt_electron_elel,"pt_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("phi_electron_elel",phi_electron_elel,"phi_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("eta_electron_elel",eta_electron_elel,"eta_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("eta_superCluster_electron_elel",eta_superCluster_electron_elel,"eta_superCluster_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("E_electron_elel",E_electron_elel,"E_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("vz_electron_elel",vz_electron_elel,"vz_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("v0_electron_elel",v0_electron_elel,"v0_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("chargedHadronIso_electron_elel",chargedHadronIso_electron_elel,"chargedHadronIso_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("neutralHadronIso_electron_elel",neutralHadronIso_electron_elel,"neutralHadronIso_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("photonIso_electron_elel",photonIso_electron_elel,"photonIso_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("pfIso_electron_elel",pfIso_electron_elel,"pfIso_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("charge_electron_elel",charge_electron_elel,"charge_electron_elel[nElectrons_elel]/I");
      myDoubleElTree->Branch("d0_electron_elel",d0_electron_elel,"d0_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("d0BeamSpot_electron_elel",d0BeamSpot_electron_elel,"d0BeamSpot_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("sigmaIEtaIEta_electron_elel",sigmaIEtaIEta_electron_elel,"sigmaIEtaIEta_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("deltaEtaIn_electron_elel",deltaEtaIn_electron_elel,"deltaEtaIn_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("deltaPhiIn_electron_elel",deltaPhiIn_electron_elel,"deltaPhiIn_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("hadronicOverEm_electron_elel",hadronicOverEm_electron_elel,"hadronicOverEm_electron_elel[nElectrons_elel]/D");
      myDoubleElTree->Branch("missingHits_electron_elel",missingHits_electron_elel,"missingHits_electron_elel[nElectrons_elel]/I");
      myDoubleElTree->Branch("passConversion_electron_elel",passConversion_electron_elel,"passConversion_electron_elel[nElectrons_elel]/O)");
      myDoubleElTree->Branch("isId_electron_elel",isId_electron_elel,"isId_electron_elel[nElectrons_elel]/O)");
      myDoubleElTree->Branch("isIso_electron_elel",isIso_electron_elel,"isIso_electron_elel[nElectrons_elel]/O)");
      myDoubleElTree->Branch("isEBEEGap_elel",isEBEEGap_elel,"isEBEEGap_elel[nElectrons_elel]/O)");
      myDoubleElTree->Branch("sf_electron_elel",sf_electron_elel,"sf_electron_elel[nElectrons_elel]/D");
	
      // electronPairs
      //      myDoubleElTree->Branch("templatevar",templatevar,"templatevar[nElectronPairs_elel]/I"); 
      myDoubleElTree->Branch("nElectronPairs_elel",&nElectronPairs_elel, "nElectronPairs_elel/I");
      myDoubleElTree->Branch("deltaVz_elel",deltaVz_elel,"deltaVz_elel[nElectronPairs_elel]/D");
      myDoubleElTree->Branch("deltaV0_elel",deltaV0_elel,"deltaV0_elel[nElectronPairs_elel]/D");
      myDoubleElTree->Branch("invMass_elel",invMass_elel,"invMass_elel[nElectronPairs_elel]/D"); 



      // muons
      //        myDoubleElTree->Branch("templatevar",templatevar,"templatevar[nMuons_elel]/D");
      myDoubleElTree->Branch("nMuons_elel",&nMuons_elel, "nMuons_elel/I");
      myDoubleElTree->Branch("pt_muon_elel",pt_muon_elel,"pt_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("phi_muon_elel",phi_muon_elel,"phi_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("eta_muon_elel",eta_muon_elel,"eta_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("E_muon_elel",E_muon_elel,"E_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("chargedHadronIso_muon_elel",chargedHadronIso_muon_elel,"chargedHadronIso_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("neutralHadronIso_muon_elel",neutralHadronIso_muon_elel,"neutralHadronIso_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("photonIso_muon_elel",photonIso_muon_elel,"photonIso_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("pfIso_muon_elel",pfIso_muon_elel,"pfIso_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("charge_muon_elel",charge_muon_elel,"charge_muon_elel[nMuons_elel]/I");
      myDoubleElTree->Branch("isId_muon_elel",isId_muon_elel,"isId_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("isIso_muon_elel",isIso_muon_elel,"isIso_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("d0_muon_elel",d0_muon_elel,"d0_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("d0BeamSpot_muon_elel",d0BeamSpot_muon_elel,"d0BeamSpot_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("sf_muon_elel",sf_muon_elel,"sf_muon_elel[nMuons_elel]/D");
      myDoubleElTree->Branch("charge_muon_elel",charge_muon_elel,"charge_muon_elel[nMuons_elel]/D");

      //

      // event related variables
       myDoubleElTree->Branch("run_num_elel",&run_num,"run_num_elel/I");
       myDoubleElTree->Branch("event_num_elel",&event_num,"event_num_elel/I");
       myDoubleElTree->Branch("lumi_num_elel",&lumi_num_elel,"lumi_num_elel/I");
       myDoubleElTree->Branch("nvtx_elel",&nvtx_elel,"nvtx_elel/I");
       myDoubleElTree->Branch("npu_elel",&npu_elel,"npu_elel/I");
       myDoubleElTree->Branch("evt_puSF_elel",&evt_puSF_elel,"evt_puSF_elel/D");
       myDoubleElTree->Branch("evt_met_elel",&evt_met_elel,"evt_met_elel/D");

      

      // eo a third tree that is filled if there is at least two electrons (el-el)




      // bo a fourth tree that is filled if there is at least two muons (mu-mu)

      TTree* myDoubleMuTree = new TTree("doubleMuTree","doubleMuTree");


      // electrons
      //	myDoubleMuTree->Branch("templatevar",templatevar,"templatevar[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("nElectrons_mumu",&nElectrons_mumu, "nElectrons_mumu/I");//                                                            
      myDoubleMuTree->Branch("pt_electron_mumu",pt_electron_mumu,"pt_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("phi_electron_mumu",phi_electron_mumu,"phi_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("eta_electron_mumu",eta_electron_mumu,"eta_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("eta_superCluster_electron_mumu",eta_superCluster_electron_mumu,"eta_superCluster_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("E_electron_mumu",E_electron_mumu,"E_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("chargedHadronIso_electron_mumu",chargedHadronIso_electron_mumu,"chargedHadronIso_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("neutralHadronIso_electron_mumu",neutralHadronIso_electron_mumu,"neutralHadronIso_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("photonIso_electron_mumu",photonIso_electron_mumu,"photonIso_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("pfIso_electron_mumu",pfIso_electron_mumu,"pfIso_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("charge_electron_mumu",charge_electron_mumu,"charge_electron_mumu[nElectrons_mumu]/I");
      myDoubleMuTree->Branch("d0_electron_mumu",d0_electron_mumu,"d0_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("d0BeamSpot_electron_mumu",d0BeamSpot_electron_mumu,"d0BeamSpot_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("sigmaIEtaIEta_electron_mumu",sigmaIEtaIEta_electron_mumu,"sigmaIEtaIEta_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("deltaEtaIn_electron_mumu",deltaEtaIn_electron_mumu,"deltaEtaIn_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("deltaPhiIn_electron_mumu",deltaPhiIn_electron_mumu,"deltaPhiIn_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("hadronicOverEm_electron_mumu",hadronicOverEm_electron_mumu,"hadronicOverEm_electron_mumu[nElectrons_mumu]/D");
      myDoubleMuTree->Branch("missingHits_electron_mumu",missingHits_electron_mumu,"missingHits_electron_mumu[nElectrons_mumu]/I");
      myDoubleMuTree->Branch("passConversion_electron_mumu",passConversion_electron_mumu,"passConversion_electron_mumu[nElectrons_mumu]/O)");
      myDoubleMuTree->Branch("isId_electron_mumu",isId_electron_mumu,"isId_electron_mumu[nElectrons_mumu]/O)");
      myDoubleMuTree->Branch("isIso_electron_mumu",isIso_electron_mumu,"isIso_electron_mumu[nElectrons_mumu]/O)");
      myDoubleMuTree->Branch("isEBEEGap_mumu",isEBEEGap_mumu,"isEBEEGap_mumu[nElectrons_mumu]/O)");
      myDoubleMuTree->Branch("sf_electron_mumu",sf_electron_mumu,"sf_electron_mumu[nElectrons_mumu]/D");
	

      // muons
      //        myDoubleMuTree->Branch("templatevar",templatevar,"templatevar[nMuons_mumu]/D");
      myDoubleMuTree->Branch("nMuons_mumu",&nMuons_mumu, "nMuons_mumu/I");
      myDoubleMuTree->Branch("pt_muon_mumu",pt_muon_mumu,"pt_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("phi_muon_mumu",phi_muon_mumu,"phi_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("eta_muon_mumu",eta_muon_mumu,"eta_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("E_muon_mumu",E_muon_mumu,"E_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("vz_muon_mumu",vz_muon_mumu,"vz_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("v0_muon_mumu",v0_muon_mumu,"v0_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("chargedHadronIso_muon_mumu",chargedHadronIso_muon_mumu,"chargedHadronIso_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("neutralHadronIso_muon_mumu",neutralHadronIso_muon_mumu,"neutralHadronIso_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("photonIso_muon_mumu",photonIso_muon_mumu,"photonIso_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("pfIso_muon_mumu",pfIso_muon_mumu,"pfIso_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("charge_muon_mumu",charge_muon_mumu,"charge_muon_mumu[nMuons_mumu]/I");
      myDoubleMuTree->Branch("isId_muon_mumu",isId_muon_mumu,"isId_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("isIso_muon_mumu",isIso_muon_mumu,"isIso_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("d0_muon_mumu",d0_muon_mumu,"d0_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("d0BeamSpot_muon_mumu",d0BeamSpot_muon_mumu,"d0BeamSpot_muon_mumu[nMuons_mumu]/D");
      myDoubleMuTree->Branch("sf_muon_mumu",sf_muon_mumu,"sf_muon_mumu[nMuons_mumu]/D");


      // muonPairs
      //	myDoubleMuTree->Branch("templatevar",templatevar,"templatevar[nMuonPairs_mumu]/D"); 
      myDoubleMuTree->Branch("nMuonPairs_mumu",&nMuonPairs_mumu,"nMuonPairs_mumu/I"); 
      myDoubleMuTree->Branch("deltaVz_mumu",deltaVz_mumu,"deltaVz_mumu[nMuonPairs_mumu]/D"); 
      myDoubleMuTree->Branch("deltaV0_mumu",deltaV0_mumu,"deltaV0_mumu[nMuonPairs_mumu]/D"); 
      myDoubleMuTree->Branch("invMass_mumu",invMass_mumu,"invMass_mumu[nMuonPairs_mumu]/D"); 


      // event related variables 
      myDoubleMuTree->Branch("run_num_mumu",&run_num,"run_num_mumu/I");
      myDoubleMuTree->Branch("event_num_mumu",&event_num,"event_num_mumu/I");
      myDoubleMuTree->Branch("lumi_num_mumu",&lumi_num_mumu,"lumi_num_mumu/I");
      myDoubleMuTree->Branch("nvtx_mumu",&nvtx_mumu,"nvtx_mumu/I");
      myDoubleMuTree->Branch("npu_mumu",&npu_mumu,"npu_mumu/I");
      myDoubleMuTree->Branch("evt_puSF_mumu",&evt_puSF_mumu,"evt_puSF_mumu/D");
      myDoubleMuTree->Branch("evt_met_mumu",&evt_met_mumu,"evt_met_mumu/D");

	
      // eo a fourth tree that is filled if there is at least two muons (mu-mu)            



      //////////////////////////////////////////////////
      // Loop on events
      /////////////////////////////////////////////////

	
      int start = 0;
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
      vector<TRootElectron*> KynElectrons;
      vector<TRootElectron*> KynIdElectrons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootElectron*> selectedLooseElectrons;
      vector<TRootElectron*> selectedExtraElectrons;

      vector<TRootMuon*> KynMuons;
      vector<TRootMuon*> KynIdMuons;
      vector<TRootMuon*> selectedMuons;
      vector<TRootMuon*> selectedLooseMuons;
      vector<TRootMuon*> selectedExtraMuons;

      vector<TRootPFJet*>    selectedJets;
      vector<TRootSubstructureJet*>    selectedFatJets;


      selectedElectrons.reserve(10);
      selectedMuons.reserve(10);



      //////////////////////////////////////
      // Begin Event Loop
      //////////////////////////////////////
	
      for (unsigned int ievt = event_start; ievt < end_d; ievt++)
	{


	  if (debug) cout << "just entered the event loop!" << endl;

	  // Set default value for evertything that goes in the Tree
	  nMuons = nElectrons = 0;
	  nMuons_pc = nElectrons_pc = 0;
	  nMuons_elel = nElectrons_elel = 0;
	  nMuons_mumu = nElectrons_mumu = 0;
	    
	    
	  double ievt_d = ievt;
	  currentfrac = ievt_d/end_d;
	  if (debug)cout << endl << endl << "Starting a new event loop!"<<endl;
	    
	  if(ievt%1000 == 0)
	    {
	      std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
	    }
	    

	  // loading event
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, init_fatjets,  mets, debug); 
	  datasets[d]->eventTree()->LoadTree(ievt);
	  string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	  int currentRun = event->runId();

	    
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
	  //-----

	  //-----



	  // setting event related variables
	  evt_met_elel = evt_met_mumu = mets[0]->Et();
	  event_num = event_num_elel = event_num_mumu =event->eventId();
	  run_num = run_num_elel = run_num_mumu =event->runId();
	  lumi_num = lumi_num_elel =lumi_num_mumu = event->lumiBlockId();
	  nvtx = nvtx_elel = nvtx_mumu = vertex.size();
	  npu = npu_elel = npu_mumu = (int)event->nTruePU();

	  
	  

	  ///////////////////////////////////////////
	  //  Trigger
	  ///////////////////////////////////////////

	  bool trigged = false;
	  bool fileChanged = false;
	  bool runChanged = false;


	  
	  if ( ! applyTriggers && previousFilename != currentFilename )
	    {
	      fileChanged = true;
	      previousFilename = currentFilename;
	      iFile++;
	      cout << "File changed!!! => iFile = " << iFile << endl;
	    }
      
	  if (applyTriggers)
	    {
	      trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTriggers);
	      trigged = trigger->checkIfFired();
        
	      if (trigged) {
		if (debug)cout << "event " << ievt << " was Trigged!!" << endl;
	      }
	      else {
		if (debug) cout << "event " << ievt << " was not trigged. Skiping event.." << endl;
		eventCount++;
		skippedEvent++;
		continue;
		
	      }
	    }  




	  ///////////////////////////////////////////
	  //  Pile up Scale Factor
	  ///////////////////////////////////////////

	  double lumiWeight = LumiWeights.ITweight( npu ); // simplest reweighting, just use reconstructed number of PV.
	  puSF = evt_puSF_elel = evt_puSF_mumu = lumiWeight;
	  if (isData) puSF =1;


	  ///////////////////////
	  // JER smearing
	  //////////////////////

	  ///////////////////////////////////////////////////////////
	  // Event selection
	  ///////////////////////////////////////////////////////////


	  // Declare selection instance
	  Run2Selection selection(init_jets, init_fatjets, init_muons, init_electrons, mets);


	  // jet selections
	  if (debug)cout<<"Getting Jets"<<endl;
	  selectedJets = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()



	    
	  // muons selections
	  if (debug)cout<<"Getting Muons"<<endl;
	  // make three collection of muons for the synch exercise
	  KynMuons = selection.GetSelectedDisplacedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, false, false); // pt, eta,
	  KynIdMuons = selection.GetSelectedDisplacedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, false, true); // pt, eta, id
	  selectedMuons = selection.GetSelectedDisplacedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, true, true); // pt, eta, id, iso
	  selectedLooseMuons = selection.GetSelectedDisplacedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, true, true); // pt, eta, id, iso
	  //	    selectedMuons = init_muons;

	  //selectedMuons = selection.GetSelectedDisplacedMuons(30., 2., 0.15, true, true); // pt, eta, iso // run normally

	  //	    selectedMuons = selection.GetSelectedMuons();
	  //selectedMuons = init_muons;


	  // electrons selecttions
	  // make a new collections of electrons
	  if (debug)cout<<"Getting Electrons"<<endl;

	  // make three collection of muons for the synch exercise
	  KynElectrons = selection.GetSelectedDisplacedElectrons(el_pt_cut, el_eta_cut, false, false);// pt, eta
	  KynIdElectrons = selection.GetSelectedDisplacedElectrons(el_pt_cut, el_eta_cut, false, true);// pt, eta, id
	  selectedElectrons = selection.GetSelectedDisplacedElectrons(el_pt_cut, el_eta_cut, true, true);// pt, eta, id, iso
	  selectedLooseElectrons = selection.GetSelectedDisplacedElectrons(el_pt_cut, el_eta_cut, true, true);// pt, eta, id, iso
	  //selectedElectrons = init_electrons;
	  //selectedElectrons = selection.GetSelectedDisplacedElectrons();// pt, eta
	  //selectedElectrons = selection.GetSelectedElectrons();

	  //	  mets.TRootParticle();
	  

	  // fill TLorentz vector to allow easier calculation 	  

	  // electrons
	  vector<TLorentzVector> selectedElectronsTLV;
	  for(int iele=0; iele<selectedElectrons.size(); iele++)
	    {
	      selectedElectronsTLV.push_back(*selectedElectrons[iele]);
	    }
	  
	  // muons                                     
	  vector<TLorentzVector> selectedMuonsTLV;
	  for(int imuo=0; imuo<selectedMuons.size(); imuo++)
	    {
	      selectedMuonsTLV.push_back(*selectedMuons[imuo]);
	    }



	  int JetCut =0;
	  int nMu, nEl, nLooseIsoMu;


	  if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
	  // Apply primary vertex selection
	  bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);


	  eventCount++;

	  if (debug) cout <<"Initial number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  " << endl << init_muons.size() <<" "  <<init_electrons.size()<<" "<< selectedJets.size()   <<" "   <<" "<<JetCut  <<" "<<endl;


	  if (debug){
	    cout <<"Selected number of selected Electrons is " << selectedElectrons.size() << " and the number of selected muon is " << selectedMuons.size() << endl;
	    if (selectedElectrons.size()+ selectedMuons.size() > 0) cout << "At least one selected Lepton!!!!" << endl;
	  }


	  if (debug)	cout <<" applying baseline event selection..."<<endl;

	  ///////////////////////////////////
	  // Filling histograms / plotting //
	  ///////////////////////////////////

	  //            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*globalScaleFactor);



	  // bo assigning values to the e-mu tree 

	  //////////////////////////                                                                                                                                            
	  // jet Based Plots //
	  //////////////////////////
	  /*
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
              
            }
	  */
	    


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
	      // id
	      isId_electron[nElectrons]=false;
	      if( fabs(selectedElectrons[selel]->superClusterEta()) <= 1.479
		  && selectedElectrons[selel]->sigmaIEtaIEta() < 0.0101
		  && fabs(selectedElectrons[selel]->deltaEtaIn()) < 0.00926
		  && fabs(selectedElectrons[selel]->deltaPhiIn()) < 0.0336
		  && selectedElectrons[selel]->hadronicOverEm() < 0.0597
		  && fabs(1/selectedElectrons[selel]->E() - 1/selectedElectrons[selel]->P()) < 0.012 // to be fixed faco
		  && selectedElectrons[selel]->missingHits() <= 2 // check wrt to expectedMissingInnerHits        
		  && selectedElectrons[selel]->passConversion()){
		isId_electron[nElectrons]=true;
	      }
	      else if (fabs(selectedElectrons[selel]->superClusterEta()) < 2.5
		       && selectedElectrons[selel]->sigmaIEtaIEta() < 0.0279
		       && fabs(selectedElectrons[selel]->deltaEtaIn()) < 0.00724
		       && fabs(selectedElectrons[selel]->deltaPhiIn()) < 0.0918
		       && selectedElectrons[selel]->hadronicOverEm() < 0.0615
		       && fabs(1/selectedElectrons[selel]->E() - 1/selectedElectrons[selel]->P()) < 0.00999 // to be fixed faco
		       && selectedElectrons[selel]->missingHits() <= 1 // check wrt to expectedMissingInnerHits  
		       && selectedElectrons[selel]->passConversion()){
		isId_electron[nElectrons]=true;
	      }
	      // iso 
	      isIso_electron[nElectrons]=false;
	      if( fabs(selectedElectrons[selel]->superClusterEta()) <= 1.479){
		if(selectedElectrons[selel]->relPfIso(3,0) < 0.0354) isIso_electron[nElectrons]=true;
	      }
	      else if (fabs(selectedElectrons[selel]->superClusterEta()) < 2.5){
		if(selectedElectrons[selel]->relPfIso(3,0) < 0.0646) isIso_electron[nElectrons]=true;
	      }
		
	      isEBEEGap[nElectrons]=selectedElectrons[selel]->isEBEEGap();
	      // following code found in http://cmslxr.fnal.gov/source/RecoEgamma/PhotonIdentification/src/PhotonIsolationCalculator.cc#0520
	      /*
		isEBEEGap[nElectrons] = false;
		Double_t eta =  eta_superCluster_electron[nElectrons];
		Double_t feta = fabs(eta);
		if (fabs(feta-1.479)<0.1) isEBEEGap[nElectrons] = true ;
	      */
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
	      // id
	      isId_muon[nMuons]=false;
	      if( selectedMuons[selmu]->isGlobalMuon() && selectedMuons[selmu]->isPFMuon()
		  && selectedMuons[selmu]->chi2() < 10
		  && selectedMuons[selmu]->nofTrackerLayersWithMeasurement() > 5
		  &&  selectedMuons[selmu]->nofValidMuHits() > 0
		  && selectedMuons[selmu]->nofValidPixelHits() > 0 
		  && selectedMuons[selmu]->nofMatchedStations()> 1){
		isId_muon[nMuons]=true;
	      }
	      //iso
	      isIso_muon[nMuons]=false;
	      if ( (selectedMuons[selmu]->chargedHadronIso(4) + max( 0.0, selectedMuons[selmu]->neutralHadronIso(4) + selectedMuons[selmu]->photonIso(4) - 0.5*selectedMuons[selmu]->puChargedHadronIso(4) ) ) / selectedMuons[selmu]->Pt() < mu_iso_cut ) isIso_muon[nMuons]=true;
	      pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
	      charge_muon[nMuons]=selectedMuons[selmu]->charge();
	      sf_muon[nMuons]=muonSFWeightIso_TT->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID_T->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons << " and pt equals to " << pt_muon[nMuons] << endl;
	      nMuons++;

	    }

	  // eo assigning values to the e-mu tree 



	  // bo assigning values to the precut tree

	  // Precut Tree
	    
	  //////////////////////
	  // Jets Based Plots //
	  //////////////////////
	  if (debug) cout << "before jets loop" << endl;
	    
	  for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
	    {
	    }


	  //////////////////////////
	  // Electron Based Plots //
	  //////////////////////////
	  if (debug) cout << "before electrons pc loop" << endl;

	    
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
	      isEBEEGap_pc[nElectrons_pc]=init_electrons[initel]->isEBEEGap();

	      // id
	      isId_electron_pc[nElectrons_pc]=false;
	      if( fabs(init_electrons[initel]->superClusterEta()) <= 1.479){
		if (init_electrons[initel]->sigmaIEtaIEta() < 0.0101
		    && fabs(init_electrons[initel]->deltaEtaIn()) < 0.00926
		    && fabs(init_electrons[initel]->deltaPhiIn()) < 0.0336
		    && init_electrons[initel]->hadronicOverEm() < 0.0597
		    && fabs(1/init_electrons[initel]->E() - 1/init_electrons[initel]->P()) < 0.012
		    && init_electrons[initel]->missingHits() <= 2 // check wrt to expectedMissingInnerHits        
		    && init_electrons[initel]->passConversion()){
		  isId_electron_pc[nElectrons_pc]=true;
		}
	      }
	      else if (fabs(init_electrons[initel]->superClusterEta()) < 2.5){
		if ( init_electrons[initel]->sigmaIEtaIEta() < 0.0279
		     && fabs(init_electrons[initel]->deltaEtaIn()) < 0.00724
		     && fabs(init_electrons[initel]->deltaPhiIn()) < 0.0918
		     && init_electrons[initel]->hadronicOverEm() < 0.0615
		     && fabs(1/init_electrons[initel]->E() - 1/init_electrons[initel]->P()) < 0.00999
		     && init_electrons[initel]->missingHits() <= 1 // check wrt to expectedMissingInnerHits  
		     && init_electrons[initel]->passConversion()){
		  isId_electron_pc[nElectrons_pc]=true;
		}
	      }
	      // iso to be checked!!! Make sure what is this function getting! faco
	      isIso_electron_pc[nElectrons_pc]=false;
	      if( fabs(init_electrons[initel]->superClusterEta()) <= 1.479){
		if(init_electrons[initel]->relPfIso(3,0) < 0.0354) isIso_electron_pc[nElectrons_pc]=true;
	      }
	      else if (fabs(init_electrons[initel]->superClusterEta()) < 2.5){
		if(init_electrons[initel]->relPfIso(3,0) < 0.0646) isIso_electron_pc[nElectrons_pc]=true;
	      }
		
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
	      // id
	      isId_muon_pc[nMuons_pc]=false;
	      if( init_muons[initmu]->isGlobalMuon() && init_muons[initmu]->isPFMuon()
		  && init_muons[initmu]->chi2() < 10
		  && init_muons[initmu]->nofTrackerLayersWithMeasurement() > 5
		  &&  init_muons[initmu]->nofValidMuHits() > 0
		  && init_muons[initmu]->nofValidPixelHits() > 0 
		  && init_muons[initmu]->nofMatchedStations()> 1){
		isId_muon_pc[nMuons_pc]=true;
	      }
	      //iso
	      isIso_muon_pc[nMuons_pc]=false;
	      if ( (init_muons[initmu]->chargedHadronIso(4) + max( 0.0, init_muons[initmu]->neutralHadronIso(4) + init_muons[initmu]->photonIso(4) - 0.5*init_muons[initmu]->puChargedHadronIso(4) ) ) / init_muons[initmu]->Pt() < mu_iso_cut ) isIso_muon_pc[nMuons_pc]=true;

	      //	      sf_muon_pc[nMuons_pc]=muonSFWeightID_T->at(init_muons[initmu]->Eta(), init_muons[initmu]->Pt(), 0)* muonSFWeightIso_TT->at(init_muons[initmu]->Eta(), init_muons[initmu]->Pt(), 0);
	      //	      sf_muon_pc[nMuons_pc]=muonSFWeightID_T->at(init_muons[initmu]->Eta(), init_muons[initmu]->Pt(), 0);//* muonSFWeightIso_TT->at(init_muons[initmu]->Eta(), init_muons[initmu]->Pt(), 0);
	      

	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons_pc << " and pt equals to " << pt_muon_pc[nMuons_pc] << endl;
	      nMuons_pc++;
	      
	      
	    }
	    
	  // eo  assigning values to the precut tree
	    

	    
	  // bo assigning values to the elel tree

	  //////////////////////////
	  ///// jet Based Plots ////
	  //////////////////////////
	    
	  /*
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
              
            }
	  */
	    

	  //////////////////////////
	  // Electron Based Plots //
	  //////////////////////////
	  if (debug) cout << "before electrons loop" << endl;

	    
	  nElectrons_elel=0;
	  for (Int_t selel =0; selel < selectedLooseElectrons.size() && selel < 10; selel++ )
	    {
	      pt_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->Pt();
	      phi_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->Phi();
	      eta_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->Eta();
	      eta_superCluster_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->superClusterEta();
	      E_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->E();
	      vz_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->vz();
	      v0_electron_elel[nElectrons_elel]=sqrt(pow (selectedLooseElectrons[selel]->vx(),2 )+pow(selectedLooseElectrons[selel]->vy(),2) );
	      d0_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->d0();
	      d0BeamSpot_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->d0BeamSpot();
	      chargedHadronIso_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->chargedHadronIso(3);
	      neutralHadronIso_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->neutralHadronIso(3);
	      photonIso_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->photonIso(3);
	      pfIso_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->relPfIso(3,0);
	      charge_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->charge();
	      sigmaIEtaIEta_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->sigmaIEtaIEta();
	      deltaEtaIn_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->deltaEtaIn();
	      deltaPhiIn_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->deltaPhiIn();
	      hadronicOverEm_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->hadronicOverEm();
	      missingHits_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->missingHits();
	      passConversion_electron_elel[nElectrons_elel]=selectedLooseElectrons[selel]->passConversion();
	      // id
	      isId_electron_elel[nElectrons_elel]=false;
	      if( fabs(selectedLooseElectrons[selel]->superClusterEta()) <= 1.479
		  && selectedLooseElectrons[selel]->sigmaIEtaIEta() < 0.0101
		  && fabs(selectedLooseElectrons[selel]->deltaEtaIn()) < 0.00926
		  && fabs(selectedLooseElectrons[selel]->deltaPhiIn()) < 0.0336
		  && selectedLooseElectrons[selel]->hadronicOverEm() < 0.0597
		  && fabs(1/selectedLooseElectrons[selel]->E() - 1/selectedLooseElectrons[selel]->P()) < 0.012
		  && selectedLooseElectrons[selel]->missingHits() <= 2 // check wrt to expectedMissingInnerHits        
		  && selectedLooseElectrons[selel]->passConversion()){
		isId_electron_elel[nElectrons_elel]=true;
	      }
	      else if (fabs(selectedLooseElectrons[selel]->superClusterEta()) < 2.5
		       && selectedLooseElectrons[selel]->sigmaIEtaIEta() < 0.0279
		       && fabs(selectedLooseElectrons[selel]->deltaEtaIn()) < 0.00724
		       && fabs(selectedLooseElectrons[selel]->deltaPhiIn()) < 0.0918
		       && selectedLooseElectrons[selel]->hadronicOverEm() < 0.0615
		       && fabs(1/selectedLooseElectrons[selel]->E() - 1/selectedLooseElectrons[selel]->P()) < 0.00999
		       && selectedLooseElectrons[selel]->missingHits() <= 1 // check wrt to expectedMissingInnerHits  
		       && selectedLooseElectrons[selel]->passConversion()){
		isId_electron_elel[nElectrons_elel]=true;
	      }
	      // iso to be checked!!! Make sure what is this function getting! faco
	      isIso_electron_elel[nElectrons_elel]=false;
	      if( fabs(selectedLooseElectrons[selel]->superClusterEta()) <= 1.479){
		if(selectedLooseElectrons[selel]->relPfIso(3,0) < 0.0354) isIso_electron_elel[nElectrons_elel]=true;
	      }
	      else if (fabs(selectedLooseElectrons[selel]->superClusterEta()) < 2.5){
		if(selectedLooseElectrons[selel]->relPfIso(3,0) < 0.0646) isIso_electron_elel[nElectrons_elel]=true;
	      }
		
	      isEBEEGap_elel[nElectrons_elel]=selectedLooseElectrons[selel]->isEBEEGap();
	      sf_electron_elel[nElectrons_elel]=electronSFWeight_->at(selectedLooseElectrons[selel]->Eta(),selectedLooseElectrons[selel]->Pt(),0);
	      if (debug) cout << "in electrons loops, nelectrons equals to " << nElectrons << " and pt equals to " << pt_electron_elel[nElectrons_elel] << endl;
	      nElectrons_elel++;
            }

	  /////////////////////////
	  // ElectronPairs Plots //
	  /////////////////////////
	  // to count distinguishable pairs of the same object we need the lowest pt to go from the last object to the second
	  // then we compare the hight pt object from the first to the second-1
	  // example with 5 object "e"
	  //  (e0;e4) (e1;e4) (e2;e4) (e3;e4)
	  //  (e0;e3) (e1;e3) (e2;e3)
	  //  (e0;e2) (e1;e2) 
	  //  (e0;e1)
	    
	    
	  nElectronPairs_elel=0;
	  for (Int_t secondEl = nElectrons_elel-1; secondEl > 0; secondEl-- )
	    {
	      if (debug) cout << "secondEl is " << secondEl << endl;
	      for (Int_t firstEl = 0; firstEl < secondEl ; firstEl++ )
		{
		  deltaVz_elel[nElectronPairs_elel]=abs(vz_electron_elel[firstEl]-vz_electron_elel[secondEl]);
		  deltaV0_elel[nElectronPairs_elel]=abs(v0_electron_elel[firstEl]-v0_electron_elel[secondEl]);
		  invMass_elel[nElectronPairs_elel]=(selectedElectronsTLV[firstEl] + selectedElectronsTLV[secondEl]).M();
		  nElectronPairs_elel++;
		    
		  // debug statement
		  if (debug){
		    cout << "firstEl is " << firstEl << endl;
		    cout << "vz_electron_elel[firstEl] is" << vz_electron_elel[firstEl] << endl;
		    cout << "vz_electron_elel[secondEl] is" << vz_electron_elel[secondEl] << endl;
		    cout << "abs(vz_electron_elel[firstEl]-vz_electron_elel[secondEl] is" << abs(vz_electron_elel[firstEl]-vz_electron_elel[secondEl]) << endl;
		  }
		
		}
	    }



	  //////////////////////
	  // Muon Based Plots //
	  //////////////////////
	    
	  nMuons_elel=0;
	  for (Int_t selmu =0; selmu < selectedLooseMuons.size() && selmu < 10; selmu++ )
            {
	      pt_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->Pt();
	      phi_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->Phi();
	      eta_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->Eta();
	      E_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->E();
	      d0_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->d0();
	      d0BeamSpot_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->d0BeamSpot();
	      chargedHadronIso_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->chargedHadronIso(4);
	      neutralHadronIso_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->neutralHadronIso(4);
	      photonIso_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->photonIso(4);
	      // id
	      isId_muon_elel[nMuons_elel]=false;
	      if( selectedLooseMuons[selmu]->isGlobalMuon() && selectedLooseMuons[selmu]->isPFMuon()
		  && selectedLooseMuons[selmu]->chi2() < 10
		  && selectedLooseMuons[selmu]->nofTrackerLayersWithMeasurement() > 5
		  &&  selectedLooseMuons[selmu]->nofValidMuHits() > 0
		  && selectedLooseMuons[selmu]->nofValidPixelHits() > 0 //no more in the twiki (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon)
		  && selectedLooseMuons[selmu]->nofMatchedStations()> 1){
		isId_muon_elel[nMuons_elel]=true;
	      }
	      //iso
	      isIso_muon_elel[nMuons_elel]=false;
	      if ( (selectedLooseMuons[selmu]->chargedHadronIso(4) + max( 0.0, selectedLooseMuons[selmu]->neutralHadronIso(4) + selectedLooseMuons[selmu]->photonIso(4) - 0.5*selectedLooseMuons[selmu]->puChargedHadronIso(4) ) ) / selectedLooseMuons[selmu]->Pt() < mu_iso_cut ) isIso_muon_elel[nMuons_elel]=true;
	      pfIso_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->relPfIso(4,0);
	      charge_muon_elel[nMuons_elel]=selectedLooseMuons[selmu]->charge();
	      sf_muon_elel[nMuons_elel]=muonSFWeightIso_TT->at(selectedLooseMuons[selmu]->Eta(), selectedLooseMuons[selmu]->Pt(), 0)* muonSFWeightID_T->at(selectedLooseMuons[selmu]->Eta(), selectedLooseMuons[selmu]->Pt(), 0);
	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons << " and pt equals to " << pt_muon_elel[nMuons_elel] << endl;
	      nMuons_elel++;

	    }

	  // eo assigning values to the elel tree



	  // bo assigning values to the mumu tree

	  //////////////////////////
	  ///// jet Based Plots ////
	  //////////////////////////
	    
	  /*
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
              
            }
	  */
	    

	  //////////////////////////
	  // Electron Based Plots //
	  //////////////////////////
	  if (debug) cout << "before electrons loop" << endl;

	    
	  nElectrons_mumu=0;
	  for (Int_t selel =0; selel < selectedLooseElectrons.size() && selel < 10; selel++ )
	    {
	      pt_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->Pt();
	      phi_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->Phi();
	      eta_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->Eta();
	      eta_superCluster_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->superClusterEta();
	      E_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->E();
	      d0_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->d0();
	      d0BeamSpot_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->d0BeamSpot();
	      chargedHadronIso_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->chargedHadronIso(3);
	      neutralHadronIso_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->neutralHadronIso(3);
	      photonIso_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->photonIso(3);
	      pfIso_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->relPfIso(3,0);
	      charge_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->charge();
	      sigmaIEtaIEta_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->sigmaIEtaIEta();
	      deltaEtaIn_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->deltaEtaIn();
	      deltaPhiIn_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->deltaPhiIn();
	      hadronicOverEm_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->hadronicOverEm();
	      missingHits_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->missingHits();
	      passConversion_electron_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->passConversion();
	      // id
	      isId_electron_mumu[nElectrons_mumu]=false;
	      if( fabs(selectedLooseElectrons[selel]->superClusterEta()) <= 1.479
		  && selectedLooseElectrons[selel]->sigmaIEtaIEta() < 0.0101
		  && fabs(selectedLooseElectrons[selel]->deltaEtaIn()) < 0.00926
		  && fabs(selectedLooseElectrons[selel]->deltaPhiIn()) < 0.0336
		  && selectedLooseElectrons[selel]->hadronicOverEm() < 0.0597
		  && fabs(1/selectedLooseElectrons[selel]->E() - 1/selectedLooseElectrons[selel]->P()) < 0.012
		  && selectedLooseElectrons[selel]->missingHits() <= 2 // check wrt to expectedMissingInnerHits        
		  && selectedLooseElectrons[selel]->passConversion()){
		isId_electron_mumu[nElectrons_mumu]=true;
	      }
	      else if (fabs(selectedLooseElectrons[selel]->superClusterEta()) < 2.5
		       && selectedLooseElectrons[selel]->sigmaIEtaIEta() < 0.0279
		       && fabs(selectedLooseElectrons[selel]->deltaEtaIn()) < 0.00724
		       && fabs(selectedLooseElectrons[selel]->deltaPhiIn()) < 0.0918
		       && selectedLooseElectrons[selel]->hadronicOverEm() < 0.0615
		       && fabs(1/selectedLooseElectrons[selel]->E() - 1/selectedLooseElectrons[selel]->P()) < 0.00999
		       && selectedLooseElectrons[selel]->missingHits() <= 1 // check wrt to expectedMissingInnerHits  
		       && selectedLooseElectrons[selel]->passConversion()){
		isId_electron_mumu[nElectrons_mumu]=true;
	      }
	      // iso to be checked!!! Make sure what is this function getting! faco
	      isIso_electron_mumu[nElectrons_mumu]=false;
	      if( fabs(selectedLooseElectrons[selel]->superClusterEta()) <= 1.479){
		if(selectedLooseElectrons[selel]->relPfIso(3,0) < 0.0354) isIso_electron_mumu[nElectrons_mumu]=true;
	      }
	      else if (fabs(selectedLooseElectrons[selel]->superClusterEta()) < 2.5){
		if(selectedLooseElectrons[selel]->relPfIso(3,0) < 0.0646) isIso_electron_mumu[nElectrons_mumu]=true;
	      }
		
	      isEBEEGap_mumu[nElectrons_mumu]=selectedLooseElectrons[selel]->isEBEEGap();
	      // following code found in http://cmslxr.fnal.gov/source/RecoEgamma/PhotonIdentification/src/PhotonIsolationCalculator.cc#0520
	      /*
		isEBEEGap_mumu[nElectrons_mumu] = false;
		Double_t eta =  eta_superCluster_electron_mumu[nElectrons_mumu];
		Double_t feta = fabs(eta);
		if (fabs(feta-1.479)<0.1) isEBEEGap_mumu[nElectrons_mumu] = true ;
	      */
	      sf_electron_mumu[nElectrons_mumu]=electronSFWeight_->at(selectedLooseElectrons[selel]->Eta(),selectedLooseElectrons[selel]->Pt(),0);
	      if (debug) cout << "in electrons loops, nelectrons equals to " << nElectrons << " and pt equals to " << pt_electron_mumu[nElectrons_mumu] << endl;
	      nElectrons_mumu++;
            }


	  //////////////////////
	  // Muon Based Plots //
	  //////////////////////
	    
	  nMuons_mumu=0;
	  for (Int_t selmu =0; selmu < selectedLooseMuons.size() && selmu < 10; selmu++ )
            {
	      pt_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->Pt();
	      phi_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->Phi();
	      eta_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->Eta();
	      E_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->E();
	      vz_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->vz();
	      v0_muon_mumu[nMuons_mumu]=sqrt( pow(selectedLooseMuons[selmu]->vx(), 2) + pow(selectedLooseMuons[selmu]->vy(), 2) );
	      d0_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->d0();
	      d0BeamSpot_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->d0BeamSpot();
	      chargedHadronIso_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->chargedHadronIso(4);
	      neutralHadronIso_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->neutralHadronIso(4);
	      photonIso_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->photonIso(4);
	      // id
	      isId_muon_mumu[nMuons_mumu]=false;
	      if( selectedLooseMuons[selmu]->isGlobalMuon() && selectedLooseMuons[selmu]->isPFMuon()
		  && selectedLooseMuons[selmu]->chi2() < 10
		  && selectedLooseMuons[selmu]->nofTrackerLayersWithMeasurement() > 5
		  &&  selectedLooseMuons[selmu]->nofValidMuHits() > 0
		  && selectedLooseMuons[selmu]->nofValidPixelHits() > 0 //no more in the twiki (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon)
		  && selectedLooseMuons[selmu]->nofMatchedStations()> 1){
		isId_muon_mumu[nMuons_mumu]=true;
	      }
	      //iso
	      isIso_muon_mumu[nMuons_mumu]=false;
	      if ( (selectedLooseMuons[selmu]->chargedHadronIso(4) + max( 0.0, selectedLooseMuons[selmu]->neutralHadronIso(4) + selectedLooseMuons[selmu]->photonIso(4) - 0.5*selectedLooseMuons[selmu]->puChargedHadronIso(4) ) ) / selectedLooseMuons[selmu]->Pt() < mu_iso_cut ) isIso_muon_mumu[nMuons_mumu]=true;
	      pfIso_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->relPfIso(4,0); // check wrt formula!!! faco
	      charge_muon_mumu[nMuons_mumu]=selectedLooseMuons[selmu]->charge();
	      sf_muon_mumu[nMuons_mumu]=muonSFWeightIso_TT->at(selectedLooseMuons[selmu]->Eta(), selectedLooseMuons[selmu]->Pt(), 0)* muonSFWeightID_T->at(selectedLooseMuons[selmu]->Eta(), selectedLooseMuons[selmu]->Pt(), 0);
	      if (debug) cout << "in muons loops, nmuons equals to " << nMuons << " and pt equals to " << pt_muon_mumu[nMuons_mumu] << endl;
	      nMuons_mumu++;

	    }


	  ///////////////////////////
	  // MuonPairs Based Plots //
	  ///////////////////////////

	  nMuonPairs_mumu=0;
	  for (Int_t secondMu = nMuons_mumu-1; secondMu > 0; secondMu-- )
	    {
	      if (debug) cout << "secondMu is " << secondMu << endl;
	      for (Int_t firstMu = 0; firstMu < secondMu ; firstMu++ )
		{
		  deltaVz_mumu[nMuonPairs_mumu]=abs(vz_muon_mumu[firstMu]-vz_muon_mumu[secondMu]);
		  deltaV0_mumu[nMuonPairs_mumu]=abs(v0_muon_mumu[firstMu]-v0_muon_mumu[secondMu]);
		  invMass_mumu[nMuonPairs_mumu]=(selectedMuonsTLV[firstMu] + selectedMuonsTLV[secondMu]).M();
		  nMuonPairs_mumu++;
		    
		  // debug statement
		  if (debug){
		    cout << "firstMu is " << firstMu << endl;
		    cout << "vz_muon_mumu[firstMu] is" << vz_muon_mumu[firstMu] << endl;
		    cout << "vz_muon_mumu[secondMu] is" << vz_muon_mumu[secondMu] << endl;
		    cout << "abs(vz_muon_mumu[firstMu]-vz_muon_mumu[secondMu] is" << abs(vz_muon_mumu[firstMu]-vz_muon_mumu[secondMu]) << endl;
		  }
		
		}
	    }

	    


	  // eo assigning values to the elel tree

	    
	    
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






	  // -------------------------------
	  // bo filling the cutflow table --
	  // -------------------------------
 


	  // Synch cut 1: pt, eta, veto, DR                                                                                                                           
	  // Declare one bool per cut
	  /// ---
	  
	  // cut on electrons
	  Bool_t passedEcalCrackVeto = false;
	  Bool_t passedGoodEl =false;
	  Bool_t passedBlindingEl =false;

	  Bool_t passedPtEl = false;
	  Bool_t passedEtaEl = false;
	  Bool_t passedIdEl = false;
	  Bool_t passedIsoEl = false;

	  
	  //cut on muons
	  Bool_t passedGoodMu = false;
	  Bool_t passedBlindingMu =false;

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
	  
	  /// ---

	  // -------------------------------------------------------------------
	  // bo logic to set boolean that will be used to fill cutflow tables --
	  // -------------------------------------------------------------------


	  
	  
	  /*

	  // vector of bool
	  vector <Bool_t> CutFlow_oneEl_ind;
	  for (int i_tab = 1; i_tab < CutFlow_oneEl.size()-1; i_tab++){
	  CutFlow_oneEl_ind[i_tab]=false;
	  cout << CutFlow_oneEl_ind[i_tab] << endl;
	  }
	  */
	  
	  // logic for single efficiencies - electrons

	  //	  /*
	  for(int iele=0; iele<init_electrons.size(); iele++){
	    if (init_electrons[iele]->Pt() > el_pt_cut) passedPtEl=true;
	    if (!isEBEEGap_pc[iele]) passedEcalCrackVeto=true;
	    if (abs(init_electrons[iele]->Eta()) < el_eta_cut) passedEtaEl=true;
	    if (isId_electron_pc[iele]) passedIdEl=true;
	    if (isIso_electron_pc[iele]) passedIsoEl=true;
	    if (init_electrons[iele]->d0BeamSpot() < el_d0_cut) passedBlindingEl=true;
	  }
	  //	  */

	  /*
	  // logic for cumulative efficiencies - electrons
	  for(int iele=0; iele<init_electrons.size(); iele++){
	  if (init_electrons[iele]->Pt() > el_pt_cut) {
	  passedPtEl=true;
	  if (init_electrons[iele]->Pt() > el_pt_cut && !isEBEEGap_pc[iele]){
	  passedEcalCrackVeto=true;
	  if (init_electrons[iele]->Pt() > el_pt_cut && !isEBEEGap_pc[iele] && abs(init_electrons[iele]->Eta()) < el_eta_cut){
	  passedEtaEl=true;		  
	  if (init_electrons[iele]->Pt() > el_pt_cut && !isEBEEGap_pc[iele] && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele]) {
	  passedIdEl=true;
	  if (init_electrons[iele]->Pt() > el_pt_cut && !isEBEEGap_pc[iele] && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele] && isIso_electron_pc[iele]) {
	  passedIsoEl=true;
	  if (init_electrons[iele]->Pt() > el_pt_cut && !isEBEEGap_pc[iele] && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele] && isIso_electron_pc[iele] && init_electrons[iele]->d0BeamSpot() < el_d0_cut) {
	  passedBlindingEl=true;
	  }
	  }
	  }
	  }
	  }
	  }
	  }
	  */
	    
	  // filling single efficiencies - electrons

	  CutFlow_oneElTable.Fill(d,0,1);
	  if (passedPtEl)
	    CutFlow_oneElTable.Fill(d,1,1);
	  if (passedEcalCrackVeto)
	    CutFlow_oneElTable.Fill(d,2,1);
	  if (passedEtaEl)
	    CutFlow_oneElTable.Fill(d,3,1);
	  if (passedIdEl)
	    CutFlow_oneElTable.Fill(d,4,1);
	  if (passedIsoEl)
	    CutFlow_oneElTable.Fill(d,5,1);
	  if (passedBlindingEl)
	    CutFlow_oneElTable.Fill(d,6,1);

	  // reset all bool to false 
	  passedPtEl=false;
	  passedEcalCrackVeto=false;
	  passedEtaEl=false;
	  passedIdEl=false;
	  passedIsoEl=false;
	  passedBlindingEl=false;
	    

	  // logic for single efficiencies - muons
	    
	  //	    /*
	  for(int imuo=0; imuo<init_muons.size(); imuo++){
	    if (init_muons[imuo]->Pt() > mu_pt_cut) passedPtMu=true;
	    if (abs(init_muons[imuo]->Eta()) < mu_eta_cut) passedEtaMu=true;
	    if (isId_muon_pc[imuo]) passedIdMu=true;
	    if (isIso_muon_pc[imuo]) passedIsoMu=true;
	    if (init_muons[imuo]->d0BeamSpot() < mu_d0_cut) passedBlindingMu=true;
	  }
	  //	    */


	  /*
	  // logic for cumulative efficiencies - muons
	    
	  for(int imuo=0; imuo<init_muons.size(); imuo++){
	  if (init_muons[imuo]->Pt() > mu_pt_cut) {
	  passedPtMu=true;
	  if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut){
	  passedEtaMu=true;
	  if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut && isId_muon_pc[imuo]){
	  passedIdMu=true;
	  if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut && isId_muon_pc[imuo] && isIso_muon_pc[imuo]){
	  passedIsoMu=true;
	  if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut && isId_muon_pc[imuo] && isIso_muon_pc[imuo] && init_muons[imuo]->d0BeamSpot() < mu_d0_cut){
	  passedBlindingMu=true;
	  }
	  }
	  }
	  }
	  }
	  }
	  */

	  // filling single efficiencies - muons

	  CutFlow_oneMuTable.Fill(d,0,1);
	  if (passedPtMu)
	    CutFlow_oneMuTable.Fill(d,1,1);
	  if (passedEtaMu)
	    CutFlow_oneMuTable.Fill(d,2,1);
	  if (passedIdMu)
	    CutFlow_oneMuTable.Fill(d,3,1);
	  if (passedIsoMu)
	    CutFlow_oneMuTable.Fill(d,4,1);
	  if (passedBlindingMu)
	    CutFlow_oneMuTable.Fill(d,5,1);

	  // reset all bool to false
	  passedPtMu=false;
	  passedEtaMu=false;
	  passedIdMu=false;
	  passedIsoMu=false;
	  passedBlindingMu=false;


	  /*		      

	  vector<TLorentzVector> electronsTLV;
	  for(int iele=0; iele<init_electrons.size() && nElectrons<10; iele++)
	    {
	      init_electronsTLV.push_back(*init_electrons[iele]);
	    }
	  
	  vector<TLorentzVector> init_muonsTLV;
	  for(int imuo=0; imuo<init_muons.size() && nMuons<10; imuo++)
	    {
	      init_muonsTLV.push_back(*init_muons[imuo]);
	    }
	  
	  
	  vector<TLorentzVector> postCut_electronsTLV;
	  for(int iele=0; iele<selectedElectrons.size(); iele++)
	    {
	      postCut_electronsTLV.push_back(*selectedElectrons[iele]);
	    }
	  
	  // muons                                     
	  vector<TLorentzVector> postCut_muonsTLV;
	  for(int imuo=0; imuo<selectedMuons.size(); imuo++)
	    {
	      postCut_muonsTLV.push_back(*selectedMuons[imuo]);//fill a new vector with the muons that passed the cuts
	    }
	  */



	  // logic for full Synch cutflow


	  for(int iele=0; iele<init_electrons.size(); iele++){
	    if (init_electrons[iele]->Pt() > el_pt_cut){ // pt
	      passedPtEl=true;
	      //if (!isEBEEGap[iele]){//ecal crack veto
	      if (1){
		passedEcalCrackVeto=true;
		if (init_electrons[iele]->Pt() > el_pt_cut && abs(init_electrons[iele]->Eta()) < el_eta_cut){ // eta
		  passedEtaEl=true;
		  // try different why to, supposedly, exactly the same either by getting the value of the tree or by using the function of the collection
		  if (init_electrons[iele]->Pt() > el_pt_cut && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele]){ // id
		    //		    if (pt_electron_pc[iele] > el_pt_cut && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele]){ // id
		    //		    if (init_electrons[iele]->Pt() > el_pt_cut && abs(eta_electron_pc[iele]) < el_eta_cut && isId_electron_pc[iele]){ // id
		    passedIdEl=true;
		    if (init_electrons[iele]->Pt() > el_pt_cut && abs(init_electrons[iele]->Eta()) < el_eta_cut && isId_electron_pc[iele] && isIso_electron_pc[iele]){ // iso
		      passedIsoEl=true;
		      for(int imuo=0; imuo<init_muons.size(); imuo++){ 
			if (init_muons[imuo]->Pt() > mu_pt_cut){ // pt
			  passedPtMu=true;
			  if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut){ // eta
			    passedEtaMu=true;
			    if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut && isId_muon_pc[imuo]){ // id
			      passedIdMu=true;
			      if (init_muons[imuo]->Pt() > mu_pt_cut && abs(init_muons[imuo]->Eta()) < mu_eta_cut && isId_muon_pc[imuo] && isIso_muon_pc[imuo]){ // iso
				passedIsoMu=true;
				if (debug) cout << "passed iso" << endl;
				if(selectedElectrons.size() == 1 ){ // extra electron veto
				  passedExtraElVeto=true;
				  if(selectedMuons.size() == 1){ // extra muon veto
				    passedExtraMuVeto=true;
				    if (debug) cout << "extra veto passed" << endl;
				    if(abs(d0BeamSpot_electron[1]) < el_d0_cut){ // el d0 blinding  // faco change index
				      passedBlindingEl=true;
				      if (debug) cout << "imuo is " << imuo << endl;
				      if (abs(d0BeamSpot_muon[0]) < mu_d0_cut){ // mu d0 blinding // faco change index
					passedBlindingMu=true;
					if (debug) cout << "passed blinding cut" << endl;
					if(charge_electron[0] * charge_muon[0] == -1){ // to be changed
					  passedElMuOS=true;
					  if(selectedElectrons[0]->DeltaR(*(selectedMuons[0])) > 0.5){ // only one el and mu in selected vector
					    passedElMuNotOverlaping=true;
					    passed++;
					    if (debug) cout << "About to fill the tree!! The number of event that have passed all the cuts is " << passed << endl;
					    myTree->Fill(); 
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
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  //	    */


	  // -------------------------------------------------------------------
	  // e ologic to set boolean that will be used to fill cutflow tables --
	  // -------------------------------------------------------------------
	    

	  // --------------------------
	  // bo filling the cut flow --
	  // --------------------------
	    
	  // Fill full synch cut flow

	  Int_t i_CutFlowTable = 0;
	  CutFlowTable.Fill(d,i_CutFlowTable,1);
	  i_CutFlowTable++;
	  if (passedPtEl){
	    CutFlowTable.Fill(d,i_CutFlowTable,1);
	    i_CutFlowTable++;
	    //	      if (passedEcalCrackVeto){
	    //		CutFlowTable.Fill(d,i_CutFlowTable,1);
	    //		i_CutFlowTable++;
	    if (passedEtaEl){
	      CutFlowTable.Fill(d,i_CutFlowTable,1);
	      i_CutFlowTable++;
	      if (passedIdEl){
		CutFlowTable.Fill(d,i_CutFlowTable,1);
		i_CutFlowTable++;
		if (passedIsoEl){
		  CutFlowTable.Fill(d,i_CutFlowTable,1);
		  i_CutFlowTable++;
		  if (passedPtMu){
		    CutFlowTable.Fill(d,i_CutFlowTable,1);
		    i_CutFlowTable++;
		    if (passedEtaMu){
		      CutFlowTable.Fill(d,i_CutFlowTable,1);
		      i_CutFlowTable++;
		      if (passedIdMu){
			CutFlowTable.Fill(d,i_CutFlowTable,1);
			i_CutFlowTable++;
			if (passedIsoMu){
			  CutFlowTable.Fill(d,i_CutFlowTable,1);
			  i_CutFlowTable++;
			  if (passedExtraElVeto){
			    CutFlowTable.Fill(d,i_CutFlowTable,1);
			    i_CutFlowTable++;
			    if (passedExtraMuVeto){
			      CutFlowTable.Fill(d,i_CutFlowTable,1);
			      i_CutFlowTable++;
			      if (passedBlindingEl){
				CutFlowTable.Fill(d,i_CutFlowTable,1);
				i_CutFlowTable++;
				if (passedBlindingMu){
				  CutFlowTable.Fill(d,i_CutFlowTable,1);
				  i_CutFlowTable++;
				  if (passedElMuOS){
				    CutFlowTable.Fill(d,i_CutFlowTable,1);
				    i_CutFlowTable++;
				    if (passedElMuNotOverlaping){
				      CutFlowTable.Fill(d,i_CutFlowTable,1);
				      i_CutFlowTable++;
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		      //	}
		    }
		  }
		}
	      }
	    }
	  }


	  // --------------------------
	  // eo filling the cut flow --
	  // --------------------------
	    

	  // logic for muons cutflow

	  // first reset all bool to false

	  passedPtMu = false;
	  passedEtaMu = false;
	  passedIdMu = false;
	  passedIsoMu = false;


	  /*
	    for(int imuo=0; imuo<init_muonsTLV.size(); imuo++){ 
	    if (init_muonsTLV[imuo].Pt() > mu_pt_cut){ // pt
	    passedPtMu=true;
	    if (init_muonsTLV[imuo].Pt() > mu_pt_cut && abs(init_muonsTLV[imuo].Eta()) < mu_eta_cut){ // eta
	    passedEtaMu=true;
	    if (KynIdMuons.size()>=1){ // id
	    passedIdMu=true;
	    if (selectedMuons.size()>=1){ // iso
	    passedIsoMu=true;
	    }
	    }
	    }
	    }
	    }
	    
	    // Fill muon cut flow 

	    CutFlow_oneMuTable.Fill(d,0,1);
	    if (passedPtMu){
	    CutFlow_oneMuTable.Fill(d,1,1);
	    if (passedEtaMu){
	    CutFlow_oneMuTable.Fill(d,2,1);
	    if (passedIdMu){
	    CutFlow_oneMuTable.Fill(d,3,1);
	    if (passedIsoMu){
	    CutFlow_oneMuTable.Fill(d,4,1);
	    }
	    }
	    }
	    }
	  */


	  // Fill muon cut flow using continue 
	  // probably does not work as all the statement after this will be skipped which definitely affect the pc_tree. Might consider to put this piece of code in a separate event loop or file 
	    

	    
	  /*
	    cout << "started filling" << endl;
	    CutFlow_oneMuTable.Fill(d,0,1);
	    if (!passedPtMu) continue;
	    CutFlow_oneMuTable.Fill(d,1,1);
	    if (!passedEtaMu) continue;
	    CutFlow_oneMuTable.Fill(d,2,1);
	    if (!passedIdMu) continue;
	    CutFlow_oneMuTable.Fill(d,3,1);
	    if (!passedIsoMu) continue;
	    CutFlow_oneMuTable.Fill(d,4,1);
	    
	    cout << "went through it!!" << endl;
	  */

		
	  // -------------------------------
	  // eo filling the cutflow table --
	  // -------------------------------




	  //////////////////
	  //Filling Trees//
	  //////////////////
	  //	    /*
	  if (debug) cout << "filling the tree, sum of leptons equals to " << nElectrons + nMuons << endl;
	  if (nElectrons_pc >= 1 && nMuons_pc >= 1){
	    myPreCutTree->Fill(); 
	    passed_pc++;
	  }

	  if (testTree && (nElectrons >= 1 || nMuons >= 1 )){
	    //	      cout << "nElectrons is " << nElectrons << endl;
	    //	      cout << "nMuons is " << nMuons <<endl;
	    myTree->Fill();
	    passed++;
	  }


	  Bool_t blindD0_elel = true;
	  Bool_t blindDvz_elel = true;
	  Bool_t blindDv0_elel = true;
	  
	  // el el final state
	  if (nElectrons_elel >= 2){

	    if (isData){

	      if (debug) cout << "Trying to pass blinding condition for data" << endl;
	      
	      // Remove the events were the D0 is too big
	      for (Int_t selel =0; selel <= nElectrons_elel; selel++)
		{
		  if (d0BeamSpot_electron_elel[selel] > 0.2){
		    blindD0_elel=false;
		  }
		}
		 
	      // Remove the events were the DeltaVz is too big
	      for (Int_t selelPairs = 0; selelPairs <= nElectronPairs_elel; selelPairs++)
		{
		  if (deltaVz_elel[selelPairs] > 0.5){
		    blindDvz_elel=false;
		  }			  
		  if (deltaV0_elel[selelPairs] > 0.01){
		    blindDv0_elel=false;
		  }			  
		}
	      if (blindD0_elel && blindDvz_elel && blindDv0_elel){
		myDoubleElTree->Fill();
		passed_elel++;
		if (debug) cout << "Blinding conditions passed!" << endl;
	      }
	      
	    }
	    // if not Data fill it anyway
	    if (!isData) {
	      myDoubleElTree->Fill();
	      passed_elel++;
	    }

	  }
	    

	  Bool_t blindD0_mumu = true;
	  Bool_t blindDvz_mumu = true;
	  Bool_t blindDv0_mumu = true;

	  // mu mu final state
	  if (nMuons_mumu >= 2){	  
	    if (isData){
	      if (debug) cout << "Trying to pass blinding condition for data" << endl;
	      
	      // Remove the events were the D0 is too big
	      for (Int_t selmu =0; selmu <= nMuons_mumu; selmu++)
		{
		  if (d0BeamSpot_muon_mumu[selmu] > 0.2){
		    blindD0_mumu=false;
		  }
		}
		 
	      // Remove the events were the DeltaVz is too big
	      for (Int_t selmuPairs = 0; selmuPairs <= nMuonPairs_mumu; selmuPairs++)
		{
		  if (deltaVz_mumu[selmuPairs] > 0.5){
		    blindDvz_mumu=false;
		  }			  
		  if (deltaV0_mumu[selmuPairs] > 0.01){
		    blindDv0_mumu=false;
		  }			  
		}
	      if (blindD0_mumu && blindDvz_mumu && blindDv0_mumu){
		myDoubleMuTree->Fill();
		passed_mumu++;
		if (debug) cout << "Blinding conditions passed!" << endl;
	      }

	    }
	      
	    // if not Data fill it anyway
	    if (!isData) {
	      myDoubleMuTree->Fill();
	      passed_mumu++;
	    }
	  }


	    
	  if (debug) cout << " DONE filling the tree, sum of leptons equals to " <<nElectrons + nMuons << endl;


	  // filling the cutflow and the histo 1 D

	  // preselection cut
	  histo1D["h_cutFlow"]->Fill(0., 1);
	  CutFlowPreselTable.Fill(d,0,lumiScale);
	  //	    if(!isEBEEGap){
	  if (1){
	    passedEcalCrackVeto = true;
	    histo1D["h_cutFlow"]->Fill(1., 1);
	    CutFlowPreselTable.Fill(d,1,lumiScale);
	    if (selectedElectrons.size() >= 1 ){
	      passedGoodEl=true;
	      histo1D["h_cutFlow"]->Fill(2., 1);
	      CutFlowPreselTable.Fill(d,2,lumiScale);
	      if (selectedMuons.size() >= 1 ){
		passedGoodMu=true;
		histo1D["h_cutFlow"]->Fill(3., 1);
		CutFlowPreselTable.Fill(d,3,lumiScale);
		if (selectedElectrons.size() == 1 ){
		  passedExtraElVeto = true;
		  histo1D["h_cutFlow"]->Fill(4., 1);
		  CutFlowPreselTable.Fill(d,4,lumiScale);
		  if (selectedMuons.size() == 1 ){
		    passedExtraMuVeto = true;
		    histo1D["h_cutFlow"]->Fill(5., 1);
		    CutFlowPreselTable.Fill(d,5,lumiScale);
		    if(abs(selectedElectrons[0]->d0BeamSpot()) < 0.01){
		      passedBlindingEl = true;
		      histo1D["h_cutFlow"]->Fill(6., 1);
		      CutFlowPreselTable.Fill(d,6,lumiScale);
		      if (abs(selectedMuons[0]->d0BeamSpot()) < 0.01){
			passedBlindingMu=true;
			histo1D["h_cutFlow"]->Fill(7., 1);
			CutFlowPreselTable.Fill(d,7,lumiScale);
			//			  if(selectedElectrons[0]->charge() * selectedMuons[0]->charge() == -1){
			if (1){
			  passedElMuOS=true;
			  histo1D["h_cutFlow"]->Fill(8., 1);
			  CutFlowPreselTable.Fill(d,8,lumiScale);
			  Double_t DeltaR = sqrt (2);// to be done
			  if (1){
			    //			  if(selectedElectrons[0]->DeltaR(selectedMuons[0]) > 0.5){
			    passedElMuNotOverlaping=true;
			    histo1D["h_cutFlow"]->Fill(9., 1);
			    CutFlowPreselTable.Fill(d,9,lumiScale);
			    //			      passed++;
			    ///			      if (debug) cout << "About to fill the tree!! The number of event that have passed all the cuts is " << passed << endl;
			    //			      myTree->Fill(); 
			    //			      cout << "Filling Tree!!!" << endl;
			    if (debug) {
			      cout << "npu vtx " << nvtx << endl;
			      cout << "npu is " << npu << endl;
			      cout << "puSF is" << puSF << endl;
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




	  /*
	  // single el 
	  CutFlow_oneElTable.Fill(d,0,lumiScale);
	  if(!isEBEEGap){
	  passedEcalCrackVeto = true;
	  CutFlow_oneElTable.Fill(d,1,lumiScale);
	  if (selectedElectrons.size() >= 1 ){
	  passedGoodEl=true;
	  CutFlow_oneElTable.Fill(d,2,lumiScale);
	  }
	  }


	  // single el 
	  CutFlow_oneMuTable.Fill(d,0,lumiScale);
	  if (selectedMuons.size() >= 1 ){
	  passedGoodMu=true;
	  CutFlow_oneMuTable.Fill(d,1,lumiScale);
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
	  */

	  //cout << "ievt is " << ievt << " and end_d is " << end_d << endl;            
	  //            if( run_num > 10000){//data
	  //	      isdata=1;
	  //            }
            

	    

        } //End Loop on Events

      fout->cd();
      ///Write histogram
      for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
	  cout << "1D Plot: " << it->first << endl;
	  TCanvas *ctemp =  new TCanvas();
	  ctemp->cd();
	  TH1F *temp = it->second;
	  string name = it->first;
	  cout << name << endl; 
	  temp->Draw();
	  delete ctemp;
	  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
	  //c1->SaveAs(pathPNG+"1DPlot/"+name modeString[mode] + "_" + cutLabel + ".png");
	  //	    temp->Write(fout, name, true, pathPNG+"1DPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
	}
	
      if (debug) cout << "Done writing the Tree" << endl;
      fout->Write();   
      fout->Close();
      cout <<"n events having at least two leptons with no id and no iso requirement is  =  "<< passed_pc <<endl;
      cout <<"n events after all the cuts for the e-mu final state is  =  "<< passed <<endl;
      cout <<"n events with at least two id and iso electrons with pt > "+el_pt_cut_str+" GeV is  =  "<< passed_elel <<endl;
      cout <<"n events with at least two id and iso muons with pt > "+mu_pt_cut_str+" GeV is  =  "<< passed_mumu <<endl;
      cout << "Event Count: " << eventCount << endl;
      cout << "Number of event skipped because they failled the trigger requirement is " << skippedEvent << endl;
      cout << "The trigger efficiency is " << eventCount-skippedEvent << "/" << eventCount << " = "  << 100.0*(eventCount-skippedEvent)/eventCount << " % " << endl;
      cout << "Weight Count: " << weightCount << endl;
      //important: free memory
      treeLoader.UnLoadDataset();
    } //End Loop on Datasets

  //  eventlist.close();

  /////////////
  // Writing //
  /////////////

  cout << " - Writing outputs to the files ..." << endl;

  //////////////////////
  // Selection tables //
  //////////////////////

  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool 
  CutFlowPreselTable.TableCalculator(  true, true, true, true, true);
  CutFlowTable.TableCalculator(  true, true, true, true, true);
  CutFlow_oneElTable.TableCalculator(  true, true, true, true, true);
  CutFlow_oneMuTable.TableCalculator(  true, true, true, true, true);

  //    CutFlowExample.Write( filename, WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false) )
  if (strJobNum != "0")
    {
      //	CutFlowPreselTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Presel_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,true,false,false,false);
      CutFlowTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,false,false,false,false);
      CutFlow_oneElTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneEl_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,false,false,false,false);
      CutFlow_oneMuTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneMu_Table"+channelpostfix+"_"+strJobNum+".tex",    false,true,true,false,false,false,false);
    }
  else 
    {
      //	CutFlowPreselTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Presel_Table"+channelpostfix+".tex",    false,true,true,false,false,false,false);
      CutFlowTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"Table"+channelpostfix+".tex",    false,true,true,false,false,false,false);
      CutFlow_oneElTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneEl_Table"+channelpostfix+".tex",    false,true,true,false,false,false,false);
      CutFlow_oneMuTable.Write(  outputDirectory+"/DisplacedTop"+postfix+"OneMu_Table"+channelpostfix+".tex",    false,true,true,false,false,false,false);
    }


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
