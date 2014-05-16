 /*|----------------------------------------------------------------------------------------------------------------------|
  |   Modified Ntuple.cc from Freya intended to be an example analysis macro for new users of the T  th1dir->cd();opBrussels framework  |
  |----------------------------------------------------------------------------------------------------------------------|*/

// ROOT includes
#include "TStyle.h"
#include <cmath>
#include "TNtuple.h"
//Standard C++ library includes (both classes declared within the std namespace)
#include <fstream> //Stream class to both read and write from/to files.
#include <sstream> //Stream class to both read and write from/to strings.
#include <sys/stat.h>

// Include header files to make the interface (user code) available
#include "TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Selection/interface/DisplacedSelection.h"
#include "../TopTreeAnalysisBase/Selection/interface/Selection.h"
#include "../TopTreeProducer/interface/TRootPFJet.h"
#include "../TopTreeProducer/interface/TRootJet.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
//#include "../TopTreeAnalysis/macros/Style.C"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;
//Passing arguments to your program, for example from command line, when a program is invoked
//argc:Nof arguments given at command line (including executable name),argv:array of char*s of size argc
int main (int argc, char *argv[])
{
  //Returns the processor time consumed by the program.
  clock_t start = clock();
     
  //SetStyle if needed, this is defined in the ../TopTreeAnalysis/macros/Style.C
  //  setTDRStyle();

  /*|----------------|
    | Configuration  |
    |----------------|*/
  
  //xml file
  string xmlFileName ="myhiggsconfig.xml";
  //additional xmlfile argument 
    if (argc > 1)
        xmlFileName = (string)argv[1];
    
    const char *xmlfile = xmlFileName.c_str();
    
    std::cout << "1) Using the '" << xmlfile << "' file as input" << std::endl;

    /*|------------------------------|
      | Configuration output format  |
      |------------------------------|*/

    //Create a ROOT TTree with name, title
    TTree *configTree = new TTree("configTree","configuration Tree");
    //Create an array that can initially hold 1000 Dataset objects -> to retrieve the ith object of the array (Dataset*) tcdatasets->At(i)
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    //Create a branch of the configTree(branchname, className=object type, &p_object=pointer to an object)
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    //Create an array that can initially hold 1000 AnalysisEnvironment objects
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    //Create a branch(branchname, className=object type, &p_object=pointer to an object)
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    

    /*|-----------------------|
      |  AnalysisEnvironment  |
      |-----------------------|*/

    //More information in http://w3.iihe.ac.be/~echabert/Documentation/TopTreeAnalysisPackage/CMSSW_35X_2/html/classAnalysisEnvironment.html    
    AnalysisEnvironment anaEnv;
    std::cout << "2) Loading environment ";
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    std::cout << "and creating AnalysisEnvironmentLoader" << std::endl;
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    int verbose = anaEnv.Verbose;
    float Luminosity=1;// = oldLuminosity// 18848.367; //pb^-1??
    float oldLuminosity = anaEnv.Luminosity; //in pb^{-1}

    /*|-----------------|
      |  Load Datasets  |
      |-----------------|*/

    //Loading the branches from the TopTrees: http://w3.iihe.ac.be/~echabert/Documentation/TopTreeAnalysisPackage/CMSSW_35X_2/html/classTTreeLoader.html
    std::cout << "3) Loading datasets ..." << std::endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    treeLoader.LoadDatasets (datasets, xmlfile);

    //Loop over datasets and check for the equivalent luminosity
    //Create booleans to be used throughout the code for data and MC                                                                                                
    bool isData = false;
    bool isMC = false;
    for (unsigned int d=0; d<datasets.size(); d++) {new ((*tcdatasets)[d]) Dataset(*datasets[d]);}

    std::cout << "------------------------------------------------------------------------------------------------" <<  std::endl;
    

    string postfix = "_EventSelection_synclana_noNJets"; // to relabel the names of the output file 
    //nof selected events
    Double_t NEvtsData = 0;
    Int_t *nEvents = new Int_t[datasets.size()];
    Int_t nAllEvents=0;
    Int_t nSkimEvents=0;
    Int_t nRejEvents=0;
    Int_t nLumi270 =0;


    //Output ROOT file                                                                                                                                                        
    string rootFileName ("Checking.root");
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

    //Normal plots
    map<string,TH1F*> histo1D;
    /// MultiSamplePlot                                                                                                                                               
    map<string,MultiSamplePlot*> MSPlot;
    MSPlot["NbOfVertices"] = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
    //MSPlot["MuonPt"]       = new MultiSamplePlot(datasets, "MuonPt", 35, 0, 300, "PT_{#mu}");

    histo1D["scaleFactor"] = new TH1F("scaleFactor","scaleFactor;scaleFactor;#events",100,0,4);

    vector<string> CutFlow;
    CutFlow.push_back(string("Total"));
    CutFlow.push_back(string("Trigger"));
    CutFlow.push_back(string("Skimming"));
    //    CutFlow.push_back(string("$\\geq$1 isScrapingVeto event"));
    CutFlow.push_back(string("$\\geq$1 isGood primary vertex"));
    CutFlow.push_back(string("$\\geq$1 electron with abs(eta) $<$ 2.5"));
    CutFlow.push_back(string("electron ECAL crack veto"));
    CutFlow.push_back(string("$\\geq$1 electron with pt $>$ 25"));
    CutFlow.push_back(string("$\\geq$1 electron with mvaNonTrigHtoZZto4l $>$ 0"));
    CutFlow.push_back(string("electron conversion rejection"));
    CutFlow.push_back(string("$\\geq$1 electron with relPFrhoIso $<$ 0.1"));
    CutFlow.push_back(string("$\\geq$1 muon with abs(eta) $<$ 2.5"));
    CutFlow.push_back(string("$\\geq$1 muon with pt $>$ 25"));
    CutFlow.push_back(string("$\\geq$1 muon with tightIDdisplaced $>$ 0"));
    CutFlow.push_back(string("$\\geq$1 muon with relPFdBetaIso $<$ 0.12"));
    CutFlow.push_back(string("extra electron veto"));
    CutFlow.push_back(string("extra muon veto"));
    CutFlow.push_back(string("$\\geq$1 electron-muon pair with chargeProduct $<$ 0"));
    CutFlow.push_back(string("$\\geq$1 electron-muon pair with deltaR $>$ 0.5"));
    CutFlow.push_back(string("$\\geq$0 jets with pt $>$ 30"));
    CutFlow.push_back(string("electron near jet veto"));
    CutFlow.push_back(string("muon near jet veto"));
    CutFlow.push_back(string("$\\geq$1 electron with abs(correctedD0) $<$ 2"));
    CutFlow.push_back(string("$\\geq$1 muon with abs(correctedD0) $<$ 2"));
    CutFlow.push_back(string("Weighted Yield"));


    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutFlow.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutFlow.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutFlow.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutFlow.push_back(string(LabelNJets));

    if (verbose > 0)
      cout << " - CutsSelectionTable instantiated ..." << endl;
    SelectionTable CutFlowTable(CutFlow, datasets);
    CutFlowTable.SetLuminosity(Luminosity);
    if (verbose > 0)
      cout << " - SelectionTable instantiated ..." << endl;

    
    /*|---------------------|                                                                                                                                         
      |  PileUp Reweighting |
      |---------------------|*/
    
    LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;
       LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");                                                                                  
    LumiWeightsUp = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");                                                                                   
    LumiWeightsDown = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");                                                                                
    std::cout << "*************** Initialized PileUpReWeighting ******************" << std::endl;           
    //LeptonTools* leptonTools = new LeptonTools(false);
    //leptonTools->readMuonSF("LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "LeptonSF/MuonEfficiencies_Run_2012C_53X.root");
    // leptonTools->readElectronSF();                

    /*|----------------------|
      |	 Loop over datasets  | 
      |----------------------|*/

    std::cout << "------------------------------------------------------------------------------------------------" <<  std::endl;
    std::cout << "                Looping over " << datasets.size () << " datasets in total" << std::endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
      if((datasets[d]->Name()).find("Data") == 0 || (datasets[d]->Name()).find("data") == 0 || (datasets[d]->Name()).find("DATA") == 0){
        isData = true;
      }else{
        isMC = true;
      }
      //old and new lumi?? question! what does it mean?
      float oldLuminosity = anaEnv.Luminosity;
      //float Luminosity = oldLuminosity;

      string previousFilename = "";
      int iFile = -1;
      string dataSetName = datasets[d]->Name();
      if(Luminosity > datasets[d]->EquivalentLumi())
	{
          Luminosity = datasets[d]->EquivalentLumi();
        }
      std::cout << "------------------------------------------------------------------------------------------------" <<  std::endl;
      std::cout << "Dataset " << d << ": " << datasets[d]->Name() << " with title " << datasets[d]->Title() << std::endl;
      std::cout <<"xsection: " <<  datasets[d]->Xsection() << ", luminosity: " << datasets[d]->EquivalentLumi() << 
	",  normFactor = " << datasets[d]->NormFactor() << ", containing " << datasets[d]->NofEvtsToRunOver() << " events." << std::endl;
      std::cout << "------------------------------------------------------------------------------------------------" <<  std::endl;
      
      /*|---------------------|
        |  Make the ROOT Tree |
        |---------------------|*/
      
      string roottreename = datasets[d]->Name();
      roottreename+="_tree.root";
      std::cout << "creating tree in file " << roottreename << std::endl;
      TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
      fileout->cd();
      
      Int_t eventID, runID, lumiBlockID;
      Int_t nElectrons;
      //Define an array of doubles with length 10
      Double_t pX_electron[10];
      Double_t pY_electron[10];
      Double_t pZ_electron[10];
      Double_t E_electron[10];
      Double_t pfIso_electron[10];
      Double_t d0_electron[10];
      Int_t charge_electron[10];

      Int_t nMuons;
      Double_t pX_muon[10];
      Double_t pY_muon[10];
      Double_t pZ_muon[10];
      Double_t E_muon[10];
      Double_t pfIso_muon[10];
      Double_t d0_muon[10];
      Int_t charge_muon[10];

      Int_t nJets;
      Double_t pX_jet[10];
      Double_t pY_jet[10];
      Double_t pZ_jet[10];
      Double_t E_jet[10];
      Double_t missingEt;
      
      // Various weights and scale factors
      Double_t pu_weight;
      
      TTree* myTree = new TTree("tree","tree");
      
      myTree->Branch("eventID",&eventID);
      myTree->Branch("runID",&runID);
      myTree->Branch("lumiBlockID",&lumiBlockID);      
      //      myTree->Branch("isdata",&isdata,"isdata/I");
      
      myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
      myTree->Branch("pX_electron",pX_electron,"pX_electron[nElectrons]/D");
      myTree->Branch("pY_electron",pY_electron,"pY_electron[nElectrons]/D");
      myTree->Branch("pZ_electron",pZ_electron,"pZ_electron[nElectrons]/D");
      myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
      myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
      myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
      myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
      
      myTree->Branch("nMuons",&nMuons, "nMuons/I");
      myTree->Branch("pX_muon",pX_muon,"pX_muon[nMuons]/D");
      myTree->Branch("pY_muon",pY_muon,"pY_muon[nMuons]/D");
      myTree->Branch("pZ_muon",pZ_muon,"pZ_muon[nMuons]/D");
      myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
      myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
      myTree->Branch("d0_muon",d0_muon,"d0_electron[nMuons]/D");
      myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
      
      myTree->Branch("nJets",&nJets, "nJets/I");
      myTree->Branch("pX_jet",pX_jet,"pX_jet[nJets]/D");
      myTree->Branch("pY_jet",pY_jet,"pY_jet[nJets]/D");
      myTree->Branch("pZ_jet",pZ_jet,"pZ_jet[nJets]/D");
      myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
      
      myTree->Branch("missingEt",&missingEt,"missingEt/D");
      myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
      
      //myTree->Print();
      
      //open files and load
      std::cout << "4) Loading Event" << std::endl;
      treeLoader.LoadDataset (datasets[d], anaEnv);

      //string Ntupname = "test"+ date_str  +"test" + dataSetName +postfix + ".root";                                                                                        
      string TNtuptitle = "event_counters";
      //TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");                                                                                                            
      //TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"nJets");                                                                                            
      TNtuple * tup = new TNtuple(TNtuptitle.c_str(), TNtuptitle.c_str(),"nAllEvents:nRejEvents:nSkimEvents:nMuons:nElectrons:nJets");

      /*|-------------------------|
	|  Initialize JEC factors |
	|-------------------------|*/
      
      vector<JetCorrectorParameters> vCorrParam;
      /*
      JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
      JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
      JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
      
      //Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
      vCorrParam.push_back(*L1JetPar);
      vCorrParam.push_back(*L2JetPar);
      vCorrParam.push_back(*L3JetPar);
      
            if(isData == true) { // DATA!
	cout << isData << endl;
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
	vCorrParam.push_back(*ResJetCorPar);
	}*/
      
      JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
      // true means redo also the L1
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
        
      bool tightIDdisplacedMuons(TRootMuon* mu);
      bool selectedDisplacedMuons(TRootMuon* mu);
      bool mvaNonTrig_HtoZZto4l(TRootElectron* el);
      double electronRelPFrhoIso(TRootElectron* el, TRootEvent* ev);
      bool jetIDLoose();
      bool isScrapingVeto(TRootEvent* ev);
      nEvents[d] = 0;
      int previousRun = -1;

      std::vector<TRootVertex*> vertex;
      std::vector<TRootMuon*> muons_;
      std::vector<TRootMuon*> muons;
      std::vector <TRootElectron*> electrons_;
      std::vector <TRootElectron*> electrons;
      std::vector <TRootJet*> jets_;
      std::vector <TRootPFJet*> PFJet;
      std::vector <TRootMET* > mets;
      std::vector <TRootGenJet*> genjets;
      
      /*|------------------| 
	|  Loop on events  |                                                                                                  
	|------------------|*/
      
      //      for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++) // event loop
      for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++) // event loop      
	{

          if(ievt%1000 == 0)
            {
	      std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected" << flush <<"\r";
            }
	  
	  /*|-------------|
	    |  Load Event |
	    |-------------|*/
	  TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, muons_, electrons_, jets_, mets);
	  if(isMC==true) {
	    genjets = treeLoader.LoadGenJet(ievt,false);
	    sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
	  }

	  eventID = event->eventId(); //404331002
          runID = event->runId();
          lumiBlockID = event->lumiBlockId();

	  
	  /*|--------------------------------|    
	    |  Determine Event Scale Factor  |
	    |--------------------------------|*/
          
	  // scale factor for the event
	  float scaleFactor = 1.;
	  Double_t lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );          

	  //A simple filter to remove events which arise from the interaction LHC beam with the beampipe known as "beam scraping events".     
	  if(isData==true) {  
	    //if(isScrapingVeto(event)==false) continue;
	  }
	  else {
	    //double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	    //double lumiWeightOLD=lumiWeight;
	  }
	    // PU REWEIGHTING
	  //	    if(isMC==true){
	      lumiWeight=1;
	      //scaleFactor = scaleFactor*lumiWeight;

	  // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
	  // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
          
	  //pu_weight=lumiWeight;          
	  //	  scaleFactor = scaleFactor*lumiWeight;
	  //}
          
	  /*|-----------------------|
            |  Trigger Requirement  |
	    |-----------------------|*/
          
	  string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	  if(previousFilename != currentFilename){
	    previousFilename = currentFilename;
	    iFile++;
	    //std::cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << std::endl;
	  }

          
	  int currentRun = event->runId();          
	  if(previousRun != currentRun) previousRun = currentRun;
	  //cout << "currentRun " << currentRun << " iFile " << iFile;   
	  /////////////////////////////////////////////////////////////////////////////
	  // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
	  /////////////////////////////////////////////////////////////////////////////
          
	  //if(isData==true)	    
	  //jetTools->correctJetJER(jets_, genjets, mets[0], "nominal",false);
	  
	  /*|-----------------|
	    | Event Selection |
	    |-----------------|*/
	  if ((lumiBlockID==270) && (runID==191830)){
	    //                cout << "lumiblock: " << lumiBlockID << ", runnumber: " << runID << ", eventnumber: " << eventID << endl;
	    //std::cout<<" root file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << std::endl;
          }else{
            continue;
          }

	  nAllEvents++;
          nEvents[d]++;
	  //Declare selection instance (See definitions in ../TopTreeAnalysisBase/Selection/src/DisplacedSelection.cc)
	  //DisplacedSelection displacedselection(jets_, muons_, electrons_, event->kt6PFJets_rho());

	  Int_t trigger_v1, trigger_v2, trigger_v3, trigger_v4, trigger_v5, trigger_v6, trigger_v7;
          Bool_t trigged = false;
	  trigger_v3 = treeLoader.iTrigger(string("HLT_Mu22_Photon22_CaloIdL_v3"),currentRun, iFile);//returns the integer corresponding to the trigger name
	  trigger_v4 = treeLoader.iTrigger(string("HLT_Mu22_Photon22_CaloIdL_v4"),currentRun, iFile);
	  trigger_v5 = treeLoader.iTrigger(string("HLT_Mu22_Photon22_CaloIdL_v5"),currentRun, iFile);
	  trigger_v6 = treeLoader.iTrigger(string("HLT_Mu22_Photon22_CaloIdL_v6"),currentRun, iFile);
	  trigger_v7 = treeLoader.iTrigger(string("HLT_Mu22_Photon22_CaloIdL_v7"),currentRun, iFile);//HLT_Mu22_Photon22_CaloIdL_v
	  //EventTrigged -> returns true if event is trigged
	  if(treeLoader.EventTrigged (trigger_v3) || treeLoader.EventTrigged (trigger_v4) || treeLoader.EventTrigged (trigger_v5)|| treeLoader.EventTrigged (trigger_v6) || treeLoader.EventTrigged (trigger_v7)){trigged=true;}

	  Selection displacedselection(jets_, muons_, electrons_, mets, event->kt6PFJets_rho());
	  Bool_t isGoodPV = displacedselection.isPVSelected(vertex, 4, 24, 2.);
	  Bool_t skim_electrons = false;
          Bool_t skim_muons = false;	  
	  
	  if (isScrapingVeto(event) == false)
            {
              //cout << " WARNING: scraping veto failed, skipping this event:  " << isScrapingVeto(event) << endl;
              //continue;
            }
	  if (trigged == false) 
	    {
	      //cout << " WARNING: specified trigger was not found in this event, skipping event:  " << trigged << endl;
	      //continue;
	    }

          if(isGoodPV == false) 
	    { 
	      //cout << " WARNING: no good primary vertex found, skipping event:  " << isGoodPV << endl; 
	      //continue;
	    }
	  
	  electrons.clear(); // always remember to clear your vector for each event!  
	  for(unsigned int i=0;i<electrons_.size();i++){
	    if(trigged==true){
	      cout << eventID << " ---------------------------> electrons_ pt,eta " << electrons_[i]->Pt() << " " << abs(electrons_[i]->Eta()) << endl;
	      if (fabs(electrons_[i]->Eta())<3 && electrons_[i]->Pt()>20)
		{
		  electrons.push_back(electrons_[i]);
		  skim_electrons=true;
		}
	    }
	  }
	
	  muons.clear();	  
          for(unsigned int i=0;i<muons_.size();i++){
	    if(trigged==true){
	      cout << eventID << " ---------------------------> muons_ pt,eta " << muons_[i]->Pt() << " " << abs(muons_[i]->Eta()) << endl;
              if (fabs(muons_[i]->Eta())<3 && muons_[i]->Pt()>20)
                {
                  muons.push_back(muons_[i]);
                  skim_muons=true;
                }
            }	  
	  }

	  Bool_t skimming = (isScrapingVeto(event) && trigged && isGoodPV && skim_electrons && skim_muons);
          if(skimming == true) {	    
	    nSkimEvents++;
	  }else{
	    nRejEvents++;
	  }

	  histo1D["scaleFactor"]->Fill(scaleFactor);
	  MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
	  //MSPlot["MuonPt"]->Fill(muons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	  //these loops over the objects are anw done in the Selection.cc and are not needed but                                                   
	  //for now it's better to be able to use them here for our analysis and our selections                                               
          //Electron Selection                                                                                                                                                      
	  Bool_t eta_electrons=false;
	  Bool_t crack_electrons=false;
	  Bool_t pt_electrons=false;
	  Bool_t mva_electrons=false;
	  Bool_t conversion_electrons=false;
	  Bool_t iso_electrons=false;

	  std::vector<TRootElectron*> displacedelectrons;
          for(unsigned int i=0;i<electrons.size();i++) {
	    if (eventID==403622270 || eventID==403685107 || eventID==403838216 || eventID==4038479 || eventID==404350395 || eventID==403867523 || eventID==404907940
		){cout << "initial electron vector: "<<electrons_.size() << "skimmed electron vector "<< electrons.size() << endl;
	      cout << eventID << " ---------------------------> electron pt, eta" << electrons[i]->Pt() << " " << abs(electrons[i]->Eta()) << endl;
	    }	    
            if (abs(electrons[i]->Eta()) <2.5) {
	      eta_electrons=true;
              if(fabs(electrons[i]->superClusterEta()) < 1.4442 || fabs(electrons[i]->superClusterEta()) > 1.5660){
		crack_electrons=true;
		if(electrons[i]->Pt()>25){
		  pt_electrons=true;
                  if(mvaNonTrig_HtoZZto4l(electrons[i])) {
		    mva_electrons=true;
                    if(electrons[i]->passConversion() && electrons[i]->missingHits()==0){
		      conversion_electrons=true;
                      if ((electronRelPFrhoIso(electrons[i], event))<0.1){
			iso_electrons=true;
                        displacedelectrons.push_back(electrons[i]);
                      }
                    }
                  }
                }
              }
            }
          }
	  cout << "eta_electrons bool: " << eta_electrons<< endl;
	  //Muon Selection                                                                                                                                                   
	  Bool_t eta_muons=false;
	  Bool_t pt_muons=false;
	  Bool_t tightID_muons=false;
	  Bool_t iso_muons=false;
	  std::vector<TRootMuon*> displacedmuons;
	  //if (displacedelectrons.size()>0){
          for(unsigned int i=0;i<muons.size();i++) {
	    if (eventID==403622270 || eventID==403685107 || eventID==403838216 || eventID==4038479 || eventID==404350395 || eventID==403867523){
              cout << eventID << " ---------------------------> muon pt, eta " << muons[i]->Pt() << " " << abs(muons[i]->Eta()) << endl;
            }
	    if (fabs(muons[i]->Eta())<2.5){
	      eta_muons=true;
              if (muons[i]->Pt()>25){
		pt_muons=true;
                if (tightIDdisplacedMuons(muons[i])){
		  tightID_muons=true;
                  if(selectedDisplacedMuons(muons[i])){
		    iso_muons=true;
                    displacedmuons.push_back(muons[i]);
                  }
                }
              }
            }
          }
	  //}
	  //Jet Selection                                                                                                                                             
	  for(unsigned int i=0;i<jets_.size();i++) {
            TRootJet* pfjets_ = (TRootJet*) jets_[i];
            if(pfjets_->jetType() == 2){
              const TRootPFJet* PFJet = static_cast<const TRootPFJet*>(pfjets_);
              if(jetIDLoose() == true){
                if (PFJet->Pt() >25 && fabs(PFJet->Eta()) < 2.4 ){
                }
              }
            }
          }

          vector<TRootJet*> selectedJets= displacedselection.GetSelectedJets(25, 2.4, true);
	  if (displacedelectrons.size()>0 && displacedmuons.size()>0){cout << "I have both displaced objects" << endl;}


	  if ((lumiBlockID==270) && (runID==191830)){
	    nLumi270++;
	    CutFlowTable.Fill(d,0,scaleFactor*lumiWeight);
	    //RUN LUMIBLOCK EVENTNUMBER TRIGGERED
	    cout << runID << " " << lumiBlockID << " " << eventID << " " << trigged << endl;
	    //cout << "nAllEvents " << nAllEvents << ", nRejEvents " << nRejEvents << ", nSkimEvents " << nSkimEvents << ", nMuons " << nMuons << ", nElectrons " << nElectrons << ", NofEvtsToRunOver "<<datasets[d]->NofEvtsToRunOver() << endl;
	    //cout << isScrapingVeto(event) << trigged << isGoodPV << skim_electrons << skim_muons << endl;
	  if(trigged==true) {
	    //if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJets;
	    CutFlowTable.Fill(d,1,scaleFactor*lumiWeight);
	    //if(skimming==true){
	    if(isScrapingVeto(event)==true){
	      CutFlowTable.Fill(d,2,scaleFactor*lumiWeight);
	      if (isGoodPV==true) {
		//cout << isScrapingVeto(event) << trigged << isGoodPV << skim_electrons << skim_muons << endl;
		CutFlowTable.Fill(d,3,scaleFactor*lumiWeight);
		if(eta_electrons==true) {
		  CutFlowTable.Fill(d,4,scaleFactor*lumiWeight);
		  if(crack_electrons==true) {
		    CutFlowTable.Fill(d,5,scaleFactor*lumiWeight);
		    if (pt_electrons==true) {
		      CutFlowTable.Fill(d,6,scaleFactor*lumiWeight);
		      CutFlowTable.Fill(d,7,scaleFactor*lumiWeight);
		      if(mva_electrons==true) {
			CutFlowTable.Fill(d,8,scaleFactor*lumiWeight);
			if(conversion_electrons==true) {
			  CutFlowTable.Fill(d,9,scaleFactor*lumiWeight);
			  if(iso_electrons==true) {
			    CutFlowTable.Fill(d,10,scaleFactor*lumiWeight);
			    if(displacedelectrons.size()>0) {
			      CutFlowTable.Fill(d,11,scaleFactor*lumiWeight);
			      if(eta_muons==true){
				CutFlowTable.Fill(d,12,scaleFactor*lumiWeight);		
				if(pt_muons==true){
				  CutFlowTable.Fill(d,13,scaleFactor*lumiWeight);
				  if(tightID_muons==true){
				    CutFlowTable.Fill(d,14,scaleFactor*lumiWeight);
				    if(iso_muons==true){
				      CutFlowTable.Fill(d,15,scaleFactor*lumiWeight);
				      if(displacedelectrons.size()==1){
					CutFlowTable.Fill(d,16,scaleFactor*lumiWeight);
					if(displacedmuons.size()==1){
					  CutFlowTable.Fill(d,17,scaleFactor*lumiWeight);
					  //			  cout << displacedelectrons[0]->charge() << displacedmuons[0]->charge()<< endl;
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
	  // cout << nLumi270 << endl;
	  

	  //Declare selection instance
	  //See definitions in ../TopTreeAnalysisBase/Selection/src/DisplacedSelection.cc
	  //DisplacedSelection displacedselection(jets_, displacedmuons, electrons_, event->kt6PFJets_rho());
	  //displacedselection.setDisplacedMuonCuts(25,2.5,0.4,0.3,1,5,0); // standard mu selection but with looser iso
	  //displacedselection.setDisplacedElectronCuts(10,2.5,0.4,0.5,0.3,0); // standard ele selection but with looser iso
	  //displacedselection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); // standard TOP jet selection

	  /*vector<TRootMuon*> selectedMuons = displacedselection.GetSelectedMuons(25, 2.5, 0.12);
	  vector<TRootElectron*> selectedElectrons = displacedselection.GetSelectedElectrons(25, 2.5, 0.1, 0);


	  nElectrons=0;
	  for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
	    pX_electron[nElectrons]=selectedElectrons[iele]->Px();
	    pY_electron[nElectrons]=selectedElectrons[iele]->Py();
	    pZ_electron[nElectrons]=selectedElectrons[iele]->Pz();
	    E_electron[nElectrons]=selectedElectrons[iele]->E();
	    d0_electron[nElectrons]=selectedElectrons[iele]->d0();
	    nElectrons++;
	  }
	  
	  nMuons=0;
	  for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
	    pX_muon[nMuons]=selectedMuons[imuo]->Px();
	    pY_muon[nMuons]=selectedMuons[imuo]->Py();
	    pZ_muon[nMuons]=selectedMuons[imuo]->Pz();
	    E_muon[nMuons]=selectedMuons[imuo]->E();
	    d0_muon[nMuons]=selectedMuons[imuo]->d0();
	    charge_muon[nMuons]=selectedMuons[imuo]->charge();
	    nMuons++;
	  }

	  nJets=0;
	  for(int ijet=0; ijet<selectedJets.size() && nJets<10; ijet++){
	    pX_jet[nJets]=selectedJets[ijet]->Px();
	    pY_jet[nJets]=selectedJets[ijet]->Py();
	    pZ_jet[nJets]=selectedJets[ijet]->Pz();
	    E_jet[nJets]=selectedJets[ijet]->E();
	    nJets++;
	  }
	  */
	  missingEt=mets[0]->Pt();







	  if(nElectrons+nMuons>0) {
	    //testing super-lite TNtuples (to be compared with ttree)                                                                                                           
	    myTree->Fill();
	    //std::cout << "found " << nMuons << " muons and " << nElectrons << " electrons!" << std::endl;
	  }
	  tup->Fill(nAllEvents, nRejEvents, nSkimEvents, nMuons, nElectrons, nJets);          
	}//end of loop on event
      std::cout<<std::endl;
              
      /*|--------------|
	| Jet Cleaning |
	|--------------|*/
      
      if (jecUnc) delete jecUnc;
      if (jetTools) delete jetTools;
      
      tup->Write();
      //      tupfile->Close();

      myTree->Write();
      fileout->Write();
      fileout->Close();
      // delete myTree;
      delete fileout;
      
      //important: free memory
      treeLoader.UnLoadDataset();
      
	}//end of loop on datasets

    //Selection tables      
    //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)                                                                                        
    CutFlowTable.TableCalculator(false, true, true, true, true);
string selectiontable = "SelectionTable_table.tex";
//O  string rootFileName ("FourTop"+postfix+channelpostfix+".root");ptions : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false) 
 CutFlowTable.Write(selectiontable.c_str(), true,true,true,true,false,true,true);

    if (verbose > 0)    
    //Once everything is filled ...
      std::cout << " We ran over all the data" << std::endl;
    
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    //    if (verbose > 0)
    //        std::cout << "Treating the special plots." << std::endl;
    
    /////////////////////////                                                                                                                        
    // Write out the plots //                                                                                                         
    /////////////////////////                                                                                                                                             
    
    mkdir("TestingPlots",0777);
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
	string name = it->first;
	MultiSamplePlot *temp = it->second;
	temp->Draw(name,0,false,false,false,5); //label, RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal                       
	temp->Write(fout, name, false, "TestingPlots/");//file, label, savePNG, pathPNG, "png"
      }
    /*
    TDirectory* th1dir = fout->mkdir("Histos1D");
    th1dir->cd();
    foxr(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
      {
	string name = it->first;
        TH1F *temp = it->second;
	temp->Write();
      }
    

    */
    delete fout;
    delete tcdatasets;
    delete tcAnaEnv;
    //   delete configTree;
    
    std::cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << std::endl;
    
    std::cout << "********************************************" << std::endl;
    std::cout << "           End of the program !!            " << std::endl;
    std::cout << "********************************************" << std::endl;
    
    return 0;
}

//Scraping veto                                                                                                                                                        
bool isScrapingVeto(TRootEvent* ev){
  bool Yes = false;
  Yes = (ev->nTracks() > 10 && (((float) ev->nHighPurityTracks()) / ((float) ev->nTracks() ) > 0.25));
  return Yes;
}

bool tightIDdisplacedMuons(TRootMuon* mu){
  bool Yes=false;
  Yes=(
       mu->isGlobalMuon() 
       && mu->isPFMuon() 
       && mu->chi2()<10 
       && mu->nofValidPixelHits()>0 
       && mu->nofMatchedStations()>1 
       && mu->nofValidPixelHits()>0 
       && mu->nofTrackerLayersWithMeasurement()>5
       );
  return Yes;
}

bool selectedDisplacedMuons(TRootMuon* mu){
  bool Yes=false;
  float muonRelPFdBetaIso=(mu->chargedHadronIso(4)+max(0.0,mu->neutralHadronIso(4)+mu->photonIso(4)-0.5*mu->puChargedHadronIso(4)))/mu->Pt(); //dBeta corrected
  Yes=(
       tightIDdisplacedMuons(mu) && muonRelPFdBetaIso<0.12
       );
  return Yes;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Non_triggering_MVA
bool mvaNonTrig_HtoZZto4l(TRootElectron* el){
  bool Yes=false;
  if(fabs(el->superClusterEta())<0.8 && el->mvaNonTrigId()>-0.34)
    Yes = true;
  else if(fabs(el->superClusterEta())<1.479 && el->mvaNonTrigId()>-0.65)
    Yes = true;
  else if(fabs(el->superClusterEta())<2.5 && el->mvaNonTrigId()>0.60)
    Yes = true;
    return Yes;
}

double electronRelPFrhoIso(TRootElectron* el, TRootEvent* ev) {
   double EffectiveArea=0.;
   double isoCorr=0;
   double electronRelPFrhoIso=0;
   //HCP2012 for el conesize = 0.3 http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
   if(fabs(el->superClusterEta()) >= 0.0   && fabs(el->superClusterEta()) < 1.0   ) EffectiveArea = 0.130;
   if(fabs(el->superClusterEta()) >= 1.0   && fabs(el->superClusterEta()) < 1.479 ) EffectiveArea = 0.137;
   if(fabs(el->superClusterEta()) >= 1.479 && fabs(el->superClusterEta()) < 2.0   ) EffectiveArea = 0.067;
   if(fabs(el->superClusterEta()) >= 2.0   && fabs(el->superClusterEta()) < 2.2   ) EffectiveArea = 0.089;
   if(fabs(el->superClusterEta()) >= 2.2   && fabs(el->superClusterEta()) < 2.3   ) EffectiveArea = 0.107;
   if(fabs(el->superClusterEta()) >= 2.3   && fabs(el->superClusterEta()) < 2.4   ) EffectiveArea = 0.110;
   if(fabs(el->superClusterEta()) >= 2.4) EffectiveArea = 0.138;
   isoCorr = ev->kt6PFJets_rho()*EffectiveArea;
   electronRelPFrhoIso=(el->chargedHadronIso() + max(0.0, el->neutralHadronIso() + el->photonIso()-isoCorr))/ el->Pt();
   return electronRelPFrhoIso;
}

bool jetIDLoose(){
  bool Yes=false;
  TRootPFJet* jet;
  Yes=(
       (jet->nConstituents() > 1)
       && (jet->neutralHadronEnergyFraction() < 0.99)
       && (jet->neutralEmEnergyFraction() < 0.99)
       && (fabs(jet->Eta()) >= 2.4 || jet->chargedEmEnergyFraction() < 0.99)
       && (fabs(jet->Eta()) >= 2.4 || jet->chargedHadronEnergyFraction() > 0)
       && (fabs(jet->Eta()) >= 2.4 || jet->chargedMultiplicity() > 0)
       );
  return Yes;
}
