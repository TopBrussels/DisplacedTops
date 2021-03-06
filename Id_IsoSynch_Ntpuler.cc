///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Id_IsoSynch_Ntpuler.cc: This macro is intended to be used to make a synchronisation exercise with people from Ohio  /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{
    
    
    clock_t start = clock();
    
    
    /////////////////////
    // Configuration
    /////////////////////
    
    //xml file
    //    string xmlFileName ="myhiggsconfig.xml";
    string xmlFileName ="testconfig.xml";

    if (argc > 1)
        xmlFileName = (string)argv[1];
    
    const char *xmlfile = xmlFileName.c_str();
    
    cout << "********************************************************" << endl;
    cout << "used config file: " << xmlfile << endl;
    cout << "********************************************************" << endl;
    
    

    cout << "********************************************************" << endl;
    cout<<"creating datasets ..."<<endl;
    cout << "********************************************************" << endl;
    
    //Configuration output format
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
    ////////////////////////////////////
    /// AnalysisEnvironment
    ////////////////////////////////////
    
    AnalysisEnvironment anaEnv;
    cout << "********************************************************" << endl;
    cout<<"Loading environment ..."<<endl;
    cout << "********************************************************" << endl;
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    
    cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
    
    int verbose = 2;//anaEnv.Verbose;

    
    cout << "now done creating AnalysisEnvironmentLoader" << endl;
    cout << "********************************************************" << endl;
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    verbose = anaEnv.Verbose;
    float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
    
    Double_t muoneffaverage[2]={0,0};
    Double_t jeteffaverage[2]={0,0};
    Double_t eleeffaverage[2]={0,0};
    
    TTreeLoader treeLoader;
    //    cout << " - Load datasets ..." << endl;
    vector < Dataset* > datasets;
    
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout << "now loaded " << datasets.size() << " datasets" << endl;
    for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    
    float Luminosity = oldLuminosity;
    
    cout << "********************************************************" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
        
        string dataSetName = datasets[d]->Name();
        cout << "datasets: " << dataSetName << endl;
    }
    cout << "********************************************************" << endl;
    

     
     // start a new table
     vector<string> CutFlowElectron;
     CutFlowElectron.push_back(string("Total"));
     CutFlowElectron.push_back(string("At least one electron with pt $>$ 25"));
     CutFlowElectron.push_back(string("At least one electron with abs(eta) $<$ 2.5"));
     CutFlowElectron.push_back(string("At least one electron with iso $<$ 0.10 "));
     CutFlowElectron.push_back(string("At least one electron with pt $>$ 30  "));
     CutFlowElectron.push_back(string("At least one electron with abs(etaSC) $<$ 1.479"));
     CutFlowElectron.push_back(string("At least one electron with abs(deltaEtaIn) $<$ 0.005"));
     CutFlowElectron.push_back(string("At least one electron with abs(deltaPhiIn) $<$ 0.02 "));
     CutFlowElectron.push_back(string("At least one electron with hadronicOverEm $<$ 0.03"));
     CutFlowElectron.push_back(string("At least one electron with d0 $<$ 0.05"));
     CutFlowElectron.push_back(string("At least one electron with dz $<$ 0.2 "));
     CutFlowElectron.push_back(string("At least one electron with 1/E - 1/P $<$ 0.15"));
     CutFlowElectron.push_back(string("At least one electron with passConversion"));
     CutFlowElectron.push_back(string("At least one electron with missingHits $<=$ 1)"));

     SelectionTable CutFlowTableElectron(CutFlowElectron, datasets);
     CutFlowTableElectron.SetLuminosity(Luminosity);

     // start a table
     vector<string> CutFlowMuon;
     CutFlowMuon.push_back(string("Total"));
     CutFlowMuon.push_back(string("At least one muon with pt $>$ 25"));
     CutFlowMuon.push_back(string("At least one muon with abs(eta)$<$ 2.5"));
     CutFlowMuon.push_back(string("At least one muon with Iso$<$ 0.12"));
     CutFlowMuon.push_back(string("At least one muon with GlobalId "));
     CutFlowMuon.push_back(string("At least one muon with pt $>$ 26 "));
     CutFlowMuon.push_back(string("At least one muon with abs(eta) $<$ 2.1 "));
     CutFlowMuon.push_back(string("At least one muon with chi2 $<$ 10"));
     CutFlowMuon.push_back(string("At least one muon with nofTrackerLayersWithMeasurement() $>$ 5"));
     CutFlowMuon.push_back(string("At least one muon with nofValidMuHits() $>$ 0"));
     CutFlowMuon.push_back(string("At least one muon with d0 $<$ 0.2"));
     CutFlowMuon.push_back(string("At least one muon with dz $<$ 0.5"));
     CutFlowMuon.push_back(string("At least one muon with nofValidPixelHits $>$ 0 "));
     CutFlowMuon.push_back(string("At least one muon with nofMatchedStations $>$ 1 "));

     SelectionTable CutFlowTableMuon(CutFlowMuon, datasets);
     CutFlowTableMuon.SetLuminosity(Luminosity);


     if (verbose > 0){
       cout << " - CutsSelectionTable instantiated ..." << endl;
       cout << " - SelectionTable instantiated ..." << endl;
     }
	 
     

    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    ////////////////////////////////////
    //	Loop on datasets
    ////////////////////////////////////
    
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
        if (verbose > 1)
            std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
//         make root tree file name
        string roottreename = datasets[d]->Name();
        roottreename+="_tree.root";
	//        cout << "creating tree in file " << roottreename << endl;
        
        // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        //////////////////////////////
        // My tree - variables      //
        //////////////////////////////
	Int_t nElectronsPostCut;
        Int_t nElectrons;
        Double_t pX_electron[10];
        Double_t pY_electron[10];
        Double_t pZ_electron[10];
        Double_t E_electron[10];
        Double_t d0_electron[10];
	Double_t chargedHadronIso_electron[10];
	Double_t neutralHadronIso_electron[10];
	Double_t photonIso_electron[10];
        Double_t pfIso_electron[10];
        Int_t charge_electron[10];
        
        Int_t nMuonsPostCut;
        Int_t nMuons;
        Double_t pX_muon[10];
        Double_t pY_muon[10];
        Double_t pZ_muon[10];
        Double_t E_muon[10];
        Double_t d0_muon[10];
	Double_t chargedHadronIso_muon[10];
	Double_t neutralHadronIso_muon[10];
	Double_t photonIso_muon[10];
        Double_t pfIso_muon[10];
        Int_t charge_muon[10];
        
        Int_t nJets;
        Double_t pX_jet[10];
        Double_t pY_jet[10];
        Double_t pZ_jet[10];
        Double_t E_jet[10];
        Double_t dx_jet[10];
        Double_t dy_jet[10];
        Double_t missingEt;
        Int_t isdata;
        // various weights
        Double_t pu_weight;
        
        // define the output tree
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");

        
	//        myTree->Branch("nElectronsPostCut",&nElectronsPostCut, "nElectronsPostCut/I");//
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//
        myTree->Branch("pX_electron",pX_electron,"pX_electron[nElectrons]/D");
        myTree->Branch("pY_electron",pY_electron,"pY_electron[nElectrons]/D");
        myTree->Branch("pZ_electron",pZ_electron,"pZ_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
	myTree->Branch("chargedHadronIsoelectron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
	myTree->Branch("neutralHadronIsoelectron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
	myTree->Branch("photonIsoelectron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        myTree->Branch("pfIsoelectron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");

	//	myTree->Branch("nMuonsPostCut",&nMuonsPostCut, "nMuonsPostCut/I");
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pX_muon",pX_muon,"pX_muon[nMuons]/D");
        myTree->Branch("pY_muon",pY_muon,"pY_muon[nMuons]/D");
        myTree->Branch("pZ_muon",pZ_muon,"pZ_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
	myTree->Branch("chargedHadronIsomuon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
	myTree->Branch("neutralHadronIsomuon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
	myTree->Branch("photonIsomuon",photonIso_muon,"photonIso_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        
        myTree->Branch("nJets",&nJets, "nJets/I");
        myTree->Branch("pX_jet",pX_jet,"pX_jet[nJets]/D");
        myTree->Branch("pY_jet",pY_jet,"pY_jet[nJets]/D");
        myTree->Branch("pZ_jet",pZ_jet,"pZ_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("dx_jet",dx_jet,"dx_jet[nJets]/D");
        myTree->Branch("dy_jet",dy_jet,"dy_jet[nJets]/D");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
        
	//        myTree->Print();

        
        
        
        //open files and load
        treeLoader.LoadDataset (datasets[d], anaEnv);
        
        
        //////////////////////////////////////////////////
        /// Initialize Jet energy correction factors   ///
        //////////////////////////////////////////////////
   	   
        vector<JetCorrectorParameters> vCorrParam;
        
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
        
        // true means redo also the L1 corrections (see CMS documentation to learn what this means)
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
        
        
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
        
        // some bookkeeping variables
        nEvents[d] = 0;
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
        
        // some printout
        cout << "running over " << datasets[d]->NofEvtsToRunOver() << endl;
        
	//int N_loop = datasets[d]->NofEvtsToRunOver();
	int N_loop = 10;
	
	
        // start event loop
	for (unsigned int ievt = 0; ievt < N_loop; ievt++) // event loop
        {
            
            // the objects loaded in each event
            vector < TRootVertex* > vertex;
            vector < TRootMuon* > init_muons;
            vector < TRootMuon* > postCut_muons;
	    vector < TRootElectron* > init_electrons;
	    vector < TRootElectron* > postCut_electrons;
            vector < TRootJet* > init_jets_corrected;
            vector < TRootJet* > init_jets;
            vector < TRootMET* > mets;
            vector < TRootGenJet* > genjets;
            
            nEvents[d]++;
            
            if(ievt%1000 == 0)
                std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected" << flush<<"\r";
            
            ////////////////
            // LOAD EVENT //
            ////////////////
            
            TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
            
            
            
            // determine if this is data from the data set name (watch out)
            isdata=0;
            if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
                genjets = treeLoader.LoadGenJet(ievt,false);
                //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
            }
            else{
                isdata=1;
            }
            
            
            /////////////////////////////////
            // DETERMINE EVENT SCALEFACTOR //
            /////////////////////////////////
            
            // scale factor for the event
            float scaleFactor = 1.;
            





            // PU reweighting
            
            double lumiWeight = 1 ; //LumiWeights.ITweight( (int)event->nTruePU() ); // currently no pile-up reweigting applied

            if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
                lumiWeight=1;
            
    // filled into output file
            pu_weight=lumiWeight;
            
            scaleFactor = scaleFactor*lumiWeight;
            
            /////////////////////////////////
            // print when you change file  //
            /////////////////////////////////
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename){
                previousFilename = currentFilename;
                iFile++;
                cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << endl;
            }
            
            // get run number
            int currentRun = event->runId();
            
            if(previousRun != currentRun)
                previousRun = currentRun;
            
            
            /////////////////////////////////////////////////////////////////////////////
            // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
            /////////////////////////////////////////////////////////////////////////////
            // not applied during CSA14/PHYS14
            //            if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
            //
            //                jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
            
            /////////////////////
            // EVENT SELECTION //
            /////////////////////
            
            //Declare selection instance
            Selection selection(init_jets_corrected, init_muons, init_electrons, mets);

            // the default selection is fine for normal use - if you want a special selection you can use the functions here
            //selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); //  void setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon);
	    //            selection.setMuonCuts(25,2.5,0.12,0.2,0.3,1,0.5,5,0); // void setMuonCuts(float Pt, float Eta, float RelIso, float d0, float DRJets, int NMatchedStations, float Dz, int NTrackerLayersWithMeas, int NValidPixelHits);
	    selection.setLooseMuonCuts();
	    //            selection.setElectronCuts(25,2.5,0.1,0.02,0.,0.3,0); // void setElectronCuts(float Pt, float Eta, float RelIso, float d0, float MVAId, float DRJets, int MaxMissingHits);
	    selection.setLooseElectronCuts();

	    //            faco
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
	    //            if(!isGoodPV)
	    //                continue;
            
            missingEt=mets[0]->Pt();

            // get the 'good' objects from the selection object
            vector<TRootJet*> selectedJets= selection.GetSelectedJets(true);
            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(selectedJets);
            
            // bookkeeping
            eleeffaverage[0]+=init_electrons.size()*scaleFactor;
            eleeffaverage[1]+=selectedElectrons.size()*scaleFactor;
            muoneffaverage[0]+=init_muons.size()*scaleFactor;
            muoneffaverage[1]+=selectedMuons.size()*scaleFactor;
            jeteffaverage[0]+=init_jets_corrected.size()*scaleFactor;
            jeteffaverage[1]+=selectedJets.size()*scaleFactor;
            
	    /*            
            // loop over electrons
            nElectrons=0;
            for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
                pX_electron[nElectrons]=selectedElectrons[iele]->Px();
                pY_electron[nElectrons]=selectedElectrons[iele]->Py();
                pZ_electron[nElectrons]=selectedElectrons[iele]->Pz();
                E_electron[nElectrons]=selectedElectrons[iele]->E();
                d0_electron[nElectrons]=selectedElectrons[iele]->d0();
		chargedHadronIso_electron[nElectrons]=selectedElectrons[iele]->chargedHadronIso(3);
		neutralHadronIso_electron[nElectrons]=selectedElectrons[iele]->neutralHadronIso(3);
		photonIso_electron[nElectrons]=selectedElectrons[iele]->photonIso(3);
                pfIso_electron[nElectrons]=selectedElectrons[iele]->relPfIso(3,0);
                charge_electron[nElectrons]=selectedElectrons[iele]->charge();
                nElectrons++;
            }
            // loop over muons
            nMuons=0;
            for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
                pX_muon[nMuons]=selectedMuons[imuo]->Px();
                pY_muon[nMuons]=selectedMuons[imuo]->Py();
                pZ_muon[nMuons]=selectedMuons[imuo]->Pz();
                E_muon[nMuons]=selectedMuons[imuo]->E();
                d0_muon[nMuons]=selectedMuons[imuo]->d0();
		chargedHadronIso_muon[nMuons]=selectedMuons[imuo]->chargedHadronIso(4);
		neutralHadronIso_muon[nMuons]=selectedMuons[imuo]->neutralHadronIso(4);
		photonIso_muon[nMuons]=selectedMuons[imuo]->photonIso(4);
                pfIso_muon[nMuons]=selectedMuons[imuo]->relPfIso(4,0);
                charge_muon[nMuons]=selectedMuons[imuo]->charge();
                nMuons++;
            }
            // loop over jets
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<10; ijet++){
                pX_jet[nJets]=selectedJets[ijet]->Px();
                pY_jet[nJets]=selectedJets[ijet]->Py();
                pZ_jet[nJets]=selectedJets[ijet]->Pz();
                E_jet[nJets]=selectedJets[ijet]->E();
                dx_jet[nJets]=selectedJets[ijet]->vx();
                dy_jet[nJets]=selectedJets[ijet]->vy();
                nJets++;
            }
	    */

	    Bool_t trigged = true;
	    //	    Bool_T 
	    

	    // put the vector in TLorenzt Vector
	    vector<TLorentzVector> init_muonsTLV, init_electronsTLV, postCut_muonsTLV, postCut_electronsTLV;
	   
	    for(int iele=0; iele<init_electrons.size() && nElectrons<10; iele++)
	      {
		init_electronsTLV.push_back(*init_electrons[iele]);
	      }

	    for(int imuo=0; imuo<init_muons.size() && nMuons<10; imuo++)
	      {
		init_muonsTLV.push_back(*init_muons[imuo]);
	      }	  



	    // Synch cut on electron : pt, eta, Iso, ID                                                                                      
            Bool_t passedPtEl = false;
            Bool_t passedEtaEl = false;
            Bool_t passedIsoEl = false;
            Bool_t passedNewPtEl = false;
            Bool_t passedNewEtaEl = false;
            Bool_t passedDeltaEtaInEl = false;
            Bool_t passedDeltaPhiInEl = false;
            Bool_t passedHadronicOverEmEl = false;
            Bool_t passedD0El = false;
            Bool_t passedDzEl = false;
            Bool_t passedEtoPFractionEl = false;
            Bool_t passedConversionEl= false;
            Bool_t passedMissingHits= false;

	    // electrons
	    for(int iele=0; iele<init_electrons.size(); iele++){
              if (init_electrons[iele]->Pt() > 25){
		passedPtEl = true;
                if (abs(init_electrons[iele]->Eta()) < 2.5){
		  passedEtaEl =true;
		  if (init_electrons[iele]->relPfIso(3, 0.5) < 0.10){
		    passedIsoEl = true;
		    if (init_electrons[iele]->Pt() > 30){
		      passedNewPtEl = true;
		      //		      if (fabs(init_electrons[iele]->superClusterEta()) <= 1.479){
		      if (1){
			passedNewEtaEl = true;
			if (fabs(init_electrons[iele]->deltaEtaIn()) < 0.005){
			  passedDeltaEtaInEl = true;
			  if (fabs(init_electrons[iele]->deltaPhiIn()) < 0.02){
			    passedDeltaPhiInEl = true;
			    if (init_electrons[iele]->hadronicOverEm() < 0.03){
			      passedHadronicOverEmEl = true;
			      if (init_electrons[iele]->d0() < 0.05){
				passedD0El = true;
				if (init_electrons[iele]->dz() < 0.2){
				  passedDzEl = true;
				  if (fabs(1/init_electrons[iele]->E() - 1/init_electrons[iele]->P()) < 0.15){
				    passedEtoPFractionEl = true;
				    if (init_electrons[iele]->passConversion()) {
				      passedConversionEl = true;
				      if (init_electrons[iele]->missingHits() <= 1) {
					passedMissingHits = true;
					//		  postCut_electrons.push_back(*init_electrons[iele]);
					//	    postCut_electronsTLV.push_back(*init_electrons[iele]);// fill a new vector with the electrons that passed the cuts
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


	    // specific value
	    /*
	    // electrons
	    for(int iele=0; iele<init_electrons.size(); iele++){
              if (init_electrons[iele]->Pt() > 25){
		passedPtEl = true;
                if (abs(init_electrons[iele]->Eta()) < 2.5){
		  passedEtaEl =true;
		  if (init_electrons[iele]->relPfIso(3, 0.5) < 0.074355){
		    passedIsoEl = true;
		    if (init_electrons[iele]->Pt() > 30){
		      passedNewPtEl = true;
		      if (fabs(init_electrons[iele]->superClusterEta()) <= 1.479){
			passedNewEtaEl = true;
			if (fabs(init_electrons[iele]->deltaEtaIn()) < 0.006574){
			  passedDeltaEtaInEl = true;
			  if (fabs(init_electrons[iele]->deltaPhiIn()) < 0.022868){
			    passedDeltaPhiInEl = true;
			    if (init_electrons[iele]->hadronicOverEm() < 0.037553){
			      passedHadronicOverEmEl = true;
			      if (init_electrons[iele]->d0() < 0.009924){
				passedD0El = true;
				if (init_electrons[iele]->dz() < 0.015310){
				  passedDzEl = true;
				  if (fabs(1/init_electrons[iele]->E() - 1/init_electrons[iele]->P()) < 0.131191){
				    passedEtoPFractionEl = true;
				    if (init_electrons[iele]->passConversion()) {
				      passedConversionEl = true;
				      if (init_electrons[iele]->missingHits() <= 1) {
					passedMissingHits = true;
					//		  postCut_electrons.push_back(*init_electrons[iele]);
					//	    postCut_electronsTLV.push_back(*init_electrons[iele]);// fill a new vector with the electrons that passed the cuts
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
	    */




	    // Synch cut on muon: pt, eta, Iso, ID
	    // Declare one bool per cut
	    Bool_t passedPtMu = false;
	    Bool_t passedEtaMu = false;
	    Bool_t passedIsoMu = false;
	    Bool_t passedGlobalIdMu = false;
	    Bool_t passedNewPtMu = false;
	    Bool_t passedNewEtaMu = false;
	    Bool_t passedChi2Mu = false;
	    Bool_t passedNofTrackerLayersWithMeasurement = false;
	    Bool_t passedNofValidMuHits = false;
	    Bool_t passedD0Mu = false;
	    Bool_t passedDzMu = false;
	    Bool_t passedNofValidPixelHitsMu = false;
	    Bool_t passedNofMatchedStationsMu = false;


	    // muons
	    for(int imuo=0; imuo<init_muons.size(); imuo++){
	      if (abs(init_muons[imuo]->Pt()) > 25){
		passedPtMu = true;
		if (abs(init_muons[imuo]->Eta()) < 2.5){
		  passedEtaMu = true;
		  float muRelIso = (init_muons[imuo]->chargedHadronIso(4) + max( 0.0, init_muons[imuo]->neutralHadronIso(4) + init_muons[imuo]->photonIso(4) - 0.5*init_muons[imuo]->puChargedHadronIso(4) ) )  / init_muons[imuo]->Pt();
		  //		  if (1){ // faco
		  if (muRelIso  < 0.12 ){ 
		    passedIsoMu = true;
		    if (init_muons[imuo]->isGlobalMuon()){
		      passedGlobalIdMu = true;
		      if( init_muons[imuo]->Pt()>26){
			passedNewPtMu = true;
			if (fabs(init_muons[imuo]->Eta())<2.1){
			  passedNewEtaMu = true;
			  if ( init_muons[imuo]->chi2() < 10){
			    passedChi2Mu = true; 
			    //			    if (1){
			    if (init_muons[imuo]->nofTrackerLayersWithMeasurement() > 5){
			      passedNofTrackerLayersWithMeasurement =true;
			      //			      if (1){
			      if (init_muons[imuo]->nofValidMuHits() > 0) {
				passedNofValidMuHits =true;
				if (fabs(init_muons[imuo]->d0()) < 0.2){
				  passedD0Mu = true;
				  if (fabs(init_muons[imuo]->dz()) < 0.5){
				    passedDzMu =true;
				    if(init_muons[imuo]->nofValidPixelHits() > 0){
				      passedNofValidPixelHitsMu =true;
				      if(init_muons[imuo]->nofMatchedStations() > 1){
					passedNofMatchedStationsMu =true;
					//					postCut_muons.push_back(*init_muons[imuo]);
					//					postCut_muonsTLV.push_back(*init_muons[imuo]);//fill a new vector with the muons that passed the cuts    
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






	    /*

	    for(int iele=0; iele<init_electronsTLV.size(); iele++){
	      if (abs(init_electronsTLV[iele].Pt()) > 25){
		passedPtEl = true;
		if (abs(init_electronsTLV[iele].Eta()) < 2.5){
		  passedEtaEl = true;
		  for(int imuo=0; imuo<init_muonsTLV.size(); imuo++){
		    if (abs(init_muonsTLV[imuo].Pt()) > 25){
		      passedPtMu=true;
		      if (abs(init_muonsTLV[imuo].Eta()) < 2.5){
			passedEtaMu=true;
			if(postCut_electronsTLV.size() == 1 ){                                                                                passedExtraElVeto=true;   
			  if(postCut_muonsTLV.size() == 1){
			    passedExtraMuVeto=true;
			    if(init_electrons[iele]->charge() * init_muons[imuo]->charge() == -1){ // to do! make a new el/muon collection
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


	    //	    cout << "checing the loop... Event number is " << ievt << endl;


	    //making cut flow for electron
	    CutFlowTableElectron.Fill(d,0,scaleFactor*lumiWeight);
	    if(passedPtEl) {
	      CutFlowTableElectron.Fill(d,1,scaleFactor*lumiWeight);
	      if (passedEtaEl){
		CutFlowTableElectron.Fill(d,2,scaleFactor*lumiWeight);
		if (passedIsoEl){
		  CutFlowTableElectron.Fill(d,3,scaleFactor*lumiWeight);
		  if (passedNewPtEl){
		    CutFlowTableElectron.Fill(d,4,scaleFactor*lumiWeight);
		    if (passedNewEtaEl){
		      CutFlowTableElectron.Fill(d,5,scaleFactor*lumiWeight);
		      if (passedDeltaEtaInEl){
			CutFlowTableElectron.Fill(d,6,scaleFactor*lumiWeight);
			if(passedDeltaPhiInEl){
			  CutFlowTableElectron.Fill(d,7,scaleFactor*lumiWeight);
			  if(passedHadronicOverEmEl){
			    CutFlowTableElectron.Fill(d,8,scaleFactor*lumiWeight);
			    if (passedD0El) {
			      CutFlowTableElectron.Fill(d,9,scaleFactor*lumiWeight);
			      if(passedDzEl){
				CutFlowTableElectron.Fill(d,10,scaleFactor*lumiWeight);
				if (passedEtoPFractionEl) {
				  CutFlowTableElectron.Fill(d,11,scaleFactor*lumiWeight);
				  if (passedConversionEl) {
				    CutFlowTableElectron.Fill(d,12,scaleFactor*lumiWeight);
				    if (passedMissingHits) {
				      CutFlowTableElectron.Fill(d,13,scaleFactor*lumiWeight);
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





	    //making cut flow for muon
	    CutFlowTableMuon.Fill(d,0,scaleFactor*lumiWeight);
	    if(passedPtMu) {
	      CutFlowTableMuon.Fill(d,1,scaleFactor*lumiWeight);
	      if (passedEtaMu){
		CutFlowTableMuon.Fill(d,2,scaleFactor*lumiWeight);
		if (passedIsoMu){
		  CutFlowTableMuon.Fill(d,3,scaleFactor*lumiWeight);
		  if (passedGlobalIdMu){
		    CutFlowTableMuon.Fill(d,4,scaleFactor*lumiWeight);
		    if (passedNewPtMu){
		      CutFlowTableMuon.Fill(d,5,scaleFactor*lumiWeight);
		      if ( passedNewEtaMu){
			CutFlowTableMuon.Fill(d,6,scaleFactor*lumiWeight);
			if( passedChi2Mu){
			  CutFlowTableMuon.Fill(d,7,scaleFactor*lumiWeight);
			  if(passedNofTrackerLayersWithMeasurement){
			    CutFlowTableMuon.Fill(d,8,scaleFactor*lumiWeight);
			    if (passedNofValidMuHits) {
			      CutFlowTableMuon.Fill(d,9,scaleFactor*lumiWeight);
			      if( passedD0Mu){
				CutFlowTableMuon.Fill(d,10,scaleFactor*lumiWeight);
				if (passedDzMu) {
				  CutFlowTableMuon.Fill(d,11,scaleFactor*lumiWeight);
				  if (passedNofValidPixelHitsMu) {
				    CutFlowTableMuon.Fill(d,12,scaleFactor*lumiWeight);
				    if (passedNofMatchedStationsMu) {
				      CutFlowTableMuon.Fill(d,13,scaleFactor*lumiWeight);
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


	    //	    Filling the tree
	    /*
	    if (passedExtraElVeto == false)
	      break;
	    if (passedExtraMuVeto == false)
	      break;
	    if(passedElMuOS == false)
	      break;
	    if(passedElMuNotOverlaping == false)
	      break;
	    else{
	    */
	    
		// loop over electrons

	    nElectrons=0;
	    nMuons=0;
	    nJets=0;

	    //	    if (passedExtraElVeto == true && passedExtraMuVeto == true && passedElMuOS == true && passedElMuNotOverlaping ==true )
	    if (1)
	      {
		for(int iele=0; iele<postCut_electronsTLV.size() ; iele++){
		  //	    for(int iele=0; iele<init_electrons.size() ; iele++){
		  pX_electron[nElectrons]=postCut_electronsTLV[iele].Px();
		  pY_electron[nElectrons]=postCut_electronsTLV[iele].Py();
		  pZ_electron[nElectrons]=postCut_electronsTLV[iele].Pz();
		  E_electron[nElectrons]=postCut_electronsTLV[iele].E();
		  //		d0_electron[nElectrons]=init_electrons[iele]->d0();
		  //		chargedHadronIso_electron[nElectrons]=init_electrons[iele]->chargedHadronIso(3);
		  //		neutralHadronIso_electron[nElectrons]=init_electrons[iele]->neutralHadronIso(3);
		  //		photonIso_electron[nElectrons]=init_electrons[iele]->photonIso(3);
		  //		pfIso_electron[nElectrons]=init_electrons[iele]->relPfIso(3,0);
		  //		charge_electron[nElectrons]=init_electrons[iele]->charge();
		  nElectrons++;
		}
		//		/*
		// loop over muons
		for(int imuo=0; imuo<postCut_muonsTLV.size(); imuo++){
		  pX_muon[nMuons]=postCut_muonsTLV[imuo].Px();
		  pY_muon[nMuons]=postCut_muonsTLV[imuo].Py();
		  pZ_muon[nMuons]=postCut_muonsTLV[imuo].Pz();
		  E_muon[nMuons]=postCut_muonsTLV[imuo].E();
		  /*
		  d0_muon[nMuons]=init_muons[imuo]->d0();
		  chargedHadronIso_muon[nMuons]=init_muons[imuo]->chargedHadronIso(4);
		  neutralHadronIso_muon[nMuons]=init_muons[imuo]->neutralHadronIso(4);
		  photonIso_muon[nMuons]=init_muons[imuo]->photonIso(4);
		  pfIso_muon[nMuons]=init_muons[imuo]->relPfIso(4,0);
		  charge_muon[nMuons]=init_muons[imuo]->charge();
		  */
		  nMuons++;
		}
	      }
		// loop over jets
		for(int ijet=0; ijet<init_jets.size() && nJets<10; ijet++){
		  pX_jet[nJets]=init_jets[ijet]->Px();
		  pY_jet[nJets]=init_jets[ijet]->Py();
		  pZ_jet[nJets]=init_jets[ijet]->Pz();
		  E_jet[nJets]=init_jets[ijet]->E();
		  dx_jet[nJets]=init_jets[ijet]->vx();
		  dy_jet[nJets]=init_jets[ijet]->vy();
		  nJets++;
		}
		//}

	    
	    // test cut: exactly 2 electrons
	    Bool_t passedNelectrons = false;

	    if (nElectrons == 2){
	      passedNelectrons=true;	      
	    }
	    

	    // Print info

	    
	    
	    if (1){
	      cout << "New EVENT!!!!" << endl;
	      cout << "the event number is " << ievt << endl; 

	      string ws = "    ";


	      //beam spot position
	      double bx, by, bz;
	      bx=0.243996;
	      by=0.392837;
	      bz=0.381877;

	      // declare variables usefull to calculate d0 by hand
	      double px, py, pt, vx, vy, vz, myd0BeamSpot;

	      cout << "there is " << init_electrons.size() << " electron(s) in that event!" << endl;
	      
	      //loop on electron
	      for(int iele=0; iele<init_electrons.size(); iele++){
		px=init_electrons[iele]->Px();
		py=init_electrons[iele]->Py();
		pt=init_electrons[iele]->Pt();
		vx=init_electrons[iele]->vx();
		vy=init_electrons[iele]->vy();
		myd0BeamSpot = (-(vx-bx)*py + (vy-by)*px)/pt;

		cout << ws << "(-(vx-bx)*py + (vy-by)*px)/pt = " << myd0BeamSpot << endl;
		cout << ws << "d0BeamSpot() is: " << init_electrons[iele]->d0BeamSpot() << endl;
		cout << ws << "dzBeamSpot() is: " << init_electrons[iele]->dzBeamSpot() << endl;
		cout << ws << "chargedHadronIso(3) is: " << init_electrons[iele]->chargedHadronIso(3) << endl;
		cout << ws << "neutralHadronIso(3) is: " << init_electrons[iele]->neutralHadronIso(3) << endl;
		cout << ws << "photonIso(3) is: " << init_electrons[iele]->photonIso(3) << endl;
		cout << ws << "puChargedHadronIso(3) is: " << init_electrons[iele]->puChargedHadronIso(3) << endl;
		cout << ws << "relPfIso(3, 0.5) is: " << init_electrons[iele]->relPfIso(3, 0.5) << endl;
		cout << ws << "relPfIso(3, 0) is: " << init_electrons[iele]->relPfIso(3, 0) << endl;
		cout << ws << "hadronicOverEm() is: " << init_electrons[iele]->hadronicOverEm() << endl;

		cout << ws << "1/ET -1/pt is: " << fabs(1/init_electrons[iele]->E() - 1/init_electrons[iele]->P()) << endl;
		}



	      cout << "there is " << init_muons.size() << " muon(s) in that event!" << endl;
	      
	      // loop on muons
	      for(int imuo=0; imuo<init_muons.size(); imuo++){
		px=init_muons[imuo]->Px();
		py=init_muons[imuo]->Py();
		pt=init_muons[imuo]->Pt();
		vx=init_muons[imuo]->vx();
		vy=init_muons[imuo]->vy();
		myd0BeamSpot = (-(vx-bx)*py + (vy-by)*px)/pt;

		cout << ws << "(-(vx-bx)*py + (vy-by)*px)/pt = " << myd0BeamSpot << endl;
		cout << ws << "d0BeamSpot() is: " << init_muons[imuo]->d0BeamSpot() << endl;
		cout << ws << "dzBeamSpot() is: " << init_muons[imuo]->dzBeamSpot() << endl;
		cout << ws << "chargedHadronIso(4) is: " << init_muons[imuo]->chargedHadronIso(4) << endl;
		cout << ws << "neutralHadronIso(4) is: " << init_muons[imuo]->neutralHadronIso(4) << endl;
		cout << ws << "photonIso(4) is: " << init_muons[imuo]->photonIso(4) << endl;
		cout << ws << "puChargedHadronIso(4) is: " << init_muons[imuo]->puChargedHadronIso(4) << endl;
		cout << ws << "relPfIso(4, 0) is: " << init_muons[imuo]->relPfIso(4, 0) << endl;
		}
	      
	      cout << "End of EVENT" << endl << endl;

	      // loop over electrons                                                                           
	      /*
	      for(int iele=0; iele<selectedElectrons.size() ; iele++){                         
		cout << "    electron # " << iele << ":" << endl; 
                cout << selectedElectrons[iele]->Pt();  
		cout << selectedElectrons[iele]->d0();
		pY_electron[nElectrons]=selectedElectrons[iele]->Py();                                       
                pZ_electron[nElectrons]=selectedElectrons[iele]->Pz();                                       
                E_electron[nElectrons]=selectedElectrons[iele]->E();                                         
                d0_electron[nElectrons]=selectedElectrons[iele]->d0();                                       
                chargedHadronIso_electron[nElectrons]=selectedElectrons[iele]->chargedHadronIso(3);          
                neutralHadronIso_electron[nElectrons]=selectedElectrons[iele]->neutralHadronIso(3);          
                photonIso_electron[nElectrons]=selectedElectrons[iele]->photonIso(3);                        
                pfIso_electron[nElectrons]=selectedElectrons[iele]->relPfIso(3,0);                           
                charge_electron[nElectrons]=selectedElectrons[iele]->charge();                               
                nElectrons++;                                                                                
	      }                                                                                                
	      // loop over muons                                                                               
	      nMuons=0;                                                                                        
	      for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){                                 
                pX_muon[nMuons]=selectedMuons[imuo]->Px();                                                   
                pY_muon[nMuons]=selectedMuons[imuo]->Py();                                                   
                pZ_muon[nMuons]=selectedMuons[imuo]->Pz();                                                   
                E_muon[nMuons]=selectedMuons[imuo]->E();                                                     
                d0_muon[nMuons]=selectedMuons[imuo]->d0();                                                   
                chargedHadronIso_muon[nMuons]=selectedMuons[imuo]->chargedHadronIso(4);                      
                neutralHadronIso_muon[nMuons]=selectedMuons[imuo]->neutralHadronIso(4);                      
                photonIso_muon[nMuons]=selectedMuons[imuo]->photonIso(4);                                    
                pfIso_muon[nMuons]=selectedMuons[imuo]->relPfIso(4,0);                                       
                charge_muon[nMuons]=selectedMuons[imuo]->charge();                                           
                nMuons++;                                                                                    
	      } 
	      */


	    }
		  
	    
	    
	    myTree->Fill();

	}


	    
	
		


        
        cout << "number of (weighted) selected jets: " << jeteffaverage[1] ; if(jeteffaverage[0]>0) cout << " at efficiency of " << jeteffaverage[1]/jeteffaverage[0];
        cout << endl;
        cout << "number of (weighted) selected electrons: " << eleeffaverage[1] ; if(eleeffaverage[0]>0) cout << " at efficiency of " << eleeffaverage[1]/eleeffaverage[0];
        cout << endl;
        cout << "number of (weighted) selected muons: " << muoneffaverage[1] ; if(muoneffaverage[0]>0) cout << " at efficiency of " << muoneffaverage[1]/muoneffaverage[0];
        cout << endl;


        
	cout<<endl;
        
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        myTree->Write();
        fileout->Write();
        fileout->Close();
	//        delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    

    CutFlowTableElectron.TableCalculator(false, true, true, true, true);
    string selectiontableElectron = "SelectionTableElectron_table.tex";
    CutFlowTableElectron.Write(selectiontableElectron.c_str(), true,true,true,true,true,true,false);


    CutFlowTableMuon.TableCalculator(false, true, true, true, true);
    string selectiontableMuon = "SelectionTableMuon_table.tex";
    CutFlowTableMuon.Write(selectiontableMuon.c_str(), true,true,true,true,true,true,false);


    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
