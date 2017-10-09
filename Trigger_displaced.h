#ifndef INTERFACE_TRIGGER_H
#define INTERFACE_TRIGGER_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include <iostream>
#include <map>

//TopTree classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"

//TopTreeAnalysisBase classes
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"


class Trigger{
	public:
		Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep);
		~Trigger();
		void bookTriggers(bool isData);
		void checkAvail(int currentRunTrig, vector<Dataset*> datasets, unsigned int d, TTreeLoader* treeLoader, TRootEvent* event, bool verbose);
		int checkIfFired();

	  std::vector<std::string> triggerList;
	  std::map<std::string,std::pair<int,bool> > triggermap;

	private:
		bool muon;
		bool electron;
		bool singleLep;
		bool doubleLep;
		bool fullHadr;
		bool trigged;
		bool redotrigmap;
		int currentRunTrig;
		int previousRunTrig;
		string currentFilenameTrig;
		string previousFilenameTrig;
		int iFileTrig;
		int treeNumberTrig;
};

#endif
