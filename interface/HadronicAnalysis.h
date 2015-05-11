#ifndef FORMAT_H
#define FORMAT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TChain.h>
#include <TTree.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//

class HadronicAnalysis : public edm::EDAnalyzer{
	public:
		explicit HadronicAnalysis(const edm::ParameterSet&);
		~HadronicAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;

		edm::Service<TFileService> fs;
		TChain*            chain_;
 
		//const double jetPtMin_;
		//const double bJetPtMin_;
		//const double bJetCSVDiscMin_;
		//const double bJetCSVDiscMax_;

		//const edm::ParameterSet jetSelParams_ ; 
		//const edm::ParameterSet evtSelParams_ ; 

		//EvtInfoBranches    EvtInfo;
		//VertexInfoBranches VtxInfo;
		//GenInfoBranches    GenInfo;
		//JetInfoBranches    GenJetInfo;
		//JetInfoBranches    JetInfo;
		//LepInfoBranches    LepInfo;

		//TH1InfoClass<TH1D> h1;

		//bool   McFlagana;
		//int RunNo;
		//long int EvtNo;
		//int LumiNo;

};

#endif
