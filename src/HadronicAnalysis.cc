// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/HadronicAnalysis.h" 

//
// constructors and destructor
//
HadronicAnalysis::HadronicAnalysis(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	reportEvery_(iConfig.getParameter<int>("ReportEvery")),
	inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
	inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles"))
	//jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
	//bJetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
	//bJetCSVDiscMin_(iConfig.getParameter<double>("BJetCSVDiscMin")),
	//bJetCSVDiscMax_(iConfig.getParameter<double>("BJetCSVDiscMax")),
	//jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	//evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
{}

HadronicAnalysis::~HadronicAnalysis()
{ 
	delete chain_;
}

// ------------ Other function -------------
// ------------ method called once each job just before starting event loop  ------------
void HadronicAnalysis::beginJob()
{ 
	//h1.CreateTH1(fs); h1.Sumw2();

	chain_  = new TChain(inputTTree_.c_str());

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		f->Close();
	}

	//EvtInfo.Register(chain_);
	//VtxInfo.Register(chain_);
	//GenInfo.Register(chain_);
	//GenJetInfo.Register(chain_,"GenJetInfo");
	//JetInfo.Register(chain_,"JetInfo");
	//FatJetInfo.Register(chain_,"FatJetInfo");
	//SubJetInfo.Register(chain_,"SubJetInfo");
	//LepInfo.Register(chain_);

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	return;  
}

// ------------ method called for each event  ------------
void HadronicAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0 ) return;

	//JetSelector jetSelAK5(jetSelParams_) ; 
	//pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

	for(int entry=0; entry<maxEvents_; entry++){
		if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_;

		//// Event variables 
		//int nGoodVtxs(0);
		//bool isdata(0);
		//double evtwt(1); 
		//double puwt(1); 
  	} //// entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void HadronicAnalysis::endJob(){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HadronicAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HadronicAnalysis);
