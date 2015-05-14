#ifndef HADRONICANALYSIS_H 
#define HADRONICANALYSIS_H 

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1Info.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h"
//
// class declaration
//

class SemiLeptanicAnalysis : public edm::EDAnalyzer{
	public:
		explicit SemiLeptanicAnalysis(const edm::ParameterSet&);
		~SemiLeptanicAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		
		template<class TH1>
		void setCutFlow( TH1* h, std::string channel="lj" );	
		template<class TH1>
		void setObservableHist(TH1* h, string ob="O" );
		template<class TH1>
		void setLeptonSelHist( TH1* h );
		template <class Object>
		void get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 );

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;

		double Owrt_;
		bool Debug_;

		edm::Service<TFileService> fs;
		TChain*            chain_;
 
		//const double jetPtMin_;

		EvtInfoBranches    EvtInfo;
		VertexInfoBranches VtxInfo;
		JetInfoBranches    JetInfo;
		LepInfoBranches    LepInfo;

		TH1InfoClass<TH1D> h1;

};

#endif
