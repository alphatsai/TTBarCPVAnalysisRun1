#ifndef SEMILEPTANICANALYSIS_H 
#define SEMILEPTANICANALYSIS_H 

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TVector3.h>
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

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1Info.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicTreeBranches.h"
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
            void setCutFlow( TH1* h );
        template<class TH1>
            void setLeptonSelHist( TH1* h );

        // ----------member data ---------------------------

        //// Configurables 

        int                             maxEvents_; 
        const int                       reportEvery_; 
        const std::vector<std::string>  inputFiles_;
        const std::vector<std::string>  inputJsons_;
        const std::string               inputTTree_;
        const std::string               file_JESUncs_;
        const std::string               file_PUDistMC_;
        const std::string               file_PUDistData_;
        const std::string               hist_PUDistMC_;
        const std::string               hist_PUDistData_;

        vector<int> HLT_MuChannel_;
        vector<int> HLT_ElChannel_;
        edm::ParameterSet selPars_Vertex_;
        edm::ParameterSet selPars_Jet_;
        edm::ParameterSet selPars_BJet_;
        edm::ParameterSet selPars_NonBJet_;
        edm::ParameterSet selPars_LooseLepton_;
        edm::ParameterSet selPars_TightMuon_;
        edm::ParameterSet selPars_TightElectron_;
        const double dR_IsoLeptonFromJets_;
        const double dR_rmElelectronOverlapeMuon_;
        const double maxChi2_;
        const double minChi2_;
        const double Owrt_;
        const unsigned int NJets_;
        const int Shift_JER_;
        const int Shift_JES_;
        const int Shift_BTagSF_;
        const int Shift_TopPtReWeight_;
        const int Shift_TightMuonIDSF_;
        const int Shift_TightMuonIsoSF_;
        const int Shift_TightElectronIDSF_;
        bool  Debug_;
        bool  isSkim_;
        bool  doSaveTree_;

        long long int b_EvtNo_;
        int b_BxNo_;
        int b_LumiNo_;
        int b_RunNo_;
        int b_isMuonEvt_;
        int b_isEleEvt_;
        double b_O2_;
        double b_O3_;
        double b_O4_;
        double b_O7_;
        double b_minChi2_;
        double b_WrtObs_;
        double b_WrtEvt_;

        edm::Service<TFileService> fs;
        TChain* chain_;
        TTree*  newtree_;
        SemiLeptanicTreeBranches newAnaBranches_;

        EvtInfoBranches    EvtInfo;
        GenInfoBranches    GenInfo;
        VertexInfoBranches VtxInfo;
        JetInfoBranches    JetInfo;
        LepInfoBranches    LepInfo;

        TH1InfoClass<TH1D> h1;
        TH2InfoClass<TH2D> h2;

        edm::LumiReWeighting LumiWeights_; 
};

#endif
