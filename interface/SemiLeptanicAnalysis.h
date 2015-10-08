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
        const double maxChi2_;
        const double Owrt_;
        const unsigned int NJets_;
        const int Shift_TopPtReWeight_;
        bool  Debug_;
        bool  isSkim_;
        bool  doSaveTree_;

        long long int EvtNo_;
        int BxNo_;
        int LumiNo_;
        int RunNo_;
        int isMuonEvt_;
        int isEleEvt_;
        double O2_;
        double O3_;
        double O4_;
        double O7_;
        double WrtObs_;
        double WrtEvt_;

        edm::Service<TFileService> fs;
        TChain* chain_;
        TTree*  newtree_;

        //const double jetPtMin_;

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
