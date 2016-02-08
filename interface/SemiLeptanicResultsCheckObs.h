#ifndef SEMILEPTANICRESULTSCHECK_H 
#define SEMILEPTANICRESULTSCHECK_H 

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

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1Info.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicTreeBranches.h"
//
// class declaration
//

class SemiLeptanicResultsCheckObs : public edm::EDAnalyzer{
    public:
        explicit SemiLeptanicResultsCheckObs(const edm::ParameterSet&);
        ~SemiLeptanicResultsCheckObs();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        void checkObsChange( TH1D* h, double recoO, double genO, double wrt=1 );

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
        const double maxChi2Cut_;
        const double minChi2Cut_;
        const double Owrt_;
        const double GenACP_;
        bool doWrtEvt_;
        bool Debug_;

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
        double minChi2_;
        double WrtObs_;
        double WrtEvt_;

        edm::Service<TFileService> fs;
        TChain* chain_;
        TTree*  newtree_;

        EvtInfoBranches    EvtInfo;
        GenInfoBranches    GenInfo;
        //VertexInfoBranches VtxInfo;
        JetInfoBranches    JetInfo;
        LepInfoBranches    LepInfo;
        newBranchesLep     IsoLepInfo; 

        TH1InfoClass<TH1D> h1;
        TH2InfoClass<TH2D> h2;

};

#endif
