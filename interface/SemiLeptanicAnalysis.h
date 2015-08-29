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

        std::string int2str( int i=3 );
        template<class TH1>
            void setCutFlow( TH1* h );	
        template<class TH1>
            void fillAsym( TH1* h, double value );
        template<class TH1>
            void setLeptonSelHist( TH1* h );
        template <class Object>
            void get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 );

        bool isIsoLeptonFromJets( Lepton lepton, vector<Jet> jetCol, double dR=0.5 );
 
        float getChi2( Jet jet1, Jet jet2, Jet bjet, float M_top=172.5, float Wth_top=16.3, float M_W=82.9, float Wth_W=9.5 );

        double Obs2( TVector3 isoLep, TVector3 hardJet, TVector3 b,   TVector3 bbar );
        double Obs7( TVector3 beam,   TVector3 b,       TVector3 bbar               );

        // ----------member data ---------------------------

        //// Configurables 

        int                             maxEvents_; 
        const int                       reportEvery_; 
        const std::string               inputTTree_;
        const std::vector<std::string>  inputFiles_;

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
        const double Owrt_;
        const unsigned int NJets_;
        bool  Debug_;
        bool  isSkim_;
        bool  doSaveTree_;

        long long int EvtNo_;
        int BxNo_;
        int LumiNo_;
        int RunNo_;
        int isMuonEvt_;
        int isEleEvt_;

        edm::Service<TFileService> fs;
        TChain* chain_;
        TTree*  newtree_;

        //const double jetPtMin_;

        EvtInfoBranches    EvtInfo;
        VertexInfoBranches VtxInfo;
        JetInfoBranches    JetInfo;
        LepInfoBranches    LepInfo;

        TH1InfoClass<TH1D> h1;
        TH2InfoClass<TH2D> h2;

};

#endif
