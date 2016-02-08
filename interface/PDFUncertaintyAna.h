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
//
// class declaration
//
namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

class PDFUncertaintyAna : public edm::EDAnalyzer{
    public:
        explicit PDFUncertaintyAna(const edm::ParameterSet&);
        ~PDFUncertaintyAna();

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
        std::string pdfName_;
        bool Debug_;

        double pdf_weight0_sum_ ; 
        std::vector<double> pdf_weights_sum_ ; 
        double passed_pdf_weight0_sum_ ; 
        std::vector<double> passed_pdf_weights_sum_ ; 

        long long int EvtNo_;
        int BxNo_;
        int LumiNo_;
        int RunNo_;
        int isMuonEvt_;
        int isEleEvt_;
        int isSignal_;
        int McFlag_;
        int PDFid1_;
        int PDFid2_;
        float PDFx1_;   
        float PDFx2_;   
        float PDFv1_;   
        float PDFv2_;   
        float PDFscale_;
        double TopMlb_;
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

        TH1InfoClass<TH1D> h1;
        TH2InfoClass<TH2D> h2;

};

#endif
