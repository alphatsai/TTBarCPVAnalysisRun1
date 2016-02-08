// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <map>
#include <string>
#include <typeinfo>
#include <algorithm>

// Root headers 
#include <TVector3.h>
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

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/functions.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/GenParticle.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/PDFUncertaintyAna.h" 

//
// constructors and destructor
//
PDFUncertaintyAna::PDFUncertaintyAna(const edm::ParameterSet& iConfig) : 
    maxEvents_(   iConfig.getParameter<int>("MaxEvents")), 
    reportEvery_( iConfig.getParameter<int>("ReportEvery")),
    inputTTree_(  iConfig.getParameter<std::string>("InputTTree")),
    inputFiles_(  iConfig.getParameter<std::vector<std::string> >("InputFiles")),
    maxChi2Cut_(  iConfig.getParameter<double>("MaxChi2Cut")),
    minChi2Cut_(  iConfig.getParameter<double>("MinChi2Cut")),
    pdfName_(     iConfig.getParameter<std::string>("PdfName")), 
    Debug_(       iConfig.getParameter<bool>("Debug"))
{
    LHAPDF::initPDFSet( 1, pdfName_);  
    std::cout<<"[INFO] PDF used "<<pdfName_<<" with "<<LHAPDF::numberPDF(1)<<" weights."<<std::endl;
}

PDFUncertaintyAna::~PDFUncertaintyAna()
{ 
    delete chain_;
}

// ------------ Other function -------------
// ------------ method called once each job just before starting event loop  ------------
void PDFUncertaintyAna::beginJob()
{
    // Create TH1D
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.ClearTH1Info();
    h1.addNewTH1( "EvtChi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl",    "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_El", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_Mu", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );

    h1.addNewTH1( "Evt_Events",              "",                        "",                  "Events", "", "",  1,   1,   2   );
    h1.addNewTH1( "Evt_Channel",             "",                        "",                  "Events", "", "",  4,   0,   4   );

    h1.addNewTH1( "Gen_Top_Mass",            "",                    "Mass",                  "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Gen_W_Mass",              "",                    "Mass",                  "Events", "", "", 500,  0,   500 );

    h1.CreateTH1( fs );
    h1.Sumw2();

    // Create TH2D
    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.ClearTH2Info();
    h2.CreateTH2( fs );
    h2.Sumw2();

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
    chain_  = new TChain(inputTTree_.c_str());

    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
    }

    chain_->SetBranchAddress("EvtInfo.McFlag",    &McFlag_    );
    chain_->SetBranchAddress("EvtInfo.PDFid1",    &PDFid1_    );
    chain_->SetBranchAddress("EvtInfo.PDFid2",    &PDFid2_    );
    chain_->SetBranchAddress("EvtInfo.PDFx1",     &PDFx1_     );
    chain_->SetBranchAddress("EvtInfo.PDFx2",     &PDFx2_     );
    chain_->SetBranchAddress("EvtInfo.PDFv1",     &PDFv1_     );
    chain_->SetBranchAddress("EvtInfo.PDFv2",     &PDFv2_     );
    chain_->SetBranchAddress("EvtInfo.PDFscale",  &PDFscale_  );
    chain_->SetBranchAddress("EvtInfo.O2",        &O2_        );
    chain_->SetBranchAddress("EvtInfo.O3",        &O3_        );
    chain_->SetBranchAddress("EvtInfo.O4",        &O4_        );
    chain_->SetBranchAddress("EvtInfo.O7",        &O7_        );
    chain_->SetBranchAddress("EvtInfo.TopMlb",    &TopMlb_    );
    chain_->SetBranchAddress("EvtInfo.MinChi2",   &minChi2_   );
    chain_->SetBranchAddress("EvtInfo.WrtObs",    &WrtObs_    );
    chain_->SetBranchAddress("EvtInfo.WrtEvt",    &WrtEvt_    );
    chain_->SetBranchAddress("EvtInfo.isMuonEvt", &isMuonEvt_ );
    chain_->SetBranchAddress("EvtInfo.isEleEvt",  &isEleEvt_  );
    chain_->SetBranchAddress("EvtInfo.isSignal",  &isSignal_  );

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    //LHAPDF::initPDFSet( 1, pdfName_) ;  
    //unsigned int nweights = 1 ;
    //if( LHAPDF::numberPDF(1) > 1 ){ 
    //    nweights += LHAPDF::numberPDF(1) ;
    //    std::vector<double> weights ; 
    //    weights.reserve(nweights) ; 
    //    for (unsigned int ii = 0; ii < nweights; ++ii) {
    //        weights.push_back(0) ; 
    //    }
    //    pdf_weights_sum_.push_back(weights) ; 
    //    passed_pdf_weights_sum_.push_back(weights) ; 
    //}

    return;  
}

// ------------ method called for each event  ------------
void PDFUncertaintyAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
    using namespace edm;
    using namespace std;

    if(  chain_ == 0 ) return;

    cout<<">> [INFO] Starting analysis loop with "<<maxEvents_<<" events..."<<endl;

    //// ----- * Loop evetns
    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);

        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        h1.GetTH1("Evt_Events")->Fill(1);

             if(  isMuonEvt_ &&  isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(3);
        else if(  isMuonEvt_ && !isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(2);
        else if( !isMuonEvt_ &&  isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(1);
        else if( !isMuonEvt_ && !isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(0);

        double wrtevt = WrtEvt_;

        // ---- * ACP INFO
        if( isMuonEvt_ && !isEleEvt_ )
        {
            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                h1.GetTH1("EvtChi2_O2")->Fill( O2_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O3")->Fill( O3_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O4")->Fill( O4_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O7")->Fill( O7_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O2_Mu")->Fill( O2_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O3_Mu")->Fill( O3_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O4_Mu")->Fill( O4_/WrtObs_, wrtevt );
                h1.GetTH1("EvtChi2_O7_Mu")->Fill( O7_/WrtObs_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym_Mu"), O2_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym_Mu"), O3_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym_Mu"), O4_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym_Mu"), O7_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym"),    O2_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym"),    O3_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym"),    O4_, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym"),    O7_, wrtevt );
            }
        }
        else if( !isMuonEvt_ &&  isEleEvt_ )
        {
            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                h1.GetTH1("EvtChi2_O2")->Fill( O2_/WrtObs_,      wrtevt );
                h1.GetTH1("EvtChi2_O3")->Fill( O3_/WrtObs_,      wrtevt );
                h1.GetTH1("EvtChi2_O4")->Fill( O4_/WrtObs_,      wrtevt );
                h1.GetTH1("EvtChi2_O7")->Fill( O7_/WrtObs_,      wrtevt );
                h1.GetTH1("EvtChi2_O2_El")->Fill( O2_/WrtObs_,   wrtevt );
                h1.GetTH1("EvtChi2_O3_El")->Fill( O3_/WrtObs_,   wrtevt );
                h1.GetTH1("EvtChi2_O4_El")->Fill( O4_/WrtObs_,   wrtevt );
                h1.GetTH1("EvtChi2_O7_El")->Fill( O7_/WrtObs_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym_El"), O2_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym_El"), O3_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym_El"), O4_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym_El"), O7_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym"),    O2_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym"),    O3_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym"),    O4_,   wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym"),    O7_,   wrtevt );
            }
        }

    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void PDFUncertaintyAna::endJob(){
    std::cout<<">> [INFO] End of Job!"<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PDFUncertaintyAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PDFUncertaintyAna);
