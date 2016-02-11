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
    //selectedEvts_=0;
    //selectedEvts_El_=0;
    //selectedEvts_Mu_=0;
    nselectedEvts_=0;
    nselectedEvts_El_=0;
    nselectedEvts_Mu_=0;
    LHAPDF::initPDFSet( 1, pdfName_);  
    numberPDF_=LHAPDF::numberPDF(1);
    std::cout<<"[INFO] PDF used "<<pdfName_<<" with "<<numberPDF_<<" weights."<<std::endl;
}

PDFUncertaintyAna::~PDFUncertaintyAna()
{ 
    delete chain_;
}

// ------------ Other function -------------
void PDFUncertaintyAna::CreateTH1Number( int numberPDF, TH1InfoClass<TH1D> &h1, std::string hName, std::string title, std::string xtitle, std::string ytitle, std::string xunit, std::string yunit, int bin, double min, double max )
{
    h1.addNewTH1( hName+"_PDF"+num2str(numberPDF), title, xtitle, ytitle, xunit, yunit, bin, min, max);
}
void PDFUncertaintyAna::FillTH1Number( int numberPDF, double value, TH1InfoClass<TH1D> &h1, std::string hName, double wrt  )
{
    h1.GetTH1(hName+"_PDF"+num2str(numberPDF))->Fill( value, wrt ); 
}
double PDFUncertaintyAna::IntegralTH1Number( int numberPDF, double value, TH1InfoClass<TH1D> &h1, std::string hName  )
{   
    int bin1 = 0;
    int bin2 = h1.GetTH1(hName+"_PDF"+num2str(numberPDF))->GetNbinsX()+1; 
    double sum = h1.GetTH1(hName+"_PDF"+num2str(numberPDF))->Integral( bin1, bin2 ); 
    return sum;
}

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
    h1.addNewTH1( "Evt_Events",              "",                        "",                  "Events", "", "",  1,   1,   2   );
    h1.addNewTH1( "Evt_Channel",             "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl",   "",                     "Mass",              "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_El","",                     "Mass",              "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_Mu","",                     "Mass",              "Events", "", "", 500,  0,   500 );
    for( int i=1; i<=numberPDF_; i++ )
    {
        CreateTH1Number( i, h1, "EvtChi2_Top_Leptonic_Mbl",    "", "Mass", "Events", "", "", 500,  0,   500 );
        CreateTH1Number( i, h1, "EvtChi2_Top_Leptonic_Mbl_El", "", "Mass", "Events", "", "", 500,  0,   500 );
        CreateTH1Number( i, h1, "EvtChi2_Top_Leptonic_Mbl_Mu", "", "Mass", "Events", "", "", 500,  0,   500 );
    }

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
    bool isSignal=false;
    vector<double> pdf_weight_sum(numberPDF_,0);
    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);

        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        h1.GetTH1("Evt_Events")->Fill(1);

             if(  isMuonEvt_ &&  isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(3);
        else if(  isMuonEvt_ && !isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(2);
        else if( !isMuonEvt_ &&  isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(1);
        else if( !isMuonEvt_ && !isEleEvt_ ) h1.GetTH1("Evt_Channel")->Fill(0);
       
 
        isSignal = (isSignal_==1)? true:false;
        double wrtevt = WrtEvt_;
        bool passChi2Cut = false;
        if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ ) passChi2Cut=true;

        std::vector<double> pdf_weights;
        LHAPDF::usePDFMember(1,0);
        double xpdf1 = LHAPDF::xfx( 1, PDFx1_, PDFscale_, PDFid1_);
        double xpdf2 = LHAPDF::xfx( 1, PDFx2_, PDFscale_, PDFid2_);
        double w0 = xpdf1 * xpdf2;

        for( int i=1; i <=numberPDF_; ++i )
        {
            LHAPDF::usePDFMember( 1, i);
            double xpdf1_new = LHAPDF::xfx( 1, PDFx1_, PDFscale_, PDFid1_);
            double xpdf2_new = LHAPDF::xfx( 1, PDFx2_, PDFscale_, PDFid2_);
            double weight = xpdf1_new * xpdf2_new / w0;
            pdf_weights.push_back(weight);
        }

        if( isEleEvt_ ==1 && isMuonEvt_ ==0 && passChi2Cut ) //Electron channel
        {
            nselectedEvts_++;
            nselectedEvts_El_++;
            h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Fill( TopMlb_, wrtevt ); 
            h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_El")->Fill( TopMlb_, wrtevt ); 
        }
        if( isEleEvt_ ==0 && isMuonEvt_ ==1 && passChi2Cut ) // Muon channel
        {
            nselectedEvts_++;
            nselectedEvts_Mu_++;
            h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Fill( TopMlb_, wrtevt ); 
            h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_Mu")->Fill( TopMlb_, wrtevt ); 
        }

        if( isSignal )
        { 
            // Singal (ttbar)
            for( int i=0; i<numberPDF_; i++ ){ pdf_weight_sum[i]+=pdf_weights[i]; }
        }
        else
        {   
            // Backgrounds
            if( isEleEvt_ ==1 && isMuonEvt_ ==0 && passChi2Cut ) //Electron channel
            {
                for( int i=0; i<numberPDF_; i++ )
                {
                    FillTH1Number( i+1, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl",    wrtevt*pdf_weights[i] );
                    FillTH1Number( i+1, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_El", wrtevt*pdf_weights[i] );
                }
            }
            if( isEleEvt_ ==0 && isMuonEvt_ ==1 && passChi2Cut ) // Muon channel
            {
                for( int i=0; i<numberPDF_; i++ )
                {
                    FillTH1Number( i+1, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl",    wrtevt*pdf_weights[i] );
                    FillTH1Number( i+1, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_Mu", wrtevt*pdf_weights[i] );
                }
            }
        }
    }//// [END] entry loop 

    //// ----- * Loop evetns only for signal
    if( isSignal )
    {
        for(int entry=0; entry<maxEvents_; entry++ )
        {
            chain_->GetEntry(entry);
            if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] Signal "<<entry<<" of "<< maxEvents_<<endl; }
            if( maxChi2Cut_ <= minChi2_ || minChi2Cut_ > minChi2_ ) continue;
            if( isEleEvt_ == 0 && isMuonEvt_ ==0 ) continue;

            double wrtevt = WrtEvt_;

            LHAPDF::usePDFMember(1,0);
            double xpdf1 = LHAPDF::xfx( 1, PDFx1_, PDFscale_, PDFid1_ );
            double xpdf2 = LHAPDF::xfx( 1, PDFx2_, PDFscale_, PDFid2_ );
            double w0 = xpdf1 * xpdf2;

            for( int i=1; i <=numberPDF_; ++i )
            {
                LHAPDF::usePDFMember( 1, i);
                double xpdf1_new = LHAPDF::xfx( 1, PDFx1_, PDFscale_, PDFid1_);
                double xpdf2_new = LHAPDF::xfx( 1, PDFx2_, PDFscale_, PDFid2_);
                double pdfwrt = xpdf1_new * xpdf2_new / w0 / pdf_weight_sum[i-1] * maxEvents_;
                if( isEleEvt_ ==1 && isMuonEvt_ ==0 ) //Electron channel
                {
                    FillTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl",    wrtevt*pdfwrt );
                    FillTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_El", wrtevt*pdfwrt );
                }
                if( isEleEvt_ ==0 && isMuonEvt_ ==1 ) // Muon channel
                {
                    FillTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl",    wrtevt*pdfwrt );
                    FillTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_Mu", wrtevt*pdfwrt );
                }
            }
        }//// [END] entry loop
    } //// [END] Signal
}

// ------------ method called once each job just after ending the event loop  ------------
void PDFUncertaintyAna::endJob()
{
    int bin1 = 0;
    int bin2 = h1.GetTH1("EvtChi2_Top_Leptonic_Mbl")->GetNbinsX()+1; 
    double selectedEvts_    = h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Integral( bin1, bin2 );
    double selectedEvts_El_ = h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_El")->Integral( bin1, bin2 );
    double selectedEvts_Mu_ = h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_Mu")->Integral( bin1, bin2 );
    std::cout<<">> [INFO] End of Job!"<<std::endl;
    std::cout<<">> [INFO] Selected events 0 "<<nselectedEvts_<<", el: "<<nselectedEvts_El_<<", mu: "<<nselectedEvts_Mu_<<std::endl;
    std::cout<<">> [INFO] Selected events 1 "<<selectedEvts_<<", el: "<<selectedEvts_El_<<", mu: "<<selectedEvts_Mu_<<std::endl;
    for( int i=1; i <=numberPDF_; ++i )
    {
        double evts    = IntegralTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl"    );
        double evts_el = IntegralTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_El" );
        double evts_mu = IntegralTH1Number( i, TopMlb_, h1, "EvtChi2_Top_Leptonic_Mbl_Mu" );
        std::cout<<">> [INFO] PDF "<<i<<std::endl;
        std::cout<<">>        Evts "<<evts<<", el: "<<evts_el<<", mu: "<<evts_mu<<std::endl;
        std::cout<<">>        Effs "<<evts/selectedEvts_<<", el: "<<evts_el/selectedEvts_El_<<", mu: "<<evts_mu/selectedEvts_Mu_<<std::endl;
        std::cout<<">>        v(%) "<<(evts-selectedEvts_)/selectedEvts_*100<<", el: "<<(evts_el-selectedEvts_El_)/selectedEvts_El_*100<<", mu: "<<(evts_mu-selectedEvts_Mu_)/selectedEvts_Mu_*100<<std::endl;
    }
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
