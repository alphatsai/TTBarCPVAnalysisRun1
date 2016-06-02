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

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/functions.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/checkEvtTool.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TopCandidate.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicAnalysis.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicTreeBranches.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorElectron.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorMuon.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorJet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorVertex.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/BTagSFUtil.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/LeptonSFUtil.h" 

//
// constructors and destructor
//
SemiLeptanicAnalysis::SemiLeptanicAnalysis(const edm::ParameterSet& iConfig) : 
    maxEvents_(                   iConfig.getParameter<int>("MaxEvents")), 
    reportEvery_(                 iConfig.getParameter<int>("ReportEvery")),
    inputFiles_(                  iConfig.getParameter<std::vector<std::string>>("InputFiles")),
    inputJsons_(                  iConfig.getParameter<std::vector<std::string>>("InputJsons")),
    inputTTree_(                  iConfig.getParameter<std::string>("InputTTree")),
    file_JESUncs_(                iConfig.getParameter<std::string>("File_JESUncs")),
    file_PUDistMC_(               iConfig.getParameter<std::string>("File_PUDistMC")),
    file_PUDistData_(             iConfig.getParameter<std::string>("File_PUDistData")),
    hist_PUDistMC_(               iConfig.getParameter<std::string>("Hist_PUDistMC")),
    hist_PUDistData_(             iConfig.getParameter<std::string>("Hist_PUDistData")),
    HLT_MuChannel_(               iConfig.getParameter<std::vector<int>>("HLT_MuChannel")),
    HLT_ElChannel_(               iConfig.getParameter<std::vector<int>>("HLT_ElChannel")),
    selPars_Vertex_(              iConfig.getParameter<edm::ParameterSet>("SelPars_Vertex")),
    selPars_Jet_(                 iConfig.getParameter<edm::ParameterSet>("SelPars_Jet")),
    selPars_BJet_(                iConfig.getParameter<edm::ParameterSet>("SelPars_BJet")),
    selPars_NonBJet_(             iConfig.getParameter<edm::ParameterSet>("SelPars_NonBJet")),
    selPars_LooseLepton_(         iConfig.getParameter<edm::ParameterSet>("SelPars_LooseLepton")),
    selPars_TightMuon_(           iConfig.getParameter<edm::ParameterSet>("SelPars_TightMuon")),
    selPars_TightElectron_(       iConfig.getParameter<edm::ParameterSet>("SelPars_TightElectron")),
    dR_IsoLeptonFromJets_(        iConfig.getParameter<double>("dR_IsoLeptonFromJets")),
    dR_rmElelectronOverlapeMuon_( iConfig.getParameter<double>("dR_rmElelectronOverlapeMuon")),
    maxChi2_(                     iConfig.getParameter<double>("MaxChi2")),
    minChi2_(                     iConfig.getParameter<double>("MinChi2")),
    maxMlb_(                      iConfig.getParameter<double>("MaxMlb")),
    minMlb_(                      iConfig.getParameter<double>("MinMlb")),
    Owrt_(                        iConfig.getParameter<double>("Owrt")),
    NJets_(                       iConfig.getParameter<int>("NJets")),
    Shift_JER_(                   iConfig.getParameter<int>("Shift_JER")),
    Shift_JES_(                   iConfig.getParameter<int>("Shift_JES")),
    Shift_BTagSF_(                iConfig.getParameter<int>("Shift_BTagSF")),
    Shift_TopPtReWeight_(         iConfig.getParameter<int>("Shift_TopPtReWeight")),
    Shift_TightMuonIDSF_(         iConfig.getParameter<int>("Shift_TightMuonIDSF")),
    Shift_TightMuonIsoSF_(        iConfig.getParameter<int>("Shift_TightMuonIsoSF")),
    Shift_TightElectronIDSF_(     iConfig.getParameter<int>("Shift_TightElectronIDSF")),
    Debug_(                       iConfig.getParameter<bool>("Debug")),
    isSkim_(                      iConfig.getParameter<bool>("IsSkim")),
    doPDFTree_(                   iConfig.getParameter<bool>("DoPDFTree")),
    doSaveTree_(                  iConfig.getParameter<bool>("DoSaveTree")),
    doFullTree_(                  iConfig.getParameter<bool>("DoFullTree"))
{
    if( inputTTree_.compare("Skim/root") == 0 ){ isSkim_=true; }
    LumiWeights_ = edm::LumiReWeighting( file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_ );
}

SemiLeptanicAnalysis::~SemiLeptanicAnalysis()
{ 
    delete chain_;
}

// ------------ Other function -------------
    template<class TH1>
void SemiLeptanicAnalysis::setCutFlow( TH1* h )
{
    if( isSkim_ )
    {
        h->GetXaxis()->SetBinLabel(1,  "All"                                    );
        h->GetXaxis()->SetBinLabel(2,  "PreSelect"                              );
        h->GetXaxis()->SetBinLabel(3,  "#geq1 goodVtx"                          );
        h->GetXaxis()->SetBinLabel(4,  "HLT"                                    );
        h->GetXaxis()->SetBinLabel(5,  "#geq1 Lep"                              );
        h->GetXaxis()->SetBinLabel(6,  "1 isoLep"                               );
        h->GetXaxis()->SetBinLabel(7,  "veto(Loose #mu)"                        );
        h->GetXaxis()->SetBinLabel(8,  "veto(Loose e)"                          );
        h->GetXaxis()->SetBinLabel(9,  ("#geq"+num2str(NJets_)+" Jets").c_str() );
        h->GetXaxis()->SetBinLabel(10, "#geq2 bjets"                            );
        h->GetXaxis()->SetBinLabel(11, "=2 bjets"                               );
        if( minChi2_ == 0 )
            h->GetXaxis()->SetBinLabel(12, ("#chi^{2}<"+num2str(maxChi2_)).c_str() );
        else
            h->GetXaxis()->SetBinLabel(12, ("#chi^{2}>"+num2str(minChi2_)).c_str() );
        h->GetXaxis()->SetBinLabel(13, ("M_{lb}<"+num2str(maxMlb_)).c_str() );
    }else{
        h->GetXaxis()->SetBinLabel(1,  "All"                                    );
        h->GetXaxis()->SetBinLabel(2,  "#geq1 goodVtx"                          );
        h->GetXaxis()->SetBinLabel(3,  "HLT"                                    );
        h->GetXaxis()->SetBinLabel(4,  "#geq1 Lep"                              );
        h->GetXaxis()->SetBinLabel(5,  "1 isoLep"                               );
        h->GetXaxis()->SetBinLabel(6,  "veto(Loose #mu)"                        );
        h->GetXaxis()->SetBinLabel(7,  "veto(Loose e)"                          );
        h->GetXaxis()->SetBinLabel(8,  ("#geq"+num2str(NJets_)+" Jets").c_str() );
        h->GetXaxis()->SetBinLabel(9,  "#geq2 bjets"                            );
        h->GetXaxis()->SetBinLabel(10,  "=2 bjets"                              );
        if( minChi2_ == 0 )
            h->GetXaxis()->SetBinLabel(11, ("#chi^{2}<"+num2str(maxChi2_)).c_str() );
        else
            h->GetXaxis()->SetBinLabel(11, ("#chi^{2}>"+num2str(minChi2_)).c_str() );
        h->GetXaxis()->SetBinLabel(12, ("M_{lb}<"+num2str(maxMlb_)).c_str() );
    }
    return ;
}

    template<class TH1>
void SemiLeptanicAnalysis::setLeptonSelHist( TH1* h )
{
    h->GetXaxis()->SetBinLabel(1,"1:0:0");
    h->GetXaxis()->SetBinLabel(2,"0:1:0");
    h->GetXaxis()->SetBinLabel(3,"0:0:1");
    h->GetXaxis()->SetBinLabel(4,"1:1:0");
    h->GetXaxis()->SetBinLabel(5,"0:1:1");
    h->GetXaxis()->SetBinLabel(6,"1:0:1");
}

// ------------ method called once each job just before starting event loop  ------------
void SemiLeptanicAnalysis::beginJob()
{
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.addNewTH1( "Evt_O7_Mu",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O7_El",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O7Asym_Mu",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O7Asym_El",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O2_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O2_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O2Asym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2Asym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2Asym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_Ob",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_Ob_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_Ob_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_ObAsym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_ObAsym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_ObAsym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_Oa",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_Oa_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_Oa_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_OaAsym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_OaAsym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_OaAsym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3",                  "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O3_Mu",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O3_El",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O3Asym",              "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_Mu",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_El",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4",                  "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O4_Mu",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O4_El",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "Evt_O4Asym",              "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4Asym_Mu",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4Asym_El",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_isoLep_Pt",           "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Et",           "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_E",            "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Eta",          "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_isoLep_Phi",          "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_isoLep_Pt_Mu",        "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Et_Mu",        "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_E_Mu",         "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Eta_Mu",       "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_isoLep_Phi_Mu",       "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_isoLep_Pt_El",        "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Et_El",        "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_E_El",         "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoLep_Eta_El",       "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_isoLep_Phi_El",       "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardJet_Pt",          "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_M",           "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_E",           "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_Eta",         "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardJet_Phi",         "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardJet_BTag",        "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardJet_Pt_Mu",       "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_M_Mu",        "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_E_Mu",        "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_Eta_Mu",      "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardJet_Phi_Mu",      "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardJet_BTag_Mu",     "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardJet_Pt_El",       "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_M_El",        "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_E_El",        "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_Eta_El",      "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardJet_Phi_El",      "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardJet_BTag_El",     "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet1_Pt",      "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_M",       "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_E",       "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_Eta",     "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet1_Phi",     "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet1_BTag",    "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet1_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet1_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet1_BTag_Mu", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet1_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_M_El",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_E_El",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet1_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet1_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet1_BTag_El", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet2_Pt",      "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_M",       "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_E",       "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_Eta",     "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet2_Phi",     "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet2_BTag",    "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet2_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet2_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet2_BTag_Mu", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_HardNonBJet2_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_M_El",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_E_El",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardNonBJet2_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardNonBJet2_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardNonBJet2_BTag_El", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet1_Pt",      "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_M",       "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_E",       "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_Eta",     "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet1_Phi",     "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet1_BTag",    "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet1_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet1_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet1_BTag_Mu", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet1_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_M_El",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_E_El",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet1_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet1_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet1_BTag_El", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet2_Pt",      "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_M",       "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_E",       "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_Eta",     "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet2_Phi",     "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet2_BTag",    "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet2_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet2_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet2_BTag_Mu", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_TopNonBJet2_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_M_El",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_E_El",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_TopNonBJet2_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_TopNonBJet2_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_TopNonBJet2_BTag_El", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bJet_Pt",             "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_M",              "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_E",              "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_Eta",            "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bJet_Phi",            "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bJet_BTag",           "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bJet_Pt_Mu",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_M_Mu",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_E_Mu",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_Eta_Mu",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bJet_Phi_Mu",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bJet_BTag_Mu",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bJet_Pt_El",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_M_El",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_E_El",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet_Eta_El",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bJet_Phi_El",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bJet_BTag_El",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bbarJet_Pt",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_M",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_E",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_Eta",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bbarJet_Phi",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bbarJet_BTag",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bbarJet_Pt_Mu",       "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_M_Mu",        "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_E_Mu",        "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_Eta_Mu",      "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bbarJet_Phi_Mu",      "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bbarJet_BTag_Mu",     "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_bbarJet_Pt_El",       "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_M_El",        "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_E_El",        "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bbarJet_Eta_El",      "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bbarJet_Phi_El",      "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bbarJet_BTag_El",     "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "Evt_Ht",                  "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "Evt_Ht_Mu",               "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "Evt_Ht_El",               "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "Evt_Top_Hadronic_Chi2",    "",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Chi2_El", "",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Chi2_Mu", "",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass",    "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass_El", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass_Mu", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_MassW",   "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_MassW_El","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_MassW_Mu","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt",      "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt_El",   "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt_Mu",   "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta",     "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta_El",  "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta_Mu",  "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi",     "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi_El",  "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi_Mu",  "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt",      "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt_El",   "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt_Mu",   "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mbl",     "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mbl_El",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mbl_Mu",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi",     "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi_El",  "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi_Mu",  "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );

    h1.addNewTH1( "EvtChi2_O7",                  "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7_Mu",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7_El",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O7Asym",              "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O7Asym_Mu",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O7Asym_El",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O2Asym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2Asym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O2Asym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_Ob",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_Ob_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_Ob_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_ObAsym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_ObAsym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_ObAsym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_Oa",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_Oa_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_Oa_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_OaAsym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_OaAsym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_OaAsym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3",                  "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3_Mu",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3_El",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O3Asym",              "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3Asym_Mu",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O3Asym_El",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4",                  "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4_Mu",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4_El",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtChi2_O4Asym",              "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4Asym_Mu",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_O4Asym_El",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtChi2_isoLep_Pt",           "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Et",           "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_E",            "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Eta",          "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_isoLep_Phi",          "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_isoLep_Pt_Mu",        "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Et_Mu",        "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_E_Mu",         "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Eta_Mu",       "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_isoLep_Phi_Mu",       "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_isoLep_Pt_El",        "pT of isoLepon",             "p_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Et_El",        "ET of isoLepon",             "E_{T}(#mu)",        "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_E_El",         "Energy of isoLepon",         "Energy(#mu)",       "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_isoLep_Eta_El",       "Eta of isoLepon",            "#eta(#mu)",         "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_isoLep_Phi_El",       "Phi of isoLepon",            "#phi(#mu)",         "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardJet_Pt",          "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_M",           "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_E",           "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_Eta",         "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardJet_Phi",         "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardJet_BTag",        "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardJet_Pt_Mu",       "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_M_Mu",        "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_E_Mu",        "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_Eta_Mu",      "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardJet_Phi_Mu",      "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardJet_BTag_Mu",     "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardJet_Pt_El",       "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_M_El",        "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_E_El",        "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardJet_Eta_El",      "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardJet_Phi_El",      "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardJet_BTag_El",     "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Pt",      "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_M",       "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_E",       "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Eta",     "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Phi",     "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_BTag",    "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_BTag_Mu", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_M_El",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_E_El",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet1_BTag_El", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Pt",      "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_M",       "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_E",       "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Eta",     "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Phi",     "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_BTag",    "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_BTag_Mu", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_M_El",    "Mass of HardJet",           "Mass(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_E_El",    "Energy of HardJet",         "Energy(HardJet)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",     "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",     "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_HardNonBJet2_BTag_El", "HardJet b-tagged",          "bTag",              "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Pt",      "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_M",       "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_E",       "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Eta",     "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Phi",     "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_BTag",    "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_BTag_Mu", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_M_El",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_E_El",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet1_BTag_El", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Pt",      "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_M",       "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_E",       "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Eta",     "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Phi",     "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_BTag",    "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Pt_Mu",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_M_Mu",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_E_Mu",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Eta_Mu",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Phi_Mu",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_BTag_Mu", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Pt_El",   "pT of HardJet",             "p_{T}(HardJet)",     "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_M_El",    "Mass of HardJet",           "Mass(HardJet)",      "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_E_El",    "Energy of HardJet",         "Energy(HardJet)",    "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Eta_El",  "Eta of HardJet",            "#eta(HardJet)",      "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_Phi_El",  "Phi of HardJet",            "#phi(HardJet)",      "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_TopNonBJet2_BTag_El", "HardJet b-tagged",          "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bJet_Pt",             "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_M",              "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_E",              "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_Eta",            "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bJet_Phi",            "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bJet_BTag",           "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bJet_Pt_Mu",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_M_Mu",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_E_Mu",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_Eta_Mu",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bJet_Phi_Mu",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bJet_BTag_Mu",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bJet_Pt_El",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_M_El",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_E_El",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bJet_Eta_El",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bJet_Phi_El",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bJet_BTag_El",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bbarJet_Pt",          "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_M",           "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_E",           "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_Eta",         "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bbarJet_Phi",         "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bbarJet_BTag",        "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bbarJet_Pt_Mu",       "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_M_Mu",        "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_E_Mu",        "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_Eta_Mu",      "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bbarJet_Phi_Mu",      "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bbarJet_BTag_Mu",     "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_bbarJet_Pt_El",       "pT of b-Jet",               "p_{T}(B-tagged j)",  "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_M_El",        "Mass of b-Jet",             "Mass(B-tagged j)",   "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_E_El",        "Energy of b-Jet",           "Energy(B-tagged j)", "Events", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_bbarJet_Eta_El",      "Eta of b-Jet",              "#eta(B-tagged j)",   "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_bbarJet_Phi_El",      "Phi of b-Jet",              "#phi(B-tagged j)",   "Events", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_bbarJet_BTag_El",     "b-Jet b-tagged",            "bTag",               "Events", "",    "", 100,  0,   1   );
    h1.addNewTH1( "EvtChi2_Ht",                  "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "EvtChi2_Ht_Mu",               "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "EvtChi2_Ht_El",               "Sum(selected jets)",        "H_{T}(selected j)",  "Events", "GeV", "", 1000,  0,  1000 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Chi2",   "",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Chi2_El","",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Chi2_Mu","",                          "#chi^{2}",          "Events", "",    "", 200,  0,   200 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Mass",   "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Mass_El","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Mass_Mu","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_MassW",   "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_MassW_El","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_MassW_Mu","",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Pt",     "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Pt_El",  "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Pt_Mu",  "",                          "Pt",                "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Eta",    "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Eta_El", "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Eta_Mu", "",                          "Eta",               "Events", "",    "", 100, -5,   5   );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Phi",    "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Phi_El", "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "EvtChi2_Top_Hadronic_Phi_Mu", "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mt",     "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mt_El",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mt_Mu",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl",    "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_El", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_Mu", "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Phi",    "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Phi_El", "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Phi_Mu", "",                          "Phi",               "Events", "",    "", 65,  -3.2, 3.2 );

    h1.addNewTH1( "EvtMlb_O7",                  "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O7_Mu",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O7_El",               "O7",                        "O_{7}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O7Asym",              "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O7Asym_Mu",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O7Asym_El",           "A_{O7}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O2",                  "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O2_Mu",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O2_El",               "O2",                        "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O2Asym",              "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O2Asym_Mu",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O2Asym_El",           "A_{O2}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O3",                  "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O3_Mu",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O3_El",               "O3",                        "O_{3}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O3Asym",              "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O3Asym_Mu",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O3Asym_El",           "A_{O3}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O4",                  "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O4_Mu",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O4_El",               "O4",                        "O_{4}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_O4Asym",              "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O4Asym_Mu",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_O4Asym_El",           "A_{O4}",                    "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_Ob",                  "Ob",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_Ob_Mu",               "Ob",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_Ob_El",               "Ob",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_ObAsym",              "A_{Ob}",                  "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_ObAsym_Mu",           "A_{Ob}",                  "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_ObAsym_El",           "A_{Ob}",                  "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_Oa",                  "O",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_Oa_Mu",               "O",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_Oa_El",               "O",                      "O_{2}",              "Events", "",    "",100,  -5,   5   );
    h1.addNewTH1( "EvtMlb_OaAsym",              "A_{O}",                  "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_OaAsym_Mu",           "A_{O}",                  "",                   "Events", "",    "",  2,  -1,   1   );
    h1.addNewTH1( "EvtMlb_OaAsym_El",           "A_{O}",                  "",                   "Events", "",    "",  2,  -1,   1   );

    h1.addNewTH1( "EvtChi2_NVertexNoWrt",     "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtChi2_NVertex",          "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtChi2_NVertexNoWrt_El",  "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtChi2_NVertex_El",       "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtChi2_NVertexNoWrt_Mu",  "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtChi2_NVertex_Mu",       "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertexNoWrt",         "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertex",              "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertexNoWrt_El",      "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertex_El",           "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertexNoWrt_Mu",      "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NVertex_Mu",           "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertexNoWrt",     "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertex",          "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertexNoWrt_El",  "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertex_El",       "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertexNoWrt_Mu",  "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNJet_NVertex_Mu",       "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertexNoWrt",      "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertex",           "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertexNoWrt_El",   "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertex_El",        "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertexNoWrt_Mu",   "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtHLT_NVertex_Mu",        "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertexNoWrt",    "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertex",         "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertexNoWrt_El", "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertex_El",      "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertexNoWrt_Mu", "Num. of Vertex before wrt", "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "EvtNoCut_NVertex_Mu",      "Num. of Vertex",            "N(Vtx)",            "Events", "",    "", 50,   0,   50  );
    h1.addNewTH1( "Evt_NSelJets",             "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NSelJets_Mu",          "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NSelJets_El",          "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NBJets",               "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NBJets_Mu",            "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NBJets_El",            "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NSelJets",         "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NSelJets_Mu",      "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NSelJets_El",      "Num. of selected jets",     "N(selected j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NBJets",           "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NBJets_Mu",        "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "EvtChi2_NBJets_El",        "Num. of b-jets",            "N(B-tagged j)",     "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_CutFlow_Mu",           "",                          "",                  "Evetns", "",    "", 14,   0,   14  );
    h1.addNewTH1( "Evt_CutFlow_El",           "",                          "",                  "Evetns", "",    "", 14,   0,   14  );
    h1.addNewTH1( "Evt_MuCut",                "isoMu:looseMu:looseEl",     "",                  "Evetns", "",    "", 7,    0,   7   );
    h1.addNewTH1( "Evt_ElCut",                "isoEl:looseMu:looseEl",     "",                  "Evetns", "",    "", 7,    0,   7   );
    h1.addNewTH1( "Evt_Wrtevt_TopPt",         "",                          "",                  "Evetns", "",    "", 200,  0,   2   );
    h1.addNewTH1( "Evt_Wrtevt_BTagSF",        "",                          "",                  "Evetns", "",    "", 200,  0,   2   );
    h1.addNewTH1( "Evt_Wrtevt_TightElIDSF",   "",                          "",                  "Evetns", "",    "", 200,  0,   2   );
    h1.addNewTH1( "Evt_Wrtevt_TightMuIDSF",   "",                          "",                  "Evetns", "",    "", 200,  0,   2   );
    h1.addNewTH1( "Evt_Wrtevt_TightMuIsoSF",  "",                          "",                  "Evetns", "",    "", 200,  0,   2   );
    h1.addNewTH1( "Evt_SameChannel",          "",                          "",                  "Evetns", "",    "", 1,    0,   1   );
    h1.addNewTH1( "Evt_Triger",               "",                          "",                  "",       "",    "", 5900, 0,   5900);
    h1.addNewTH1( "Evt_NJsonExEvts",          "Num. of excluded evts",     "",                  "Events", "",    "", 1,    0,    1 );

    h1.addNewTH1( "JetNoJER_Pt",       "pT of Jet",              "p_{T}(j)",           "Yields", "GeV", "", 500,   0,   500 );
    h1.addNewTH1( "JetNoJER_Eta",      "Eta of Jet",             "#eta(j)",            "Yields", "",    "", 100,  -5,   5   );
    h1.addNewTH1( "JetNoJER_BTag",     "Jet b-tagged",           "bTag",               "Yields", "",    "", 100,   0,   1   );

    h1.CreateTH1( fs );
    h1.Sumw2();

    setCutFlow( h1.GetTH1("Evt_CutFlow"   ));
    setCutFlow( h1.GetTH1("Evt_CutFlow_Mu"));
    setCutFlow( h1.GetTH1("Evt_CutFlow_El"));
    setLeptonSelHist(h1.GetTH1("Evt_MuCut"));
    setLeptonSelHist(h1.GetTH1("Evt_ElCut"));

    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.addNewTH2("TH2_Chi2_vs_MET",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_MET_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_MET_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopLeptonicMbl",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopLeptonicMbl_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopLeptonicMbl_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopHadronicMass",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopHadronicMass_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_TopHadronicMass_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_Ht",                 "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_Ht_El",              "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2_Chi2_vs_Ht_Mu",              "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_Ht",      "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_Ht_El",   "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_Ht_Mu",   "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl",      "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl_El",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl_Mu",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopLeptonicMbl_vs_MET",      "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopLeptonicMbl_vs_MET_El",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopLeptonicMbl_vs_MET_Mu",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2_TopHadronicMassW_vs_dRj1j2",      "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.addNewTH2("TH2_TopHadronicMassW_vs_dRj1j2_El",   "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.addNewTH2("TH2_TopHadronicMassW_vs_dRj1j2_Mu",   "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_MET",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_MET_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_MET_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopHadronicMass",    "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopHadronicMass_El", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_TopHadronicMass_Mu", "",  "", "", "", "",  500,   0,   500,  200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_Ht",                 "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_Ht_El",              "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2Chi2_Chi2_vs_Ht_Mu",              "",  "", "", "", "",  1000,  0,   1000, 200,  0,   200 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_Ht",      "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_Ht_El",   "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_Ht_Mu",   "",  "", "", "", "",  1000,  0,   1000, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl",      "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl_El",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl_Mu",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopLeptonicMbl_vs_MET",      "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopLeptonicMbl_vs_MET_El",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopLeptonicMbl_vs_MET_Mu",   "",  "", "", "", "",  500,  0,   500, 500,  0,   500 );
    h2.addNewTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2",      "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.addNewTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2_El",   "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.addNewTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2_Mu",   "",  "", "", "", "",  500,  0,   500, 1000,  0,   10 );
    h2.CreateTH2( fs );
    h2.Sumw2();

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
    chain_  = new TChain(inputTTree_.c_str());

    double allEvents(0);
    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
        if( isSkim_ ) 
        {
            TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
            TH1D* h = (TH1D*)f->Get("Skim/Evt_CutFlow");
            allEvents += h->GetBinContent(1);
            f->Close();
        }
    }

    EvtInfo.Register(chain_);
    GenInfo.Register(chain_);
    VtxInfo.Register(chain_);
    JetInfo.Register(chain_,"PFJetInfo");
    LepInfo.Register(chain_,"PFLepInfo");
    
    // Do PDF tree
    if( doPDFTree_ )
    {
        pdftree_ = fs->make<TTree>("pdftree", "");
        pdftree_ -> Branch("EvtInfo.McFlag"   , &EvtInfo.McFlag  , "EvtInfo.McFlag/I"   ); 
        pdftree_ -> Branch("EvtInfo.PDFid1"   , &EvtInfo.PDFid1  , "EvtInfo.PDFid1/I"   ); 
        pdftree_ -> Branch("EvtInfo.PDFid2"   , &EvtInfo.PDFid1  , "EvtInfo.PDFid2/I"   ); 
        pdftree_ -> Branch("EvtInfo.PDFx1"    , &EvtInfo.PDFx1   , "EvtInfo.PDFx1/F"    ); 
        pdftree_ -> Branch("EvtInfo.PDFx2"    , &EvtInfo.PDFx1   , "EvtInfo.PDFx2/F"    ); 
        pdftree_ -> Branch("EvtInfo.PDFv1"    , &EvtInfo.PDFv1   , "EvtInfo.PDFv1/F"    ); 
        pdftree_ -> Branch("EvtInfo.PDFv2"    , &EvtInfo.PDFv2   , "EvtInfo.PDFv2/F"    ); 
        pdftree_ -> Branch("EvtInfo.PDFscale" , &EvtInfo.PDFscale, "EvtInfo.PDFscale/F" ); 
        pdftree_ -> Branch("EvtInfo.O2"       , &b_O2_           , "EvtInfo.O2/D"       );
        pdftree_ -> Branch("EvtInfo.Ob"       , &b_Ob_           , "EvtInfo.Ob/D"      );
        pdftree_ -> Branch("EvtInfo.Oa"       , &b_Oa_           , "EvtInfo.Oa/D"      );
        pdftree_ -> Branch("EvtInfo.O3"       , &b_O3_           , "EvtInfo.O3/D"       );
        pdftree_ -> Branch("EvtInfo.O4"       , &b_O4_           , "EvtInfo.O4/D"       );
        pdftree_ -> Branch("EvtInfo.O7"       , &b_O7_           , "EvtInfo.O7/D"       );
        pdftree_ -> Branch("EvtInfo.TopMlb"   , &b_TopMlb_       , "EvtInfo.TopMlb/D"   ); 
        pdftree_ -> Branch("EvtInfo.MinChi2"  , &b_minChi2_      , "EvtInfo.MinChi2/D"  );
        pdftree_ -> Branch("EvtInfo.WrtObs"   , &b_WrtObs_       , "EvtInfo.WrtObs/D"   );
        pdftree_ -> Branch("EvtInfo.WrtEvt"   , &b_WrtEvt_       , "EvtInfo.WrtEvt/D"   ); 
        pdftree_ -> Branch("EvtInfo.isMuonEvt", &b_isMuonEvt_    , "EvtInfo.isMuonEvt/I");
        pdftree_ -> Branch("EvtInfo.isEleEvt" , &b_isEleEvt_     , "EvtInfo.isEleEvt/I" );
        pdftree_ -> Branch("EvtInfo.isSignal" , &b_isSignal_     , "EvtInfo.isSignal/I" );
    }

    // Do analysis tree
    if( doFullTree_ )
    {
        fs->cd();
        newtree_ = chain_->CloneTree(0);
        newtree_->Branch("EvtInfo.TopMlb",    &b_TopMlb_,    "EvtInfo.TopMlb/D"    ); 
        newtree_->Branch("EvtInfo.O2",        &b_O2_,        "EvtInfo.O2/D"        );
        newtree_->Branch("EvtInfo.Ob",        &b_Ob_,        "EvtInfo.Ob/D"        );
        newtree_->Branch("EvtInfo.Oa",        &b_Oa_,        "EvtInfo.Oa/D"        );
        newtree_->Branch("EvtInfo.O3",        &b_O3_,        "EvtInfo.O3/D"        );
        newtree_->Branch("EvtInfo.O4",        &b_O4_,        "EvtInfo.O4/D"        );
        newtree_->Branch("EvtInfo.O7",        &b_O7_,        "EvtInfo.O7/D"        );
        newtree_->Branch("EvtInfo.MinChi2",   &b_minChi2_,   "EvtInfo.MinChi2/D"   );
        newtree_->Branch("EvtInfo.WrtObs",    &b_WrtObs_,    "EvtInfo.WrtObs/D"    );
        newtree_->Branch("EvtInfo.WrtEvt",    &b_WrtEvt_,    "EvtInfo.WrtEvt/D"    );
        newtree_->Branch("EvtInfo.isMuonEvt", &b_isMuonEvt_, "EvtInfo.isMuonEvt/I" );
        newtree_->Branch("EvtInfo.isEleEvt",  &b_isEleEvt_,  "EvtInfo.isEleEvt/I"  );
        newAnaBranches_.RegisterTree(newtree_);
    }
    else if( doSaveTree_ || !doFullTree_ )
    {
        newtree_ = fs->make<TTree>("tree", "");
        newtree_->Branch("EvtInfo.RunNo",     &b_RunNo_,     "EvtInfo.RunNo/I"     );
        newtree_->Branch("EvtInfo.EvtNo",     &b_EvtNo_,     "EvtInfo.EvtNo/L"     );
        newtree_->Branch("EvtInfo.BxNo",      &b_BxNo_,      "EvtInfo.BxNo/I"      );
        newtree_->Branch("EvtInfo.LumiNo",    &b_LumiNo_,    "EvtInfo.LumiNo/I"    );
        newtree_->Branch("EvtInfo.TopMlb",    &b_TopMlb_,    "EvtInfo.TopMlb/D"    ); 
        newtree_->Branch("EvtInfo.O2",        &b_O2_,        "EvtInfo.O2/D"        );
        newtree_->Branch("EvtInfo.O3",        &b_O3_,        "EvtInfo.O3/D"        );
        newtree_->Branch("EvtInfo.O4",        &b_O4_,        "EvtInfo.O4/D"        );
        newtree_->Branch("EvtInfo.O7",        &b_O7_,        "EvtInfo.O7/D"        );
        newtree_->Branch("EvtInfo.MinChi2",   &b_minChi2_,   "EvtInfo.MinChi2/D"   );
        newtree_->Branch("EvtInfo.WrtObs",    &b_WrtObs_,    "EvtInfo.WrtObs/D"    );
        newtree_->Branch("EvtInfo.WrtEvt",    &b_WrtEvt_,    "EvtInfo.WrtEvt/D"    );
        newtree_->Branch("EvtInfo.isMuonEvt", &b_isMuonEvt_, "EvtInfo.isMuonEvt/I" );
        newtree_->Branch("EvtInfo.isEleEvt",  &b_isEleEvt_,  "EvtInfo.isEleEvt/I"  );
        newAnaBranches_.RegisterTree(newtree_);
    }

    if( isSkim_ )
    { 
        h1.GetTH1("Evt_CutFlow"   )->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow_Mu")->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow_El")->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow"   )->SetBinError(1, sqrt(allEvents));
        h1.GetTH1("Evt_CutFlow_Mu")->SetBinError(1, sqrt(allEvents));
        h1.GetTH1("Evt_CutFlow_El")->SetBinError(1, sqrt(allEvents));
    }

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void SemiLeptanicAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
    using namespace edm;
    using namespace std;

    if(  chain_ == 0 ) return;

    checkEvtTool checkEvt(Debug_);
    for( unsigned i=0; i<inputJsons_.size(); ++i)
    {
        checkEvt.addJson( inputJsons_.at(i) );
    }
    checkEvt.makeJsonMap();

    cout<<">> [INFO] Starting analysis loop with "<<maxEvents_<<" events..."<<endl;

    TVector3 ax, ay, az;
    ax.SetXYZ(1, 0, 0);
    ay.SetXYZ(0, 1, 0);
    az.SetXYZ(0, 0, 1);

    //
    //// *** Create all selectors *** -----------------------------------------------------------------------------------------
    //

    SelectorVertex  VtxSelection( selPars_Vertex_, Debug_ );

    SelectorJet     JetSelection(     selPars_Jet_, Debug_ );
    SelectorJet    BJetSelection(    selPars_BJet_, Debug_ );
    SelectorJet NonBJetSelection( selPars_NonBJet_, Debug_ );

    SelectorMuon MuonSelectionTight(   selPars_TightMuon_, Debug_ );
    SelectorMuon MuonSelectionLoose( selPars_LooseLepton_, Debug_ );

    SelectorElectron ElectronSelectionLoose(   selPars_LooseLepton_, Debug_ );
    SelectorElectron ElectronSelectionTight( selPars_TightElectron_, Debug_ );

    //
    //// *** Loop evetns *** -------------------------------------------------------------------------------------------------
    //  

    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);
        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        //// *** Check evt by json *** 
        bool isdata(0);
        isdata  = EvtInfo.McFlag ? 0 : 1; 
        if( isdata && !checkEvt.isGoodEvt( EvtInfo.RunNo, EvtInfo.LumiNo ) )
        {
            h1.GetTH1("Evt_NJsonExEvts")->Fill(0); 
            continue;
        }
        h1.GetTH1("Evt_Events")->Fill(1);


        //// *** Check GenInfo ***
        double genTopPt(-1), genAntiTopPt(-1);
        for( int i=0; i<GenInfo.Size; i++)
        {
            if( GenInfo.Status[i] != 3 ) continue;
            if( GenInfo.PdgID[i] ==  6 ) genTopPt=GenInfo.Pt[i]; 
            if( GenInfo.PdgID[i] == -6 ) genAntiTopPt=GenInfo.Pt[i];
        }

        //// *** MC event weight ***
        double wrtevt(1);
        if( !isdata ) wrtevt = GenInfo.Weight;

        //// *** PU reweighting ***
        double wrtevtNoPU(1);
        double wrtevt_pu(1); 
        if( !isdata ) wrtevt_pu *= LumiWeights_.weight( EvtInfo.TrueIT[0] );

        //// *** Top pT reweighting (Only top-pair) ***
        double wrtevt_topPt(1); 
        if( !isdata && genTopPt>=0. && genAntiTopPt>=0. )
        {
            double wrt = sqrt(exp(0.156-0.00137*genTopPt)*exp(0.156-0.00137*genAntiTopPt));
            if(      Shift_TopPtReWeight_ == -1 ) wrtevt_topPt = 1;
            else if( Shift_TopPtReWeight_ ==  1 ) wrtevt_topPt = wrt*wrt;
            else wrtevt_topPt = wrt;
        }
        wrtevt *= wrtevt_pu;
        wrtevt *= wrtevt_topPt;
        wrtevtNoPU *= wrtevt_topPt;
        
        //// *** HLT selection ***
        for( int i=0; i<EvtInfo.nTrgBook; i++){
            if( int(EvtInfo.TrgBook[i]) == 1)
                h1.GetTH1("Evt_Triger")->Fill(i);
        }
        bool passMuonHLT=false;
        bool passElectronHLT=false;
        for( std::vector<int>::const_iterator ihlt = HLT_ElChannel_.begin(); ihlt != HLT_ElChannel_.end(); ++ihlt )
        {
            if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
                passElectronHLT=true;
                break;
            }
        }
        for( std::vector<int>::const_iterator ihlt = HLT_MuChannel_.begin(); ihlt != HLT_MuChannel_.end(); ++ihlt )
        {
            if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
                passMuonHLT=true;
                break;
            }
        }

        //// *** Vertex selection ***
        vector<Vertex> VxtColSelected;
        for( int idx=0; idx < VtxInfo.Size; idx++)
        {
            Vertex vtx( VtxInfo, idx );
            if( VtxSelection.isPass(vtx) ) VxtColSelected.push_back(vtx);
        }
        h1.GetTH1("EvtNoCut_NVertex"        )->Fill(VxtColSelected.size(), wrtevt     );
        h1.GetTH1("EvtNoCut_NVertexNoWrt"   )->Fill(VxtColSelected.size(), wrtevtNoPU );
        h1.GetTH1("EvtNoCut_NVertex_El"     )->Fill(VxtColSelected.size(), wrtevt     );
        h1.GetTH1("EvtNoCut_NVertexNoWrt_El")->Fill(VxtColSelected.size(), wrtevtNoPU );
        h1.GetTH1("EvtNoCut_NVertex_Mu"     )->Fill(VxtColSelected.size(), wrtevt     );
        h1.GetTH1("EvtNoCut_NVertexNoWrt_Mu")->Fill(VxtColSelected.size(), wrtevtNoPU );

        //// *** Jet selection ***
        vector<Jet> JetColSelected, BJetCol, nonBJetCol;
        for( int idx=0; idx < JetInfo.Size; idx++)
        {
            Jet jet( JetInfo, idx );
            h1.GetTH1("JetNoJER_Pt"  )->Fill( jet.Pt                 , wrtevt);
            h1.GetTH1("JetNoJER_Eta" )->Fill( jet.Eta                , wrtevt);
            h1.GetTH1("JetNoJER_BTag")->Fill( jet.CombinedSVBJetTags , wrtevt);

            if( !isdata ) jet.applyJER( Shift_JER_ );
            if( !isdata ) jet.applyJES( Shift_JES_, file_JESUncs_ );

            h1.GetTH1("Jet_Pt"  )->Fill( jet.Pt                 , wrtevt);
            h1.GetTH1("Jet_Eta" )->Fill( jet.Eta                , wrtevt);
            h1.GetTH1("Jet_BTag")->Fill( jet.CombinedSVBJetTags , wrtevt);

            if( JetSelection.isPass(jet) )
            { 
                JetColSelected.push_back(jet);
                h1.GetTH1("SelJet_Pt"  )->Fill( jet.Pt                 , wrtevt);
                h1.GetTH1("SelJet_Eta" )->Fill( jet.Eta                , wrtevt);
                h1.GetTH1("SelJet_BTag")->Fill( jet.CombinedSVBJetTags , wrtevt);
            }
            if( BJetSelection.isPass(jet) )
            { 
                BJetCol.push_back(jet);
                h1.GetTH1("bJet_Pt"  )->Fill( jet.Pt                 , wrtevt);
                h1.GetTH1("bJet_Eta" )->Fill( jet.Eta                , wrtevt);
                h1.GetTH1("bJet_BTag")->Fill( jet.CombinedSVBJetTags , wrtevt);
            }
            if( NonBJetSelection.isPass(jet) ) nonBJetCol.push_back(jet);
        }
        h1.GetTH1("Evt_NJets")->Fill( JetInfo.Size, wrtevt);
        if( JetColSelected.size() != ( nonBJetCol.size()+BJetCol.size()) ) std::cout<<">> [WARNING] JetColSelected.size() "<<JetColSelected.size()<<" != ( nonBJetCol.size() "<<nonBJetCol.size()<<" + BJetCol.size() "<<BJetCol.size()<<" )"<<std::endl;

        //// *** Lepton selection ***
        vector<Lepton> MuColTight, MuColLoose_MuChannel, ElColLoose_MuChannel;
        vector<Lepton> ElColTight, MuColLoose_ElChannel, ElColLoose_ElChannel;
        for( int idx=0; idx < LepInfo.Size; idx++)
        {
            Lepton lepton( LepInfo, idx, EvtInfo.RhoPU[0] );
            //// ** Electron selections 
            if( lepton.LeptonType == 11 )
            {
                if( ElectronSelectionLoose.isPass(lepton) ) ElColLoose_MuChannel.push_back(lepton);
                if( ElectronSelectionLoose.isPass(lepton) &&
                    lepton.Et < ElectronSelectionTight.getCut("lepEtMin")) 
                {
                    ElColLoose_ElChannel.push_back(lepton);
                }
                if( ElectronSelectionTight.isPass(lepton) &&
                    isIsoLeptonFromJets( lepton, BJetCol,    dR_IsoLeptonFromJets_ ) && 
                    isIsoLeptonFromJets( lepton, nonBJetCol, dR_IsoLeptonFromJets_ )) 
                {
                    ElColTight.push_back(lepton);
                }
            }
            //// ** Muon selections
            if( lepton.LeptonType == 13 )
            {
                if( MuonSelectionLoose.isPass(lepton) ) MuColLoose_ElChannel.push_back(lepton);
                if( MuonSelectionLoose.isPass(lepton) && 
                    lepton.Pt < MuonSelectionTight.getCut("lepPtMin"))
                {
                    MuColLoose_MuChannel.push_back(lepton);
                }
                if( MuonSelectionTight.isPass(lepton) && 
                    isIsoLeptonFromJets( lepton, BJetCol,    dR_IsoLeptonFromJets_ ) && 
                    isIsoLeptonFromJets( lepton, nonBJetCol, dR_IsoLeptonFromJets_ )) 
                {
                    MuColTight.push_back(lepton);
                }
            }
        }
        if( ElColTight.size() > 0 && MuColTight.size() > 0 ) rmElelectronOverlapeMuon( ElColTight, MuColTight, dR_rmElelectronOverlapeMuon_ );

        //
        //// *** Event selection ***
        //  
        float minChi2 = +1E10;
        float Ht = 0.;
        TopCandidate top_hadronic, top_leptonic;         
        Jet hardJet, hardNonBJet1, hardNonBJet2, TopNonBJet1, TopNonBJet2, b_jet, bbar_jet;
        Lepton isoLep;
        bool isGoodMuonEvt(false), isGoodElectronEvt(false), passChi2Cut(false), passMlbCut(false);
        
        if( isSkim_ )
        {
            h1.GetTH1("Evt_CutFlow"   )->Fill("PreSelect", 1);
            h1.GetTH1("Evt_CutFlow_El")->Fill("PreSelect", 1);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("PreSelect", 1);
        }else{
            h1.GetTH1("Evt_CutFlow"   )->Fill("All", 1);
            h1.GetTH1("Evt_CutFlow_El")->Fill("All", 1);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("All", 1);
        }

        //// *** Number of vertex selection ***
        if( VxtColSelected.size() )
        {
            h1.GetTH1("Evt_CutFlow"   )->Fill("#geq1 goodVtx", wrtevt);
            h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 goodVtx", wrtevt);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 goodVtx", wrtevt);

            Lepton isoMu, isoEl;
            LeptonSFUtil leptonSFUtil;
            bool passElectronSel(false), passMuonSel(false);
            float wrtevt_tightElID=1.;
            float wrtevt_tightMuID=1.;
            float wrtevt_tightMuIso=1.;
            //// *** Muon channel *
            if( passMuonHLT )
            {
                h1.GetTH1("EvtHLT_NVertex"        )->Fill(VxtColSelected.size(), wrtevt     );
                h1.GetTH1("EvtHLT_NVertexNoWrt"   )->Fill(VxtColSelected.size(), wrtevtNoPU );
                h1.GetTH1("EvtHLT_NVertex_Mu"     )->Fill(VxtColSelected.size(), wrtevt     );
                h1.GetTH1("EvtHLT_NVertexNoWrt_Mu")->Fill(VxtColSelected.size(), wrtevtNoPU );
                h1.GetTH1("Evt_CutFlow_Mu")->Fill("HLT", wrtevt);
                //// *** Tight muon 
                if( MuColTight.size() > 0 )
                { 
                    h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 Lep", wrtevt);
                    //// ** iso muon 
                    if( MuColTight.size() == 1 )
                    {
                        isoMu = MuColTight[0];
                        if( !isdata )
                        {
                            wrtevt_tightMuID  = leptonSFUtil.getSF_TightMuonID(  isoMu, Shift_TightMuonIDSF_  );
                            wrtevt_tightMuIso = leptonSFUtil.getSF_TightMuonIso( isoMu, Shift_TightMuonIsoSF_ );
                        }
                        h1.GetTH1("Evt_CutFlow_Mu")->Fill("1 isoLep", wrtevt*wrtevt_tightMuID*wrtevt_tightMuIso );
                        //// * muon veto 
                        if( MuColLoose_MuChannel.size() == 0 )
                        {
                            h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose #mu)", wrtevt*wrtevt_tightMuID*wrtevt_tightMuIso );
                            //// * electron veto 
                            if( ElColLoose_MuChannel.size() == 0 )
                            {
                                h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose e)", wrtevt*wrtevt_tightMuID*wrtevt_tightMuIso );
                                passMuonSel =true;
                            } ////[END] electron veto
                        } //// [END] muon veto
                    } //// [END] iso muon
                } //// [END] Tight muon
            } //// [END] Muon channel

            //// *** Electron channel *
            if( passElectronHLT )
            {
                h1.GetTH1("EvtHLT_NVertex"        )->Fill(VxtColSelected.size(), wrtevt     );
                h1.GetTH1("EvtHLT_NVertexNoWrt"   )->Fill(VxtColSelected.size(), wrtevtNoPU );
                h1.GetTH1("EvtHLT_NVertex_El"     )->Fill(VxtColSelected.size(), wrtevt     );
                h1.GetTH1("EvtHLT_NVertexNoWrt_El")->Fill(VxtColSelected.size(), wrtevtNoPU );
                h1.GetTH1("Evt_CutFlow_El")->Fill("HLT", wrtevt);
                //// *** Tight electron
                if( ElColTight.size() > 0 )
                { 
                    h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 Lep", wrtevt);
                    //// ** iso electron 
                    if( ElColTight.size() == 1 ) 
                    {
                        isoEl=ElColTight[0];
                        if( !isdata ) wrtevt_tightElID = leptonSFUtil.getSF_TightElectronID( isoEl, Shift_TightElectronIDSF_ );
                        h1.GetTH1("Evt_CutFlow_El")->Fill("1 isoLep", wrtevt*wrtevt_tightElID );
                        //// * muon veto 
                        if( MuColLoose_ElChannel.size() == 0 )
                        {
                            h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose #mu)", wrtevt*wrtevt_tightElID );
                            //// * electron veto 
                            if( ElColLoose_ElChannel.size() == 0 )
                            {
                                h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose e)", wrtevt*wrtevt_tightElID );
                                passElectronSel=true;
                            } //// [END] electron veto
                        } //// [END] muon veto
                    } //// [END] iso electron
                } //// [END] Tight electron
            } //// [END] Electron channel

            //// *** jets and b-jet cut flow *
            if( passElectronSel || passMuonSel ) 
            {
                //// *** Number of jets cuts 
                int unsigned selectedJetsSize = BJetCol.size() + nonBJetCol.size();
                if( selectedJetsSize >= NJets_ ) 
                { 
                    h1.GetTH1("EvtNJet_NVertex"        )->Fill(VxtColSelected.size(), wrtevt     );
                    h1.GetTH1("EvtNJet_NVertexNoWrt"   )->Fill(VxtColSelected.size(), wrtevtNoPU );
                    if( passMuonSel ){     
                        h1.GetTH1("Evt_CutFlow_Mu"         )->Fill( ("#geq"+num2str(NJets_)+" Jets").c_str(), wrtevt*wrtevt_tightMuID*wrtevt_tightMuIso );
                        h1.GetTH1("EvtNJet_NVertex_Mu"     )->Fill( VxtColSelected.size(), wrtevt*wrtevt_tightMuID*wrtevt_tightMuIso     );
                        h1.GetTH1("EvtNJet_NVertexNoWrt_Mu")->Fill( VxtColSelected.size(), wrtevtNoPU*wrtevt_tightMuID*wrtevt_tightMuIso );
                    }
                    if( passElectronSel ){ 
                        h1.GetTH1("Evt_CutFlow_El"         )->Fill(("#geq"+num2str(NJets_)+" Jets").c_str(), wrtevt*wrtevt_tightElID );
                        h1.GetTH1("EvtNJet_NVertex_El"     )->Fill( VxtColSelected.size(), wrtevt*wrtevt_tightElID     );
                        h1.GetTH1("EvtNJet_NVertexNoWrt_El")->Fill( VxtColSelected.size(), wrtevtNoPU*wrtevt_tightElID );
                    }
                    
                    //// ** Number of b-jets cuts, 1
                    if( BJetCol.size() >= 2 )
                    {
                        double wrtevt_btagSF(1);
                        if( !isdata )
                        {
                            BTagSFUtil BTagSF;
                            for( int b=0; b<int(BJetCol.size()); b++ )
                            {
                                wrtevt_btagSF *= BTagSF.getSF( "CSVM", BJetCol[b], Shift_BTagSF_ );
                            }
                        }
                        if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq2 bjets", wrtevt*wrtevt_btagSF*wrtevt_tightMuID*wrtevt_tightMuIso );
                        if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill("#geq2 bjets", wrtevt*wrtevt_btagSF*wrtevt_tightElID                   );
                    }

                    //// ** Number of b-jets cuts, 2
                    if( BJetCol.size() == 2 )
                    {
                        double wrtevt_btagSF(1);
                        if( !isdata )
                        {
                            BTagSFUtil BTagSF;
                            for( int b=0; b<int(BJetCol.size()); b++ )
                            {
                                wrtevt_btagSF *= BTagSF.getSF( "CSVM", BJetCol[b], Shift_BTagSF_ );
                            }
                        }
                        wrtevtNoPU *= wrtevt_btagSF;
                        wrtevt *= wrtevt_btagSF;
                        h1.GetTH1("Evt_Wrtevt_BTagSF")->Fill(wrtevt_btagSF);

                        //// * Fill cutflow hist to each channel
                        if( passMuonSel )
                        {     
                            isGoodMuonEvt = true;
                            isoLep = isoMu;
                            wrtevt *= wrtevt_tightMuID*wrtevt_tightMuIso;
                            wrtevtNoPU *= wrtevt_tightMuID*wrtevt_tightMuIso;
                            h1.GetTH1("Evt_Wrtevt_TightMuIDSF" )->Fill( wrtevt_tightMuID  );
                            h1.GetTH1("Evt_Wrtevt_TightMuIsoSF")->Fill( wrtevt_tightMuIso );
                            h1.GetTH1("Evt_CutFlow_Mu"         )->Fill("=2 bjets", wrtevt);
                        }
                        if( passElectronSel )
                        { 
                            isGoodElectronEvt = true;
                            isoLep = isoEl;
                            wrtevt *= wrtevt_tightElID;
                            wrtevtNoPU *= wrtevt_tightElID;
                            h1.GetTH1("Evt_Wrtevt_TightElIDSF")->Fill( wrtevt_tightElID );
                            h1.GetTH1("Evt_CutFlow_El"        )->Fill("=2 bjets", wrtevt);
                        }

                        //// * Lable the two hardest non_bjet
                        getHighPtObject( JetColSelected, hardJet ); 
                        get2HighPtObject( nonBJetCol, hardNonBJet1, hardNonBJet2 );

                        //// * Distinguish hadronic-top and leptonic-top's b-jet by chi^2
                        //int hadronicTopbjet(-1), leptonicTopbjet(-1); 
                        //for( int bj=0; bj<2; bj++)
                        //{
                        //    float chi2_ = getChi2( hardNonBJet1, hardNonBJet2, BJetCol[bj] );
                        //    if( chi2_ < minChi2 )
                        //    { 
                        //        hadronicTopbjet = bj;
                        //        minChi2 = chi2_;
                        //    }
                        //}
                        //TopNonBJet1 = hardNonBJet1;
                        //TopNonBJet2 = hardNonBJet2;

                        //// * Distinguish hadronic-top and leptonic-top's b-jet by chi^2
                        int topjet1(-1), topjet2(-1), hadronicTopbjet(-1), leptonicTopbjet(-1); 
                        for( int ij1=1; ij1<int(nonBJetCol.size()); ij1++){
                            //for( int ij2=ij1-1; ij2<ij1; ij2++){
                            for( int ij2=0; ij2<ij1; ij2++){
                                for( int bj=0; bj<2; bj++)
                                {
                                    float chi2_ = getChi2( nonBJetCol[ij1], nonBJetCol[ij2], BJetCol[bj] );
                                    if( chi2_ < minChi2 )
                                    { 
                                        topjet1  = ij1;
                                        topjet2  = ij2;
                                        hadronicTopbjet = bj;
                                        minChi2 = chi2_;
                                    }
                                }
                            }
                        }
                        if( nonBJetCol[topjet1].Pt > nonBJetCol[topjet2].Pt )
                        { 
                            TopNonBJet1 = nonBJetCol[topjet1]; TopNonBJet1.Index = topjet1;
                            TopNonBJet2 = nonBJetCol[topjet2]; TopNonBJet2.Index = topjet2;
                        }
                        else
                        {
                            TopNonBJet1 = nonBJetCol[topjet2]; TopNonBJet1.Index = topjet2;
                            TopNonBJet2 = nonBJetCol[topjet1]; TopNonBJet2.Index = topjet1;
                        }

                        //// * reco top  
                        leptonicTopbjet = ( hadronicTopbjet==0 )? 1:0;
                        if( isoLep.Charge < 0 ) // tbar->bbar w- ( w- > l- v )
                        {
                            b_jet    = BJetCol[hadronicTopbjet]; b_jet.Index    = hadronicTopbjet;
                            bbar_jet = BJetCol[leptonicTopbjet]; bbar_jet.Index = leptonicTopbjet;
                            top_hadronic.Fill( b_jet,    TopNonBJet1, TopNonBJet2, 0 );
                            top_leptonic.Fill( bbar_jet, isoLep, EvtInfo.PFMET, EvtInfo.PFMETPhi, 1 );
                        }
                        else if( isoLep.Charge > 0 ) //t->b w+ ( w+ > l+ v )
                        {
                            b_jet    = BJetCol[leptonicTopbjet]; b_jet.Index    = leptonicTopbjet;
                            bbar_jet = BJetCol[hadronicTopbjet]; bbar_jet.Index = hadronicTopbjet;
                            top_hadronic.Fill( bbar_jet, TopNonBJet1, TopNonBJet2, 1 );
                            top_leptonic.Fill( b_jet,    isoLep, EvtInfo.PFMET, EvtInfo.PFMETPhi, 0 );
                        }
                        else
                        { std::cout<<">> [ERROR] There an nuetral lepton!? "<<std::endl; }

                        //// * Fill Ht = scale sum of selected jets
                        for( int j=0; j<int(nonBJetCol.size()); j++ ){ Ht += nonBJetCol[j].Pt; }
                        for( int j=0; j<int(BJetCol.size());    j++ ){ Ht += BJetCol[j].Pt;    }

                        //// * chi^2 cut
                        if( maxChi2_ > minChi2 && minChi2_ <= minChi2 )
                        {
                            passChi2Cut = true;
                            if( minChi2_ == 0 )
                            {
                                if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill(("#chi^{2}<"+num2str(maxChi2_)).c_str(), wrtevt);
                                if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill(("#chi^{2}<"+num2str(maxChi2_)).c_str(), wrtevt);
                            }
                            else
                            {
                                if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill(("#chi^{2}>"+num2str(minChi2_)).c_str(), wrtevt);
                                if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill(("#chi^{2}>"+num2str(minChi2_)).c_str(), wrtevt);
                            }
                            //if( Ht < 250 )
                            //{
                            //    isGoodElectronEvt = false;
                            //    isGoodMuonEvt = false;
                            //}
                        } //// [END] Chi^2 cut

                        ////// * Mlb cut
                        if( passChi2Cut ){
                            if( maxMlb_ > top_leptonic.Massbl && minMlb_ <= top_leptonic.Massbl )
                            {
                                passMlbCut=true;
                                if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill(("M_{lb}<"+num2str(maxMlb_)).c_str(), wrtevt);
                                if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill(("M_{lb}<"+num2str(maxMlb_)).c_str(), wrtevt);
                            } //// [END] Mlb cut
                        } //// [END] Chi^2 cut
                    } //// [END] Number of b-jets cuts, 2
                } //// [END] Number of jets cut
            } //// [END] Jet and bjet cutflow

            //// *** Check events *
            if( passElectronHLT || passMuonHLT ){
                h1.GetTH1("Evt_CutFlow"           )->Fill("HLT", wrtevt);
                if( ( MuColTight.size() + ElColTight.size()) > 0 ) h1.GetTH1("Evt_CutFlow")->Fill("#geq1 Lep", wrtevt);
                if( ( MuColTight.size() + ElColTight.size()) == 1 ){ 
                    h1.GetTH1("Evt_CutFlow")->Fill("1 isoLep", wrtevt);
                    if( MuColTight.size() == 1 && MuColLoose_MuChannel.size() == 0 && ElColLoose_MuChannel.size() == 0 )
                        h1.GetTH1("Evt_MuCut")->Fill("1:0:0", 1);
                    if( MuColTight.size() == 0 && MuColLoose_MuChannel.size() == 1 && ElColLoose_MuChannel.size() == 0 )
                        h1.GetTH1("Evt_MuCut")->Fill("0:1:0", 1);
                    if( MuColTight.size() == 0 && MuColLoose_MuChannel.size() == 0 && ElColLoose_MuChannel.size() == 1 )
                        h1.GetTH1("Evt_MuCut")->Fill("0:0:1", 1);
                    if( MuColTight.size() == 1 && MuColLoose_MuChannel.size() == 1 && ElColLoose_MuChannel.size() == 0 )
                        h1.GetTH1("Evt_MuCut")->Fill("1:1:0", 1);
                    if( MuColTight.size() == 1 && MuColLoose_MuChannel.size() == 0 && ElColLoose_MuChannel.size() == 1 )
                        h1.GetTH1("Evt_MuCut")->Fill("1:0:1", 1);
                    if( MuColTight.size() == 0 && MuColLoose_MuChannel.size() == 1 && ElColLoose_MuChannel.size() == 1 )
                        h1.GetTH1("Evt_MuCut")->Fill("0:1:1", 1);
                    if( MuColTight.size() == 1 && MuColLoose_MuChannel.size() == 1 && ElColLoose_MuChannel.size() == 1 )
                        h1.GetTH1("Evt_MuCut")->Fill("1:1:1", 1);
                    if( ElColTight.size() == 1 && MuColLoose_ElChannel.size() == 0 && ElColLoose_ElChannel.size() == 0 )
                        h1.GetTH1("Evt_ElCut")->Fill("1:0:0", 1);
                    if( ElColTight.size() == 0 && MuColLoose_ElChannel.size() == 1 && ElColLoose_ElChannel.size() == 0 )
                        h1.GetTH1("Evt_ElCut")->Fill("0:1:0", 1);
                    if( ElColTight.size() == 0 && MuColLoose_ElChannel.size() == 0 && ElColLoose_ElChannel.size() == 1 )
                        h1.GetTH1("Evt_ElCut")->Fill("0:0:1", 1);
                    if( ElColTight.size() == 1 && MuColLoose_ElChannel.size() == 1 && ElColLoose_ElChannel.size() == 0 )
                        h1.GetTH1("Evt_ElCut")->Fill("1:1:0", 1);
                    if( ElColTight.size() == 1 && MuColLoose_ElChannel.size() == 0 && ElColLoose_ElChannel.size() == 1 )
                        h1.GetTH1("Evt_ElCut")->Fill("1:0:1", 1);
                    if( ElColTight.size() == 0 && MuColLoose_ElChannel.size() == 1 && ElColLoose_ElChannel.size() == 1 )
                        h1.GetTH1("Evt_ElCut")->Fill("0:1:1", 1);
                    if( ElColTight.size() == 1 && MuColLoose_ElChannel.size() == 1 && ElColLoose_ElChannel.size() == 1 )
                        h1.GetTH1("Evt_ElCut")->Fill("1:1:1", 1);
                }
            }
        } //// [END] Vxt selection
        //// [END] Event selections

        //
        //// *** Fill other events plots ***
        //  
        double O2 = 0;
        double Ob = 0;
        double Oa = 0;
        double O3 = 0;
        double O4 = 0;
        double O7 = 0;
        if( isGoodMuonEvt && isGoodElectronEvt ) h1.GetTH1("Evt_SameChannel")->Fill(0);
        if( isGoodMuonEvt || isGoodElectronEvt )
        {
            Oa = Obs4( isoLep.P3, TopNonBJet1.P3, b_jet.P3, bbar_jet.P3 );
            Ob = Obs2( isoLep.P3, TopNonBJet1.P3, b_jet.P3, bbar_jet.P3, isoLep.Charge );
            O2 = Obs2( isoLep.P3, TopNonBJet1.P3, b_jet.P3, bbar_jet.P3 );
            O3 = Obs3( isoLep.P4, TopNonBJet1.P4, b_jet.P4, bbar_jet.P4, isoLep.Charge );
            O4 = Obs4( isoLep.P3, TopNonBJet1.P3, b_jet.P3, bbar_jet.P3, isoLep.Charge );
            O7 = Obs7( az, b_jet.P3, bbar_jet.P3 );

            if( !doFullTree_ && doSaveTree_ )
            {
                b_RunNo_ = EvtInfo.RunNo;
                b_EvtNo_ = EvtInfo.EvtNo;
                b_BxNo_  = EvtInfo.BxNo;
                b_LumiNo_= EvtInfo.LumiNo;
            }
            if(  doFullTree_ || doSaveTree_ )
            {
                newAnaBranches_.fill_BJetNewBranches( b_jet ); 
                newAnaBranches_.fill_BbarJetNewBranches( bbar_jet ); 
                newAnaBranches_.fill_nonBJetColNewBranches( nonBJetCol );
                newAnaBranches_.fill_topHadronicNewBranches( top_hadronic, TopNonBJet1.Index, TopNonBJet2.Index );
                newAnaBranches_.fill_isoLepNewBranches( isoLep ); 
                b_TopMlb_=top_leptonic.Massbl;
                b_O2_  = O2;
                b_O3_  = O3;
                b_O4_  = O4;
                b_O7_  = O7;
                b_Oa_  = Oa;
                b_Ob_  = Ob;
                b_minChi2_ = minChi2;
                b_WrtObs_ = Owrt_;
                b_WrtEvt_ = wrtevt;
            }

            h2.GetTH2("TH2_Chi2_vs_MET"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
            h2.GetTH2("TH2_Chi2_vs_TopLeptonicMbl"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
            h2.GetTH2("TH2_Chi2_vs_TopHadronicMass"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
            h2.GetTH2("TH2_Chi2_vs_Ht"                       )->Fill( Ht,                  minChi2,             wrtevt );
            h2.GetTH2("TH2_TopHadronicMass_vs_Ht"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
            h2.GetTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
            h2.GetTH2("TH2_TopLeptonicMbl_vs_MET"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
            h2.GetTH2("TH2_TopHadronicMassW_vs_dRj1j2"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
            h1.GetTH1("Evt_Wrtevt_TopPt"      )->Fill( wrtevt_topPt );
            h1.GetTH1("Evt_Ht"                )->Fill( Ht                             , wrtevt     );
            h1.GetTH1("Evt_NVertex"           )->Fill( VxtColSelected.size()          , wrtevt     );
            h1.GetTH1("Evt_NVertexNoWrt"      )->Fill( VxtColSelected.size()          , wrtevtNoPU );
            h1.GetTH1("Evt_NSelJets"          )->Fill( JetColSelected.size()          , wrtevt     );
            h1.GetTH1("Evt_NBJets"            )->Fill( BJetCol.size()                 , wrtevt     );
            h1.GetTH1("Evt_bJet_Pt"           )->Fill( b_jet.Pt                       , wrtevt     );
            h1.GetTH1("Evt_bJet_Eta"          )->Fill( b_jet.Eta                      , wrtevt     );
            h1.GetTH1("Evt_bJet_Phi"          )->Fill( b_jet.Phi                      , wrtevt     );
            h1.GetTH1("Evt_bJet_E"            )->Fill( b_jet.Energy                   , wrtevt     );
            h1.GetTH1("Evt_bJet_M"            )->Fill( b_jet.Mass                     , wrtevt     );
            h1.GetTH1("Evt_bJet_BTag"         )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
            h1.GetTH1("Evt_bbarJet_Pt"        )->Fill( bbar_jet.Pt                    , wrtevt     );
            h1.GetTH1("Evt_bbarJet_Eta"       )->Fill( bbar_jet.Eta                   , wrtevt     );
            h1.GetTH1("Evt_bbarJet_Phi"       )->Fill( bbar_jet.Phi                   , wrtevt     );
            h1.GetTH1("Evt_bbarJet_E"         )->Fill( bbar_jet.Energy                , wrtevt     );
            h1.GetTH1("Evt_bbarJet_M"         )->Fill( bbar_jet.Mass                  , wrtevt     );
            h1.GetTH1("Evt_bbarJet_BTag"      )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
            h1.GetTH1("Evt_HardJet_Pt"        )->Fill( hardJet.Pt                     , wrtevt     );
            h1.GetTH1("Evt_HardJet_Eta"       )->Fill( hardJet.Eta                    , wrtevt     );
            h1.GetTH1("Evt_HardJet_Phi"       )->Fill( hardJet.Phi                    , wrtevt     );
            h1.GetTH1("Evt_HardJet_E"         )->Fill( hardJet.Energy                 , wrtevt     );
            h1.GetTH1("Evt_HardJet_M"         )->Fill( hardJet.Mass                   , wrtevt     );
            h1.GetTH1("Evt_HardJet_BTag"      )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_Pt"   )->Fill( hardNonBJet1.Pt                , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_Eta"  )->Fill( hardNonBJet1.Eta               , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_Phi"  )->Fill( hardNonBJet1.Phi               , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_E"    )->Fill( hardNonBJet1.Energy            , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_M"    )->Fill( hardNonBJet1.Mass              , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet1_BTag" )->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_Pt"   )->Fill( hardNonBJet2.Pt                , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_Eta"  )->Fill( hardNonBJet2.Eta               , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_Phi"  )->Fill( hardNonBJet2.Phi               , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_E"    )->Fill( hardNonBJet2.Energy            , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_M"    )->Fill( hardNonBJet2.Mass              , wrtevt     );
            h1.GetTH1("Evt_HardNonBJet2_BTag" )->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_Pt"    )->Fill( TopNonBJet1.Pt                 , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_Eta"   )->Fill( TopNonBJet1.Eta                , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_Phi"   )->Fill( TopNonBJet1.Phi                , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_E"     )->Fill( TopNonBJet1.Energy             , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_M"     )->Fill( TopNonBJet1.Mass               , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet1_BTag"  )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_Pt"    )->Fill( TopNonBJet2.Pt                 , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_Eta"   )->Fill( TopNonBJet2.Eta                , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_Phi"   )->Fill( TopNonBJet2.Phi                , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_E"     )->Fill( TopNonBJet2.Energy             , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_M"     )->Fill( TopNonBJet2.Mass               , wrtevt     );
            h1.GetTH1("Evt_TopNonBJet2_BTag"  )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
            h1.GetTH1("Evt_isoLep_Pt"         )->Fill( isoLep.Pt                      , wrtevt     );
            h1.GetTH1("Evt_isoLep_Et"         )->Fill( isoLep.Et                      , wrtevt     );
            h1.GetTH1("Evt_isoLep_Eta"        )->Fill( isoLep.Eta                     , wrtevt     );
            h1.GetTH1("Evt_isoLep_Phi"        )->Fill( isoLep.Phi                     , wrtevt     );
            h1.GetTH1("Evt_isoLep_E"          )->Fill( isoLep.Energy                  , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_Chi2" )->Fill( minChi2                        , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_Mass" )->Fill( top_hadronic.Mass              , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_MassW" )->Fill( top_hadronic.MassW            , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_Pt"   )->Fill( top_hadronic.Pt                , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_Eta"  )->Fill( top_hadronic.Eta               , wrtevt     );
            h1.GetTH1("Evt_Top_Hadronic_Phi"  )->Fill( top_hadronic.Phi               , wrtevt     );
            h1.GetTH1("Evt_Top_Leptonic_Mt"   )->Fill( top_leptonic.MassT             , wrtevt     );
            h1.GetTH1("Evt_Top_Leptonic_Mbl"  )->Fill( top_leptonic.Massbl            , wrtevt     );
            h1.GetTH1("Evt_Top_Leptonic_Phi"  )->Fill( top_leptonic.Phi               , wrtevt     );
            h1.GetTH1("Evt_O2"                )->Fill( O2/Owrt_                       , wrtevt     );
            h1.GetTH1("Evt_Oa"                )->Fill( Oa/Owrt_                      , wrtevt     );
            h1.GetTH1("Evt_Ob"                )->Fill( Ob/Owrt_                      , wrtevt     );
            h1.GetTH1("Evt_O3"                )->Fill( O3/Owrt_                       , wrtevt     );
            h1.GetTH1("Evt_O4"                )->Fill( O4/Owrt_                       , wrtevt     );
            h1.GetTH1("Evt_O7"                )->Fill( O7/Owrt_                       , wrtevt     );
            fillAsym( h1.GetTH1("Evt_O2Asym"),  O2,  wrtevt );
            fillAsym( h1.GetTH1("Evt_OaAsym"),  Oa, wrtevt );
            fillAsym( h1.GetTH1("Evt_ObAsym"),  Ob, wrtevt );
            fillAsym( h1.GetTH1("Evt_O3Asym"),  O3,  wrtevt );
            fillAsym( h1.GetTH1("Evt_O4Asym"),  O4,  wrtevt );
            fillAsym( h1.GetTH1("Evt_O7Asym"),  O7,  wrtevt );
            if( passChi2Cut )
            {
                if( passMlbCut )
                {
                    h1.GetTH1("EvtMlb_Ob")->Fill( Ob/Owrt_, wrtevt);
                    h1.GetTH1("EvtMlb_Oa")->Fill( Oa/Owrt_, wrtevt);
                    h1.GetTH1("EvtMlb_O3" )->Fill( O3/Owrt_, wrtevt);
                    h1.GetTH1("EvtMlb_O4" )->Fill( O4/Owrt_, wrtevt);
                    h1.GetTH1("EvtMlb_O7" )->Fill( O7/Owrt_, wrtevt);
                    fillAsym( h1.GetTH1("EvtMlb_ObAsym"), Ob, wrtevt );
                    fillAsym( h1.GetTH1("EvtMlb_OaAsym"), Oa, wrtevt );
                    fillAsym( h1.GetTH1("EvtMlb_O3Asym"), O3, wrtevt );
                    fillAsym( h1.GetTH1("EvtMlb_O4Asym"), O4, wrtevt );
                    fillAsym( h1.GetTH1("EvtMlb_O7Asym"), O7, wrtevt );
                }
                h2.GetTH2("TH2Chi2_Chi2_vs_MET"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
                h2.GetTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
                h2.GetTH2("TH2Chi2_Chi2_vs_TopHadronicMass"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
                h2.GetTH2("TH2Chi2_Chi2_vs_Ht"                       )->Fill( Ht,                  minChi2,             wrtevt );
                h2.GetTH2("TH2Chi2_TopHadronicMass_vs_Ht"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2Chi2_TopLeptonicMbl_vs_MET"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
                h2.GetTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
                h1.GetTH1("EvtChi2_Ht"                )->Fill( Ht                             , wrtevt     );
                h1.GetTH1("EvtChi2_NVertex"           )->Fill( VxtColSelected.size()          , wrtevt     );
                h1.GetTH1("EvtChi2_NVertexNoWrt"      )->Fill( VxtColSelected.size()          , wrtevtNoPU );
                h1.GetTH1("EvtChi2_NSelJets"          )->Fill( JetColSelected.size()          , wrtevt     );
                h1.GetTH1("EvtChi2_NBJets"            )->Fill( BJetCol.size()                 , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_Pt"           )->Fill( b_jet.Pt                       , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_Eta"          )->Fill( b_jet.Eta                      , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_Phi"          )->Fill( b_jet.Phi                      , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_E"            )->Fill( b_jet.Energy                   , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_M"            )->Fill( b_jet.Mass                     , wrtevt     );
                h1.GetTH1("EvtChi2_bJet_BTag"         )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_Pt"        )->Fill( bbar_jet.Pt                    , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_Eta"       )->Fill( bbar_jet.Eta                   , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_Phi"       )->Fill( bbar_jet.Phi                   , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_E"         )->Fill( bbar_jet.Energy                , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_M"         )->Fill( bbar_jet.Mass                  , wrtevt     );
                h1.GetTH1("EvtChi2_bbarJet_BTag"      )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_Pt"        )->Fill( hardJet.Pt                     , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_Eta"       )->Fill( hardJet.Eta                    , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_Phi"       )->Fill( hardJet.Phi                    , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_E"         )->Fill( hardJet.Energy                 , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_M"         )->Fill( hardJet.Mass                   , wrtevt     );
                h1.GetTH1("EvtChi2_HardJet_BTag"      )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_Pt"   )->Fill( hardNonBJet1.Pt                , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_Eta"  )->Fill( hardNonBJet1.Eta               , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_Phi"  )->Fill( hardNonBJet1.Phi               , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_E"    )->Fill( hardNonBJet1.Energy            , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_M"    )->Fill( hardNonBJet1.Mass              , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet1_BTag" )->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_Pt"   )->Fill( hardNonBJet2.Pt                , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_Eta"  )->Fill( hardNonBJet2.Eta               , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_Phi"  )->Fill( hardNonBJet2.Phi               , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_E"    )->Fill( hardNonBJet2.Energy            , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_M"    )->Fill( hardNonBJet2.Mass              , wrtevt     );
                h1.GetTH1("EvtChi2_HardNonBJet2_BTag" )->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_Pt"    )->Fill( TopNonBJet1.Pt                 , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_Eta"   )->Fill( TopNonBJet1.Eta                , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_Phi"   )->Fill( TopNonBJet1.Phi                , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_E"     )->Fill( TopNonBJet1.Energy             , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_M"     )->Fill( TopNonBJet1.Mass               , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet1_BTag"  )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_Pt"    )->Fill( TopNonBJet2.Pt                 , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_Eta"   )->Fill( TopNonBJet2.Eta                , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_Phi"   )->Fill( TopNonBJet2.Phi                , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_E"     )->Fill( TopNonBJet2.Energy             , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_M"     )->Fill( TopNonBJet2.Mass               , wrtevt     );
                h1.GetTH1("EvtChi2_TopNonBJet2_BTag"  )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("EvtChi2_isoLep_Pt"         )->Fill( isoLep.Pt                      , wrtevt     );
                h1.GetTH1("EvtChi2_isoLep_Et"         )->Fill( isoLep.Et                      , wrtevt     );
                h1.GetTH1("EvtChi2_isoLep_Eta"        )->Fill( isoLep.Eta                     , wrtevt     );
                h1.GetTH1("EvtChi2_isoLep_Phi"        )->Fill( isoLep.Phi                     , wrtevt     );
                h1.GetTH1("EvtChi2_isoLep_E"          )->Fill( isoLep.Energy                  , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_Chi2" )->Fill( minChi2                        , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_Mass" )->Fill( top_hadronic.Mass              , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_MassW" )->Fill( top_hadronic.MassW              , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_Pt"   )->Fill( top_hadronic.Pt                , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_Eta"  )->Fill( top_hadronic.Eta               , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Hadronic_Phi"  )->Fill( top_hadronic.Phi               , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Leptonic_Mt"   )->Fill( top_leptonic.MassT             , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"  )->Fill( top_leptonic.Massbl            , wrtevt     );
                h1.GetTH1("EvtChi2_Top_Leptonic_Phi"  )->Fill( top_leptonic.Phi               , wrtevt     );
                h1.GetTH1("EvtChi2_O2"                )->Fill( O2/Owrt_                       , wrtevt     );
                h1.GetTH1("EvtChi2_Oa"               )->Fill( Oa/Owrt_                      , wrtevt     );
                h1.GetTH1("EvtChi2_Ob"               )->Fill( Ob/Owrt_                      , wrtevt     );
                h1.GetTH1("EvtChi2_O3"                )->Fill( O3/Owrt_                       , wrtevt     );
                h1.GetTH1("EvtChi2_O4"                )->Fill( O4/Owrt_                       , wrtevt     );
                h1.GetTH1("EvtChi2_O7"                )->Fill( O7/Owrt_                       , wrtevt     );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym"), O2,  wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_ObAsym"), Ob, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_OaAsym"), Oa, wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym"), O3,  wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym"), O4,  wrtevt );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym"), O7,  wrtevt );
            }
            //// *** Muon channel
            if( isGoodMuonEvt )
            {
                //h1.GetTH1("Evt_CutFlow_Mu")->Fill("H_{T}>250GeV", wrtevt);
                b_isMuonEvt_ = 1;
                b_isEleEvt_  = 0;
                h2.GetTH2("TH2_Chi2_vs_MET_Mu"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_TopLeptonicMbl_Mu"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_TopHadronicMass_Mu"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_Ht_Mu"                       )->Fill( Ht,                  minChi2,             wrtevt );
                h2.GetTH2("TH2_TopHadronicMass_vs_Ht_Mu"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl_Mu")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2_TopLeptonicMbl_vs_MET_Mu"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
                h2.GetTH2("TH2_TopHadronicMassW_vs_dRj1j2_Mu"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
                h1.GetTH1("Evt_Ht_Mu"               )->Fill( Ht                             , wrtevt     );
                h1.GetTH1("Evt_NVertex_Mu"          )->Fill( VxtColSelected.size()          , wrtevt     );
                h1.GetTH1("Evt_NVertexNoWrt_Mu"     )->Fill( VxtColSelected.size()          , wrtevtNoPU );
                h1.GetTH1("Evt_NSelJets_Mu"         )->Fill( JetColSelected.size()          , wrtevt     );
                h1.GetTH1("Evt_NBJets_Mu"           )->Fill( BJetCol.size()                 , wrtevt     );
                h1.GetTH1("Evt_bJet_Pt_Mu"          )->Fill( b_jet.Pt                       , wrtevt     );
                h1.GetTH1("Evt_bJet_Eta_Mu"         )->Fill( b_jet.Eta                      , wrtevt     );
                h1.GetTH1("Evt_bJet_Phi_Mu"         )->Fill( b_jet.Phi                      , wrtevt     );
                h1.GetTH1("Evt_bJet_E_Mu"           )->Fill( b_jet.Energy                   , wrtevt     );
                h1.GetTH1("Evt_bJet_M_Mu"           )->Fill( b_jet.Mass                     , wrtevt     );
                h1.GetTH1("Evt_bJet_BTag_Mu"        )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Pt_Mu"       )->Fill( bbar_jet.Pt                    , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Eta_Mu"      )->Fill( bbar_jet.Eta                   , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Phi_Mu"      )->Fill( bbar_jet.Phi                   , wrtevt     );
                h1.GetTH1("Evt_bbarJet_E_Mu"        )->Fill( bbar_jet.Energy                , wrtevt     );
                h1.GetTH1("Evt_bbarJet_M_Mu"        )->Fill( bbar_jet.Mass                  , wrtevt     );
                h1.GetTH1("Evt_bbarJet_BTag_Mu"     )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
                h1.GetTH1("Evt_HardJet_Pt_Mu"       )->Fill( hardJet.Pt                     , wrtevt     );
                h1.GetTH1("Evt_HardJet_Eta_Mu"      )->Fill( hardJet.Eta                    , wrtevt     );
                h1.GetTH1("Evt_HardJet_Phi_Mu"      )->Fill( hardJet.Phi                    , wrtevt     );
                h1.GetTH1("Evt_HardJet_E_Mu"        )->Fill( hardJet.Energy                 , wrtevt     );
                h1.GetTH1("Evt_HardJet_M_Mu"        )->Fill( hardJet.Mass                   , wrtevt     );
                h1.GetTH1("Evt_HardJet_BTag_Mu"     )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Pt_Mu"  )->Fill( hardNonBJet1.Pt                , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Eta_Mu" )->Fill( hardNonBJet1.Eta               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Phi_Mu" )->Fill( hardNonBJet1.Phi               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_E_Mu"   )->Fill( hardNonBJet1.Energy            , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_M_Mu"   )->Fill( hardNonBJet1.Mass              , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_BTag_Mu")->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Pt_Mu"  )->Fill( hardNonBJet2.Pt                , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Eta_Mu" )->Fill( hardNonBJet2.Eta               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Phi_Mu" )->Fill( hardNonBJet2.Phi               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_E_Mu"   )->Fill( hardNonBJet2.Energy            , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_M_Mu"   )->Fill( hardNonBJet2.Mass              , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_BTag_Mu")->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Pt_Mu"   )->Fill( TopNonBJet1.Pt                 , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Eta_Mu"  )->Fill( TopNonBJet1.Eta                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Phi_Mu"  )->Fill( TopNonBJet1.Phi                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_E_Mu"    )->Fill( TopNonBJet1.Energy             , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_M_Mu"    )->Fill( TopNonBJet1.Mass               , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_BTag_Mu" )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Pt_Mu"   )->Fill( TopNonBJet2.Pt                 , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Eta_Mu"  )->Fill( TopNonBJet2.Eta                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Phi_Mu"  )->Fill( TopNonBJet2.Phi                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_E_Mu"    )->Fill( TopNonBJet2.Energy             , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_M_Mu"    )->Fill( TopNonBJet2.Mass               , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_BTag_Mu" )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("Evt_isoLep_Pt_Mu"        )->Fill( isoLep.Pt                      , wrtevt     );
                h1.GetTH1("Evt_isoLep_Et_Mu"        )->Fill( isoLep.Et                      , wrtevt     );
                h1.GetTH1("Evt_isoLep_Eta_Mu"       )->Fill( isoLep.Eta                     , wrtevt     );
                h1.GetTH1("Evt_isoLep_Phi_Mu"       )->Fill( isoLep.Phi                     , wrtevt     );
                h1.GetTH1("Evt_isoLep_E_Mu"         )->Fill( isoLep.Energy                  , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Chi2_Mu")->Fill( minChi2                        , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Mass_Mu")->Fill( top_hadronic.Mass              , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_MassW_Mu" )->Fill( top_hadronic.MassW              , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Pt_Mu"  )->Fill( top_hadronic.Pt                , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Eta_Mu" )->Fill( top_hadronic.Eta               , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Phi_Mu" )->Fill( top_hadronic.Phi               , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Mt_Mu"  )->Fill( top_leptonic.MassT             , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Mbl_Mu" )->Fill( top_leptonic.Massbl            , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Phi_Mu" )->Fill( top_leptonic.Phi               , wrtevt     );
                h1.GetTH1("Evt_O2_Mu"               )->Fill( O2/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_Oa_Mu"               )->Fill( Oa/Owrt_                      , wrtevt     );
                h1.GetTH1("Evt_Ob_Mu"               )->Fill( Ob/Owrt_                      , wrtevt     );
                h1.GetTH1("Evt_O3_Mu"               )->Fill( O3/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_O4_Mu"               )->Fill( O4/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_O7_Mu"               )->Fill( O7/Owrt_                       , wrtevt     );
                fillAsym( h1.GetTH1("Evt_O2Asym_Mu"),  O2, wrtevt );
                fillAsym( h1.GetTH1("Evt_OaAsym_Mu"),  Oa, wrtevt );
                fillAsym( h1.GetTH1("Evt_ObAsym_Mu"),  Ob, wrtevt );
                fillAsym( h1.GetTH1("Evt_O3Asym_Mu"),  O3, wrtevt );
                fillAsym( h1.GetTH1("Evt_O4Asym_Mu"),  O4, wrtevt );
                fillAsym( h1.GetTH1("Evt_O7Asym_Mu"),  O7, wrtevt );
                if( passChi2Cut )
                {
                    if( passMlbCut )
                    {
                        h1.GetTH1("EvtMlb_O2"   )->Fill( O2/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O2_Mu")->Fill( O2/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_Oa_Mu")->Fill( Oa/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_Ob_Mu")->Fill( Ob/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O3_Mu")->Fill( O3/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O4_Mu")->Fill( O4/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O7_Mu")->Fill( O7/Owrt_, wrtevt);
                        fillAsym( h1.GetTH1("EvtMlb_O2Asym"   ),  O2, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O2Asym_Mu"),  O2, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_OaAsym_Mu"), Oa, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_ObAsym_Mu"), Ob, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O3Asym_Mu"), O3, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O4Asym_Mu"), O4, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O7Asym_Mu"), O7, wrtevt );
                    }
                    h2.GetTH2("TH2Chi2_Chi2_vs_MET_Mu"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl_Mu"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_TopHadronicMass_Mu"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_Ht_Mu"                       )->Fill( Ht,                  minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMass_vs_Ht_Mu"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl_Mu")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
                    h2.GetTH2("TH2Chi2_TopLeptonicMbl_vs_MET_Mu"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2_Mu"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
                    h1.GetTH1("EvtChi2_Ht_Mu"               )->Fill( Ht                             , wrtevt     );
                    h1.GetTH1("EvtChi2_NVertex_Mu"          )->Fill( VxtColSelected.size()          , wrtevt     );
                    h1.GetTH1("EvtChi2_NVertexNoWrt_Mu"     )->Fill( VxtColSelected.size()          , wrtevtNoPU );
                    h1.GetTH1("EvtChi2_NSelJets_Mu"         )->Fill( JetColSelected.size()          , wrtevt     );
                    h1.GetTH1("EvtChi2_NBJets_Mu"           )->Fill( BJetCol.size()                 , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Pt_Mu"          )->Fill( b_jet.Pt                       , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Eta_Mu"         )->Fill( b_jet.Eta                      , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Phi_Mu"         )->Fill( b_jet.Phi                      , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_E_Mu"           )->Fill( b_jet.Energy                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_M_Mu"           )->Fill( b_jet.Mass                     , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_BTag_Mu"        )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Pt_Mu"       )->Fill( bbar_jet.Pt                    , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Eta_Mu"      )->Fill( bbar_jet.Eta                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Phi_Mu"      )->Fill( bbar_jet.Phi                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_E_Mu"        )->Fill( bbar_jet.Energy                , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_M_Mu"        )->Fill( bbar_jet.Mass                  , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_BTag_Mu"     )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Pt_Mu"       )->Fill( hardJet.Pt                     , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Eta_Mu"      )->Fill( hardJet.Eta                    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Phi_Mu"      )->Fill( hardJet.Phi                    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_E_Mu"        )->Fill( hardJet.Energy                 , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_M_Mu"        )->Fill( hardJet.Mass                   , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_BTag_Mu"     )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Pt_Mu"  )->Fill( hardNonBJet1.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Eta_Mu" )->Fill( hardNonBJet1.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Phi_Mu" )->Fill( hardNonBJet1.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_E_Mu"   )->Fill( hardNonBJet1.Energy            , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_M_Mu"   )->Fill( hardNonBJet1.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_BTag_Mu")->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Pt_Mu"  )->Fill( hardNonBJet2.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Eta_Mu" )->Fill( hardNonBJet2.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Phi_Mu" )->Fill( hardNonBJet2.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_E_Mu"   )->Fill( hardNonBJet2.Energy            , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_M_Mu"   )->Fill( hardNonBJet2.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_BTag_Mu")->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Pt_Mu"   )->Fill( TopNonBJet1.Pt                 , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Eta_Mu"  )->Fill( TopNonBJet1.Eta                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Phi_Mu"  )->Fill( TopNonBJet1.Phi                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_E_Mu"    )->Fill( TopNonBJet1.Energy             , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_M_Mu"    )->Fill( TopNonBJet1.Mass               , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_BTag_Mu" )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Pt_Mu"   )->Fill( TopNonBJet2.Pt                 , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Eta_Mu"  )->Fill( TopNonBJet2.Eta                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Phi_Mu"  )->Fill( TopNonBJet2.Phi                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_E_Mu"    )->Fill( TopNonBJet2.Energy             , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_M_Mu"    )->Fill( TopNonBJet2.Mass               , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_BTag_Mu" )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Pt_Mu"        )->Fill( isoLep.Pt                      , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Et_Mu"        )->Fill( isoLep.Et                      , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Eta_Mu"       )->Fill( isoLep.Eta                     , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Phi_Mu"       )->Fill( isoLep.Phi                     , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_E_Mu"         )->Fill( isoLep.Energy                  , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Chi2_Mu")->Fill( minChi2                        , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Mass_Mu")->Fill( top_hadronic.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_MassW_Mu" )->Fill( top_hadronic.MassW              , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Pt_Mu"  )->Fill( top_hadronic.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Eta_Mu" )->Fill( top_hadronic.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Phi_Mu" )->Fill( top_hadronic.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Mt_Mu"  )->Fill( top_leptonic.MassT             , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_Mu" )->Fill( top_leptonic.Massbl            , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Phi_Mu" )->Fill( top_leptonic.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_O2_Mu"               )->Fill( O2/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_Oa_Mu"              )->Fill( Oa/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_Ob_Mu"              )->Fill( Ob/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_O3_Mu"               )->Fill( O3/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_O4_Mu"               )->Fill( O4/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_O7_Mu"               )->Fill( O7/Owrt_                       , wrtevt     );
                    fillAsym( h1.GetTH1("EvtChi2_O2Asym_Mu" ), O2, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_OaAsym_Mu"), Oa, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_ObAsym_Mu"), Ob, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O3Asym_Mu" ), O3, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O4Asym_Mu" ), O4, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O7Asym_Mu" ), O7, wrtevt );
                }
            }
            //// *** Electron channel
            if( isGoodElectronEvt )
            {
                //h1.GetTH1("Evt_CutFlow_El")->Fill("H_{T}>250GeV", wrtevt);
                b_isMuonEvt_ = 0;
                b_isEleEvt_  = 1;
                h2.GetTH2("TH2_Chi2_vs_MET_El"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_TopLeptonicMbl_El"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_TopHadronicMass_El"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
                h2.GetTH2("TH2_Chi2_vs_Ht_El"                       )->Fill( Ht,                  minChi2,             wrtevt );
                h2.GetTH2("TH2_TopHadronicMass_vs_Ht_El"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2_TopHadronicMass_vs_TopLeptonicMbl_El")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
                h2.GetTH2("TH2_TopLeptonicMbl_vs_MET_El"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
                h2.GetTH2("TH2_TopHadronicMassW_vs_dRj1j2_El"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
                h1.GetTH1("Evt_Ht_El"               )->Fill( Ht                             , wrtevt     );
                h1.GetTH1("Evt_NVertex_El"          )->Fill( VxtColSelected.size()          , wrtevt     );
                h1.GetTH1("Evt_NVertexNoWrt_El"     )->Fill( VxtColSelected.size()          , wrtevtNoPU );
                h1.GetTH1("Evt_NSelJets_El"         )->Fill( JetColSelected.size()          , wrtevt     );
                h1.GetTH1("Evt_NBJets_El"           )->Fill( BJetCol.size()                 , wrtevt     );
                h1.GetTH1("Evt_bJet_Pt_El"          )->Fill( b_jet.Pt                       , wrtevt     );
                h1.GetTH1("Evt_bJet_Eta_El"         )->Fill( b_jet.Eta                      , wrtevt     );
                h1.GetTH1("Evt_bJet_Phi_El"         )->Fill( b_jet.Phi                      , wrtevt     );
                h1.GetTH1("Evt_bJet_E_El"           )->Fill( b_jet.Energy                   , wrtevt     );
                h1.GetTH1("Evt_bJet_M_El"           )->Fill( b_jet.Mass                     , wrtevt     );
                h1.GetTH1("Evt_bJet_BTag_El"        )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Pt_El"       )->Fill( bbar_jet.Pt                    , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Eta_El"      )->Fill( bbar_jet.Eta                   , wrtevt     );
                h1.GetTH1("Evt_bbarJet_Phi_El"      )->Fill( bbar_jet.Phi                   , wrtevt     );
                h1.GetTH1("Evt_bbarJet_E_El"        )->Fill( bbar_jet.Energy                , wrtevt     );
                h1.GetTH1("Evt_bbarJet_M_El"        )->Fill( bbar_jet.Mass                  , wrtevt     );
                h1.GetTH1("Evt_bbarJet_BTag_El"     )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
                h1.GetTH1("Evt_HardJet_Pt_El"       )->Fill( hardJet.Pt                     , wrtevt     );
                h1.GetTH1("Evt_HardJet_Eta_El"      )->Fill( hardJet.Eta                    , wrtevt     );
                h1.GetTH1("Evt_HardJet_Phi_El"      )->Fill( hardJet.Phi                    , wrtevt     );
                h1.GetTH1("Evt_HardJet_E_El"        )->Fill( hardJet.Energy                 , wrtevt     );
                h1.GetTH1("Evt_HardJet_M_El"        )->Fill( hardJet.Mass                   , wrtevt     );
                h1.GetTH1("Evt_HardJet_BTag_El"     )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Pt_El"  )->Fill( hardNonBJet1.Pt                , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Eta_El" )->Fill( hardNonBJet1.Eta               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_Phi_El" )->Fill( hardNonBJet1.Phi               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_E_El"   )->Fill( hardNonBJet1.Energy            , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_M_El"   )->Fill( hardNonBJet1.Mass              , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet1_BTag_El")->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Pt_El"  )->Fill( hardNonBJet2.Pt                , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Eta_El" )->Fill( hardNonBJet2.Eta               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_Phi_El" )->Fill( hardNonBJet2.Phi               , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_E_El"   )->Fill( hardNonBJet2.Energy            , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_M_El"   )->Fill( hardNonBJet2.Mass              , wrtevt     );
                h1.GetTH1("Evt_HardNonBJet2_BTag_El")->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Pt_El"   )->Fill( TopNonBJet1.Pt                 , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Eta_El"  )->Fill( TopNonBJet1.Eta                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_Phi_El"  )->Fill( TopNonBJet1.Phi                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_E_El"    )->Fill( TopNonBJet1.Energy             , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_M_El"    )->Fill( TopNonBJet1.Mass               , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet1_BTag_El" )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Pt_El"   )->Fill( TopNonBJet2.Pt                 , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Eta_El"  )->Fill( TopNonBJet2.Eta                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_Phi_El"  )->Fill( TopNonBJet2.Phi                , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_E_El"    )->Fill( TopNonBJet2.Energy             , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_M_El"    )->Fill( TopNonBJet2.Mass               , wrtevt     );
                h1.GetTH1("Evt_TopNonBJet2_BTag_El" )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
                h1.GetTH1("Evt_isoLep_Pt_El"        )->Fill( isoLep.Pt                      , wrtevt     );
                h1.GetTH1("Evt_isoLep_Et_El"        )->Fill( isoLep.Et                      , wrtevt     );
                h1.GetTH1("Evt_isoLep_Eta_El"       )->Fill( isoLep.Eta                     , wrtevt     );
                h1.GetTH1("Evt_isoLep_Phi_El"       )->Fill( isoLep.Phi                     , wrtevt     );
                h1.GetTH1("Evt_isoLep_E_El"         )->Fill( isoLep.Energy                  , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Chi2_El")->Fill( minChi2                        , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Mass_El")->Fill( top_hadronic.Mass              , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_MassW_El" )->Fill( top_hadronic.MassW              , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Pt_El"  )->Fill( top_hadronic.Pt                , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Eta_El" )->Fill( top_hadronic.Eta               , wrtevt     );
                h1.GetTH1("Evt_Top_Hadronic_Phi_El" )->Fill( top_hadronic.Phi               , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Mt_El"  )->Fill( top_leptonic.MassT             , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Mbl_El" )->Fill( top_leptonic.Massbl            , wrtevt     );
                h1.GetTH1("Evt_Top_Leptonic_Phi_El" )->Fill( top_leptonic.Phi               , wrtevt     );
                h1.GetTH1("Evt_O2_El"               )->Fill( O2/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_Oa_El"              )->Fill( Oa/Owrt_                      , wrtevt     );
                h1.GetTH1("Evt_Ob_El"              )->Fill( Ob/Owrt_                      , wrtevt     );
                h1.GetTH1("Evt_O3_El"               )->Fill( O3/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_O4_El"               )->Fill( O4/Owrt_                       , wrtevt     );
                h1.GetTH1("Evt_O7_El"               )->Fill( O7/Owrt_                       , wrtevt     );
                fillAsym( h1.GetTH1("Evt_O2Asym_El" ), O2, wrtevt );
                fillAsym( h1.GetTH1("Evt_OaAsym_El"), Oa, wrtevt );
                fillAsym( h1.GetTH1("Evt_ObAsym_El"), Ob, wrtevt );
                fillAsym( h1.GetTH1("Evt_O3Asym_El" ), O3, wrtevt );
                fillAsym( h1.GetTH1("Evt_O4Asym_El" ), O4, wrtevt );
                fillAsym( h1.GetTH1("Evt_O7Asym_El" ), O7, wrtevt );
                if( passChi2Cut )
                {
                    if( passMlbCut )
                    {
                        h1.GetTH1("EvtMlb_O2"   )->Fill( O2/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O2_El")->Fill( O2/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_Oa_El")->Fill( Oa/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_Ob_El")->Fill( Ob/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O3_El")->Fill( O3/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O4_El")->Fill( O4/Owrt_, wrtevt);
                        h1.GetTH1("EvtMlb_O7_El")->Fill( O7/Owrt_, wrtevt);
                        fillAsym( h1.GetTH1("EvtMlb_O2Asym"   ), O2, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O2Asym_El"), O2, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_OaAsym_El"), Oa, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_ObAsym_El"), Ob, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O3Asym_El"), O3, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O4Asym_El"), O4, wrtevt );
                        fillAsym( h1.GetTH1("EvtMlb_O7Asym_El"), O7, wrtevt );
                    }
                    h2.GetTH2("TH2Chi2_Chi2_vs_MET_El"                      )->Fill( EvtInfo.PFMET,       minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_TopLeptonicMbl_El"           )->Fill( top_leptonic.Massbl, minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_TopHadronicMass_El"          )->Fill( top_hadronic.Mass,   minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_Chi2_vs_Ht_El"                       )->Fill( Ht,                  minChi2,             wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMass_vs_Ht_El"            )->Fill( Ht,                  top_hadronic.Mass,   wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMass_vs_TopLeptonicMbl_El")->Fill( top_leptonic.Massbl, top_hadronic.Mass,   wrtevt );
                    h2.GetTH2("TH2Chi2_TopLeptonicMbl_vs_MET_El"            )->Fill( EvtInfo.PFMET,       top_leptonic.Massbl, wrtevt );
                    h2.GetTH2("TH2Chi2_TopHadronicMassW_vs_dRj1j2_El"       )->Fill( top_hadronic.MassW,  top_hadronic.dRj1j2, wrtevt );
                    h1.GetTH1("EvtChi2_Ht_El"               )->Fill( Ht                             , wrtevt     );
                    h1.GetTH1("EvtChi2_NVertex_El"          )->Fill( VxtColSelected.size()          , wrtevt     );
                    h1.GetTH1("EvtChi2_NVertexNoWrt_El"     )->Fill( VxtColSelected.size()          , wrtevtNoPU );
                    h1.GetTH1("EvtChi2_NSelJets_El"         )->Fill( JetColSelected.size()          , wrtevt     );
                    h1.GetTH1("EvtChi2_NBJets_El"           )->Fill( BJetCol.size()                 , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Pt_El"          )->Fill( b_jet.Pt                       , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Eta_El"         )->Fill( b_jet.Eta                      , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_Phi_El"         )->Fill( b_jet.Phi                      , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_E_El"           )->Fill( b_jet.Energy                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_M_El"           )->Fill( b_jet.Mass                     , wrtevt     );
                    h1.GetTH1("EvtChi2_bJet_BTag_El"        )->Fill( b_jet.CombinedSVBJetTags       , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Pt_El"       )->Fill( bbar_jet.Pt                    , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Eta_El"      )->Fill( bbar_jet.Eta                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_Phi_El"      )->Fill( bbar_jet.Phi                   , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_E_El"        )->Fill( bbar_jet.Energy                , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_M_El"        )->Fill( bbar_jet.Mass                  , wrtevt     );
                    h1.GetTH1("EvtChi2_bbarJet_BTag_El"     )->Fill( bbar_jet.CombinedSVBJetTags    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Pt_El"       )->Fill( hardJet.Pt                     , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Eta_El"      )->Fill( hardJet.Eta                    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_Phi_El"      )->Fill( hardJet.Phi                    , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_E_El"        )->Fill( hardJet.Energy                 , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_M_El"        )->Fill( hardJet.Mass                   , wrtevt     );
                    h1.GetTH1("EvtChi2_HardJet_BTag_El"     )->Fill( hardJet.CombinedSVBJetTags     , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Pt_El"  )->Fill( hardNonBJet1.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Eta_El" )->Fill( hardNonBJet1.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_Phi_El" )->Fill( hardNonBJet1.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_E_El"   )->Fill( hardNonBJet1.Energy            , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_M_El"   )->Fill( hardNonBJet1.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet1_BTag_El")->Fill( hardNonBJet1.CombinedSVBJetTags, wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Pt_El"  )->Fill( hardNonBJet2.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Eta_El" )->Fill( hardNonBJet2.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_Phi_El" )->Fill( hardNonBJet2.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_E_El"   )->Fill( hardNonBJet2.Energy            , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_M_El"   )->Fill( hardNonBJet2.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_HardNonBJet2_BTag_El")->Fill( hardNonBJet2.CombinedSVBJetTags, wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Pt_El"   )->Fill( TopNonBJet1.Pt                 , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Eta_El"  )->Fill( TopNonBJet1.Eta                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_Phi_El"  )->Fill( TopNonBJet1.Phi                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_E_El"    )->Fill( TopNonBJet1.Energy             , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_M_El"    )->Fill( TopNonBJet1.Mass               , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet1_BTag_El" )->Fill( TopNonBJet1.CombinedSVBJetTags , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Pt_El"   )->Fill( TopNonBJet2.Pt                 , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Eta_El"  )->Fill( TopNonBJet2.Eta                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_Phi_El"  )->Fill( TopNonBJet2.Phi                , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_E_El"    )->Fill( TopNonBJet2.Energy             , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_M_El"    )->Fill( TopNonBJet2.Mass               , wrtevt     );
                    h1.GetTH1("EvtChi2_TopNonBJet2_BTag_El" )->Fill( TopNonBJet2.CombinedSVBJetTags , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Pt_El"        )->Fill( isoLep.Pt                      , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Et_El"        )->Fill( isoLep.Et                      , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Eta_El"       )->Fill( isoLep.Eta                     , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_Phi_El"       )->Fill( isoLep.Phi                     , wrtevt     );
                    h1.GetTH1("EvtChi2_isoLep_E_El"         )->Fill( isoLep.Energy                  , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Chi2_El")->Fill( minChi2                        , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Mass_El")->Fill( top_hadronic.Mass              , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_MassW_El" )->Fill( top_hadronic.MassW              , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Pt_El"  )->Fill( top_hadronic.Pt                , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Eta_El" )->Fill( top_hadronic.Eta               , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Hadronic_Phi_El" )->Fill( top_hadronic.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Mt_El"  )->Fill( top_leptonic.MassT             , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_El" )->Fill( top_leptonic.Massbl            , wrtevt     );
                    h1.GetTH1("EvtChi2_Top_Leptonic_Phi_El" )->Fill( top_leptonic.Phi               , wrtevt     );
                    h1.GetTH1("EvtChi2_O2_El"               )->Fill( O2/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_Oa_El"              )->Fill( Oa/Owrt_                      , wrtevt     );
                    h1.GetTH1("EvtChi2_Ob_El"              )->Fill( Ob/Owrt_                      , wrtevt     );
                    h1.GetTH1("EvtChi2_O3_El"               )->Fill( O3/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_O4_El"               )->Fill( O4/Owrt_                       , wrtevt     );
                    h1.GetTH1("EvtChi2_O7_El"               )->Fill( O7/Owrt_                       , wrtevt     );
                    fillAsym( h1.GetTH1("EvtChi2_O2Asym_El" ), O2, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_OaAsym_El"), Oa, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_ObAsym_El"), Ob, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O3Asym_El" ), O3, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O4Asym_El" ), O4, wrtevt );
                    fillAsym( h1.GetTH1("EvtChi2_O7Asym_El" ), O7, wrtevt );
                }
            }
            newtree_->Fill();
        } //// [END] Fill plots
        if( doPDFTree_ )
        {
            //if( passChi2Cut && passMlbCut )
            if( passChi2Cut )
            { 
                b_TopMlb_=top_leptonic.Massbl;
                b_WrtEvt_=wrtevt;
            }
            else
            {
                b_TopMlb_=0;
                b_WrtEvt_=0;
            }
            if( isGoodElectronEvt && !isGoodMuonEvt )
            { 
                b_isMuonEvt_ = 0;
                b_isEleEvt_  = 1;
            }
            else if( !isGoodElectronEvt && isGoodMuonEvt )
            {
                b_isMuonEvt_ = 1;
                b_isEleEvt_  = 0;
            }
            else
            {
                b_isMuonEvt_ = 0;
                b_isEleEvt_  = 0;
            }
            if( !isdata && genTopPt>=0. && genAntiTopPt>=0. ) b_isSignal_=1;
            else b_isSignal_=0;
            b_O2_ = O2;
            b_O3_ = O3;
            b_O4_ = O4;
            b_O7_ = O7;
            b_Ob_ = Ob;
            b_Oa_ = Oa;
            b_minChi2_ = minChi2;
            b_WrtObs_  = Owrt_;
            pdftree_->Fill();
        }
    } //// [END] entry loop 
} //// [END] Class

// ------------ method called once each job just after ending the event loop  ------------
void SemiLeptanicAnalysis::endJob(){
    std::cout<<">> [INFO] End of Job!"<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SemiLeptanicAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SemiLeptanicAnalysis);
