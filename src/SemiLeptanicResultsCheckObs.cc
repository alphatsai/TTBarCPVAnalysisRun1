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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/checkEvtTool.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/GenParticle.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicResultsCheckObs.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorElectron.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorMuon.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorJet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorVertex.h" 

//
// constructors and destructor
//
SemiLeptanicResultsCheckObs::SemiLeptanicResultsCheckObs(const edm::ParameterSet& iConfig) : 
    maxEvents_(   iConfig.getParameter<int>("MaxEvents")), 
    reportEvery_( iConfig.getParameter<int>("ReportEvery")),
    inputTTree_(  iConfig.getParameter<std::string>("InputTTree")),
    inputFiles_(  iConfig.getParameter<std::vector<std::string> >("InputFiles")),
    maxChi2Cut_(  iConfig.getParameter<double>("MaxChi2Cut")),
    minChi2Cut_(  iConfig.getParameter<double>("MinChi2Cut")),
    Owrt_(        iConfig.getParameter<double>("Owrt")),
    GenACP_(      iConfig.getParameter<double>("GenACP")),
    doWrtEvt_(    iConfig.getParameter<bool>("DoWrtEvt")),
    Debug_(       iConfig.getParameter<bool>("Debug"))
{}

SemiLeptanicResultsCheckObs::~SemiLeptanicResultsCheckObs()
{ 
    delete chain_;
}

// ------------ Other function -------------
void SemiLeptanicResultsCheckObs::checkObsChange( TH1D* h, double recoO, double genO, double wrt )
{
    if( recoO > 0 && genO > 0 ) h->Fill(0.,wrt);
    if( recoO < 0 && genO < 0 ) h->Fill(1.,wrt);
    if( recoO > 0 && genO < 0 ) h->Fill(2.,wrt);
    if( recoO < 0 && genO > 0 ) h->Fill(3.,wrt);
}
// ------------ method called once each job just before starting event loop  ------------
void SemiLeptanicResultsCheckObs::beginJob()
{
    // Create TH1D
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.ClearTH1Info();
    h1.addNewTH1( "Evt_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Evt_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
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
    h1.addNewTH1( "Evt_ChargeDiff",          "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "Evt_ChargeDiff_El",       "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "Evt_ChargeDiff_Mu",       "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "EvtChi2_ChargeDiff",      "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "EvtChi2_ChargeDiff_El",   "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "EvtChi2_ChargeDiff_Mu",   "",                        "",                  "Events", "", "",  3,  -1,   2   );
    h1.addNewTH1( "Evt_ChangeO2",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO3",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO4",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO7",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO2_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO3_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO4_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO7_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO2_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO3_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO4_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_ChangeO7_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO2",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO3",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO4",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO7",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO2_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO3_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO4_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO7_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO2_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO3_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO4_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtChi2_ChangeO7_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );

    h1.addNewTH1( "Gen_Top_Mass",            "",                    "Mass",                  "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Gen_W_Mass",              "",                    "Mass",                  "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Gen_PID_DiffLep_El",      "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_DiffLep_Mu",      "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Num_wJets",               "",                        "",                  "Events", "", "",  10,   0,   10   );
    h1.addNewTH1( "Num_wLeps",               "",                        "",                  "Events", "", "",  10,   0,   10   );
    h1.addNewTH1( "GenChi2_PID_Lep_Mu",      "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "GenChi2_PID_Lep_El",      "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_Lep_Mu",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_Lep_El",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID",                 "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo1",            "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo2",            "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo1",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo2",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo1",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo2",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo1",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo2",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa1",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa2",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa1",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa2",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wJets",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wLeps",           "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_O2",                  "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O2_Mu",               "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O2_El",               "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O3",                  "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O3_Mu",               "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O3_El",               "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O4",                  "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O4_Mu",               "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O4_El",               "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O7",                  "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O7_Mu",               "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O7_El",               "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "Gen_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "Gen_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O2",              "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O2_Mu",           "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O2_El",           "O2",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O3",              "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O3_Mu",           "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O3_El",           "O3",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O4",              "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O4_Mu",           "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O4_El",           "O4",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O7",              "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O7_Mu",           "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O7_El",           "O7",                      "",                  "Events", "", "", 100,  -5,   5  );
    h1.addNewTH1( "GenChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );
    h1.addNewTH1( "GenChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1  );

    h1.CreateTH1( fs );
    h1.Sumw2();

    // Create TH2D
    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.ClearTH2Info();
    h2.addNewTH2("TH2_Gen_vs_RECO_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_RECO_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_RECO_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.CreateTH2( fs );
    h2.Sumw2();

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
    chain_  = new TChain(inputTTree_.c_str());

    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
    }

    //VtxInfo.Register(chain_);
    EvtInfo.Register(chain_);
    GenInfo.Register(chain_);
    JetInfo.Register(chain_,"PFJetInfo");
    //LepInfo.Register(chain_,"PFLepInfo");
    IsoLepInfo.Register(chain_,"isoLepton");

    chain_->SetBranchAddress("EvtInfo.O2",        &O2_        );
    chain_->SetBranchAddress("EvtInfo.O3",        &O3_        );
    chain_->SetBranchAddress("EvtInfo.O4",        &O4_        );
    chain_->SetBranchAddress("EvtInfo.O7",        &O7_        );
    chain_->SetBranchAddress("EvtInfo.MinChi2",   &minChi2_   );
    chain_->SetBranchAddress("EvtInfo.WrtObs",    &WrtObs_    );
    chain_->SetBranchAddress("EvtInfo.WrtEvt",    &WrtEvt_    );
    chain_->SetBranchAddress("EvtInfo.isMuonEvt", &isMuonEvt_ );
    chain_->SetBranchAddress("EvtInfo.isEleEvt",  &isEleEvt_  );

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void SemiLeptanicResultsCheckObs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
    using namespace edm;
    using namespace std;

    if(  chain_ == 0 ) return;

    cout<<">> [INFO] Starting analysis loop with "<<maxEvents_<<" events..."<<endl;

    TVector3 ax, ay, az;
    ax.SetXYZ(1, 0, 0);
    ay.SetXYZ(0, 1, 0);
    az.SetXYZ(0, 0, 1);

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
        double wrtobs = WrtObs_;
        if( !doWrtEvt_ ){ wrtevt=1; };
        // ---- * Checking generator info
        //int lepWpid=0;
        vector<GenParticle> particles;
        vector<GenParticle> wJets;
        vector<GenParticle> wLeps;
        GenParticle  b_quark, bbar_quark, lepton, j1, t_quark, tbar_quark;
        for( int idx=0; idx < GenInfo.Size; idx++ )
        {
            GenParticle particle( GenInfo, idx );
            if( particle.Status == 3 )
            { 
                particles.push_back( particle );
                h1.GetTH1("Gen_PID")->Fill( particle.PdgID );  
                if( particle.PdgID == -24 )
                {
                    
                    h1.GetTH1("Gen_W_Mass")->Fill( particle.Mass );   
                    h1.GetTH1("Gen_PID_wpDa1")->Fill( particle.Da1PdgID );   
                    h1.GetTH1("Gen_PID_wpDa2")->Fill( particle.Da2PdgID );   
                }
                if( particle.PdgID == 24 )
                {
                    h1.GetTH1("Gen_W_Mass")->Fill( particle.Mass );   
                    h1.GetTH1("Gen_PID_wmDa1")->Fill( particle.Da1PdgID );   
                    h1.GetTH1("Gen_PID_wmDa2")->Fill( particle.Da2PdgID );   
                }
                if( abs(particle.PdgID) == 6 )
                {
                    if( particle.PdgID > 0 ) t_quark    = particle; 
                    if( particle.PdgID < 0 ) tbar_quark = particle; 
                    h1.GetTH1("Gen_Top_Mass")->Fill( particle.Mass );   
                }
                if( abs(particle.PdgID) == 5 )
                { 
                    h1.GetTH1("Gen_PID_bMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_bMo2")->Fill( particle.Mo2PdgID );  
                    if( particle.PdgID ==  5 )    b_quark = particle; 
                    if( particle.PdgID == -5 ) bbar_quark = particle; 
                }
                if( abs(particle.PdgID)  < 5 || particle.PdgID == 21 )
                { 
                    h1.GetTH1("Gen_PID_lqMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lqMo2")->Fill( particle.Mo2PdgID );  
                    if( abs(particle.Mo1PdgID) == 24 || abs(particle.Mo2PdgID) == 24 )
                    {
                        wJets.push_back(particle);
                        h1.GetTH1("Gen_PID_wJets")->Fill( particle.PdgID );  
                    } 
                }
                if( abs(particle.PdgID) == 11 || abs(particle.PdgID) == 13 || abs(particle.PdgID) == 15)
                {
                    if( abs(particle.Mo1PdgID) == 24 || abs(particle.Mo2PdgID) == 24 )
                    {
                        //lepWpid = (abs(particle.Mo1PdgID)==24) ? particle.Mo1PdgID:particle.Mo2PdgID; 
                        wLeps.push_back(particle);        
                        h1.GetTH1("Gen_PID_wLeps")->Fill( particle.PdgID );  
                    }
                    h1.GetTH1("Gen_PID_lepMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lepMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) == 12 || abs(particle.PdgID) == 14 || abs(particle.PdgID) == 16 )
                { 
                    h1.GetTH1("Gen_PID_nuMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_nuMo2")->Fill( particle.Mo2PdgID );   
                }
            }
        }
        h1.GetTH1("Num_wJets")->Fill(wJets.size());
        h1.GetTH1("Num_wLeps")->Fill(wLeps.size());

        //// * RECO - GEN Jet matching GEN Particle
        //GenParticle  b_quark_Htop;
        //if( lepWpid ==  24 ) b_quark_Htop = b_quark;
        //if( lepWpid == -24 ) b_quark_Htop = bbar_quark;

        //int ibj=-1;
        //int ij1=-1;
        //int ij2=-1;
        //float dRbb=100000;
        //for( int i=0; i<JetInfo.Size; i++ )
        //{
        //    TVector3 P3;
        //    P3.SetPtEtaPhi( JetInfo.GenPt[i], JetInfo.GenEta[i], JetInfo.GenPhi[i] );
        //    if( JetInfo.GenPdgID[i] == b_quark_Htop.PdgID )
        //    {
        //        float tmp=b_quark_Htop.P3.DeltaR(P3);
        //        if( tmp < dRbb )
        //        {
        //            dRbb=tmp; 
        //            ibj=i;
        //        }
        //    }
        //}

        // * GEN ACP
        int charge=0;
        if( wJets.size() == 2 )
        {   
            int idx=-1;
            float pt_=0.;
            for( int i=0; i<int(wJets.size()); i++ )
            {
                if( wJets[i].Pt >= pt_ )
                {
                    idx = i;
                    pt_ = wJets[i].Pt;
                }
            }
            j1=wJets[idx];

        }else{ std::cout<<"[WARNING] More or less than 2 jets from W: "<<wJets.size()<<std::endl; }
        if( wLeps.size() == 1 )
        {
                lepton = wLeps[0];
                if( lepton.PdgID < 0 )      charge=1;
                else if( lepton.PdgID > 0 ) charge=-1;
        }else{ std::cout<<"[WARNING] More or less than 1 jets from W: "<<wLeps.size()<<std::endl; }

        double GenO2(0), GenO3(0), GenO4(0), GenO7(0); 
        GenO2 = Obs2( lepton.P3, j1.P3, b_quark.P3, bbar_quark.P3 );
        GenO3 = Obs3( lepton.P4, j1.P4, b_quark.P4, bbar_quark.P4, charge );
        GenO4 = Obs4( lepton.P3, j1.P3, b_quark.P3, bbar_quark.P3, charge );
        GenO7 = Obs7( az, b_quark.P3, bbar_quark.P3 );

        // ---- * ACP INFO
        double acpWrtO2 = (GenO2>0)? (1+GenACP_):(1-GenACP_);  
        double acpWrtO3 = (GenO3>0)? (1+GenACP_):(1-GenACP_);  
        double acpWrtO4 = (GenO4>0)? (1+GenACP_):(1-GenACP_);  
        double acpWrtO7 = (GenO7>0)? (1+GenACP_):(1-GenACP_); 

        checkObsChange( h1.GetTH1("Evt_ChangeO2"), O2_, GenO2, wrtevt );
        checkObsChange( h1.GetTH1("Evt_ChangeO3"), O3_, GenO3, wrtevt );
        checkObsChange( h1.GetTH1("Evt_ChangeO4"), O4_, GenO4, wrtevt );
        checkObsChange( h1.GetTH1("Evt_ChangeO7"), O7_, GenO7, wrtevt );
        h1.GetTH1("Evt_ChargeDiff")->Fill( IsoLepInfo.Charge*charge );
        h1.GetTH1("Evt_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
        h1.GetTH1("Evt_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
        h1.GetTH1("Evt_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
        h1.GetTH1("Evt_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
        h1.GetTH1("Gen_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
        h1.GetTH1("Gen_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
        h1.GetTH1("Gen_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
        h1.GetTH1("Gen_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
        fillAsym( h1.GetTH1("Evt_O2Asym"), O2_, wrtevt*acpWrtO2 );
        fillAsym( h1.GetTH1("Evt_O3Asym"), O3_, wrtevt*acpWrtO3 );
        fillAsym( h1.GetTH1("Evt_O4Asym"), O4_, wrtevt*acpWrtO4 );
        fillAsym( h1.GetTH1("Evt_O7Asym"), O7_, wrtevt*acpWrtO7 );
        fillAsym( h1.GetTH1("Gen_O2Asym"), GenO2, wrtevt*acpWrtO2 );
        fillAsym( h1.GetTH1("Gen_O3Asym"), GenO3, wrtevt*acpWrtO3 );
        fillAsym( h1.GetTH1("Gen_O4Asym"), GenO4, wrtevt*acpWrtO4 );
        fillAsym( h1.GetTH1("Gen_O7Asym"), GenO7, wrtevt*acpWrtO7 );
        h2.GetTH2("TH2_Gen_vs_RECO_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
        h2.GetTH2("TH2_Gen_vs_RECO_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
        h2.GetTH2("TH2_Gen_vs_RECO_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
        h2.GetTH2("TH2_Gen_vs_RECO_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
        if( isMuonEvt_ && !isEleEvt_ )
        {
            if( IsoLepInfo.Charge*charge < 0 ) h1.GetTH1("Gen_PID_DiffLep_Mu")->Fill(lepton.PdgID);
            checkObsChange( h1.GetTH1("Evt_ChangeO2_Mu"), O2_, GenO2, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO3_Mu"), O3_, GenO3, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO4_Mu"), O4_, GenO4, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO7_Mu"), O7_, GenO7, wrtevt );
            h1.GetTH1("Evt_ChargeDiff_Mu")->Fill( IsoLepInfo.Charge*charge );
            h1.GetTH1("Gen_PID_Lep_Mu")->Fill( lepton.PdgID );   
            h1.GetTH1("Evt_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
            h1.GetTH1("Evt_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
            h1.GetTH1("Evt_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
            h1.GetTH1("Evt_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
            h1.GetTH1("Gen_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
            h1.GetTH1("Gen_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
            h1.GetTH1("Gen_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
            h1.GetTH1("Gen_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
            fillAsym( h1.GetTH1("Evt_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
            fillAsym( h1.GetTH1("Evt_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
            fillAsym( h1.GetTH1("Evt_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
            fillAsym( h1.GetTH1("Evt_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
            fillAsym( h1.GetTH1("Gen_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
            fillAsym( h1.GetTH1("Gen_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
            fillAsym( h1.GetTH1("Gen_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
            fillAsym( h1.GetTH1("Gen_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
            h2.GetTH2("TH2_Gen_vs_RECO_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
            h2.GetTH2("TH2_Gen_vs_RECO_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
            h2.GetTH2("TH2_Gen_vs_RECO_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
            h2.GetTH2("TH2_Gen_vs_RECO_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
          
            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO2"),    O2_, GenO2, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO3"),    O3_, GenO3, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO4"),    O4_, GenO4, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO7"),    O7_, GenO7, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO2_Mu"), O2_, GenO2, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO3_Mu"), O3_, GenO3, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO4_Mu"), O4_, GenO4, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO7_Mu"), O7_, GenO7, wrtevt );
                h1.GetTH1("EvtChi2_ChargeDiff")->Fill( IsoLepInfo.Charge*charge );
                h1.GetTH1("EvtChi2_ChargeDiff_Mu")->Fill( IsoLepInfo.Charge*charge );
                h1.GetTH1("GenChi2_PID_Lep_Mu")->Fill( lepton.PdgID );   
                h1.GetTH1("EvtChi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("EvtChi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("EvtChi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("EvtChi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                h1.GetTH1("GenChi2_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("GenChi2_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("GenChi2_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("GenChi2_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                h1.GetTH1("EvtChi2_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("EvtChi2_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("EvtChi2_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("EvtChi2_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                h1.GetTH1("GenChi2_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("GenChi2_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("GenChi2_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("GenChi2_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym"),    O2_, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym"),    O3_, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym"),    O4_, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym"),    O7_, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("GenChi2_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("GenChi2_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("GenChi2_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("GenChi2_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("GenChi2_O2Asym"),    GenO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("GenChi2_O3Asym"),    GenO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("GenChi2_O4Asym"),    GenO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("GenChi2_O7Asym"),    GenO7, wrtevt*acpWrtO7 );
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
            }
        }
        else if( !isMuonEvt_ &&  isEleEvt_ )
        {
            if( IsoLepInfo.Charge*charge < 0 ) h1.GetTH1("Gen_PID_DiffLep_El")->Fill(lepton.PdgID);
            checkObsChange( h1.GetTH1("Evt_ChangeO2_El"), O2_, GenO2, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO3_El"), O3_, GenO3, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO4_El"), O4_, GenO4, wrtevt );
            checkObsChange( h1.GetTH1("Evt_ChangeO7_El"), O7_, GenO7, wrtevt );
            h1.GetTH1("Evt_ChargeDiff_El")->Fill( IsoLepInfo.Charge*charge );
            h1.GetTH1("Gen_PID_Lep_El")->Fill( lepton.PdgID );   
            h1.GetTH1("Evt_O2_El")->Fill( O2_/wrtobs,    wrtevt*acpWrtO2 );
            h1.GetTH1("Evt_O3_El")->Fill( O3_/wrtobs,    wrtevt*acpWrtO3 );
            h1.GetTH1("Evt_O4_El")->Fill( O4_/wrtobs,    wrtevt*acpWrtO4 );
            h1.GetTH1("Evt_O7_El")->Fill( O7_/wrtobs,    wrtevt*acpWrtO7 );
            h1.GetTH1("Gen_O2_El")->Fill( GenO2/wrtobs,  wrtevt*acpWrtO2 );
            h1.GetTH1("Gen_O3_El")->Fill( GenO3/wrtobs,  wrtevt*acpWrtO3 );
            h1.GetTH1("Gen_O4_El")->Fill( GenO4/wrtobs,  wrtevt*acpWrtO4 );
            h1.GetTH1("Gen_O7_El")->Fill( GenO7/wrtobs,  wrtevt*acpWrtO7 );
            fillAsym( h1.GetTH1("Evt_O2Asym_El"), O2_,   wrtevt*acpWrtO2 );
            fillAsym( h1.GetTH1("Evt_O3Asym_El"), O3_,   wrtevt*acpWrtO3 );
            fillAsym( h1.GetTH1("Evt_O4Asym_El"), O4_,   wrtevt*acpWrtO4 );
            fillAsym( h1.GetTH1("Evt_O7Asym_El"), O7_,   wrtevt*acpWrtO7 );
            fillAsym( h1.GetTH1("Gen_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
            fillAsym( h1.GetTH1("Gen_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
            fillAsym( h1.GetTH1("Gen_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
            fillAsym( h1.GetTH1("Gen_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
            h2.GetTH2("TH2_Gen_vs_RECO_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
            h2.GetTH2("TH2_Gen_vs_RECO_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
            h2.GetTH2("TH2_Gen_vs_RECO_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
            h2.GetTH2("TH2_Gen_vs_RECO_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO2"),    O2_, GenO2, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO3"),    O3_, GenO3, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO4"),    O4_, GenO4, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO7"),    O7_, GenO7, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO2_El"), O2_, GenO2, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO3_El"), O3_, GenO3, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO4_El"), O4_, GenO4, wrtevt );
                checkObsChange( h1.GetTH1("EvtChi2_ChangeO7_El"), O7_, GenO7, wrtevt );
                h1.GetTH1("EvtChi2_ChargeDiff")->Fill( IsoLepInfo.Charge*charge );
                h1.GetTH1("EvtChi2_ChargeDiff_El")->Fill( IsoLepInfo.Charge*charge );
                h1.GetTH1("GenChi2_PID_Lep_El")->Fill( lepton.PdgID );   
                h1.GetTH1("EvtChi2_O2")->Fill( O2_/wrtobs,      wrtevt*acpWrtO2 );
                h1.GetTH1("EvtChi2_O3")->Fill( O3_/wrtobs,      wrtevt*acpWrtO3 );
                h1.GetTH1("EvtChi2_O4")->Fill( O4_/wrtobs,      wrtevt*acpWrtO4 );
                h1.GetTH1("EvtChi2_O7")->Fill( O7_/wrtobs,      wrtevt*acpWrtO7 );
                h1.GetTH1("GenChi2_O2")->Fill( GenO2/wrtobs,    wrtevt*acpWrtO2 );
                h1.GetTH1("GenChi2_O3")->Fill( GenO3/wrtobs,    wrtevt*acpWrtO3 );
                h1.GetTH1("GenChi2_O4")->Fill( GenO4/wrtobs,    wrtevt*acpWrtO4 );
                h1.GetTH1("GenChi2_O7")->Fill( GenO7/wrtobs,    wrtevt*acpWrtO7 );
                h1.GetTH1("EvtChi2_O2_El")->Fill( O2_/wrtobs,   wrtevt*acpWrtO2 );
                h1.GetTH1("EvtChi2_O3_El")->Fill( O3_/wrtobs,   wrtevt*acpWrtO3 );
                h1.GetTH1("EvtChi2_O4_El")->Fill( O4_/wrtobs,   wrtevt*acpWrtO4 );
                h1.GetTH1("EvtChi2_O7_El")->Fill( O7_/wrtobs,   wrtevt*acpWrtO7 );
                h1.GetTH1("GenChi2_O2_El")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("GenChi2_O3_El")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("GenChi2_O4_El")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("GenChi2_O7_El")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym_El"), O2_,   wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym_El"), O3_,   wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym_El"), O4_,   wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym_El"), O7_,   wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtChi2_O2Asym"),    O2_,   wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtChi2_O3Asym"),    O3_,   wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtChi2_O4Asym"),    O4_,   wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtChi2_O7Asym"),    O7_,   wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("GenChi2_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("GenChi2_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("GenChi2_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("GenChi2_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("GenChi2_O2Asym"),    GenO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("GenChi2_O3Asym"),    GenO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("GenChi2_O4Asym"),    GenO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("GenChi2_O7Asym"),    GenO7, wrtevt*acpWrtO7 );
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                h2.GetTH2("TH2Chi2_Gen_vs_RECO_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
            }
        }

    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void SemiLeptanicResultsCheckObs::endJob(){
    std::cout<<">> [INFO] End of Job!"<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SemiLeptanicResultsCheckObs::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SemiLeptanicResultsCheckObs);
