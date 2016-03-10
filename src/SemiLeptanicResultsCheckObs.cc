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
    maxMlbCut_(   iConfig.getParameter<double>("MaxMlbCut")),
    minMlbCut_(   iConfig.getParameter<double>("MinMlbCut")),
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
bool SemiLeptanicResultsCheckObs::checkObsChange( TH1D* h, double recoO, double genO, double wrt, std::string name )
{
    bool isGot=true; 
    if( recoO > 0 && genO > 0 ){ h->Fill(0.,wrt); }
    else if( recoO < 0 && genO < 0 ){ h->Fill(1.,wrt); }
    else if( recoO > 0 && genO < 0 ){ h->Fill(2.,wrt); }
    else if( recoO < 0 && genO > 0 ){ h->Fill(3.,wrt); }
    else{ isGot=false; std::cout<<">>[WARNING] No match "<<name<<", RECO "<<recoO<<", GEN "<<genO<<", wrt "<<wrt<<endl; }
    return isGot; 
}
// ------------ method called once each job just before starting event loop  ------------
void SemiLeptanicResultsCheckObs::beginJob()
{
    // Create TH1D
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.ClearTH1Info();
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl",            "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_El",         "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2_Top_Leptonic_Mbl_Mu",         "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2Matched_Top_Leptonic_Mbl",     "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2Matched_Top_Leptonic_Mbl_El",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EvtChi2Matched_Top_Leptonic_Mbl_Mu",  "",                          "Mass",              "Events", "",    "", 500,  0,   500 );
    h1.addNewTH1( "EffMatched_Top_Leptonic_Mbl",         "",                          "Mass",              "Events", "",    "",  50,  0,   500 );
    h1.addNewTH1( "EffMatched_Top_Leptonic_Mbl_El",      "",                          "Mass",              "Events", "",    "",  50,  0,   500 );
    h1.addNewTH1( "EffMatched_Top_Leptonic_Mbl_Mu",      "",                          "Mass",              "Events", "",    "",  50,  0,   500 );

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

    h1.addNewTH1( "EvtPar_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtPar_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtPar_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtParChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtParChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_PID_b",                 "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Par_PID_bbar",              "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Par_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "Par_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "ParChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "ParChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Par_ChangeO2",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO3",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO4",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO7",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO2_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO3_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO4_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO7_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO2_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO3_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO4_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Par_ChangeO7_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO2",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO3",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO4",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO7",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO2_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO3_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO4_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO7_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO2_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO3_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO4_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "ParChi2_ChangeO7_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );

    h1.addNewTH1( "EvtBB_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBB_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtBBChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBBChi2_Flavor_j1",       "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "EvtBBChi2_MatchedPID_j1q1", "",                        "",                  "Evetns", "", "",  2,   0,   1 );
    h1.addNewTH1( "EvtBBChi2_dR_j1q1",         "",                        "",                  "Evetns", "", "", 1000, 0,   10 );
    h1.addNewTH1( "GenBB_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBB_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBB_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenBBChi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenBBChi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtBB_ChangeO2",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO3",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO4",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO7",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO2_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO3_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO4_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO7_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO2_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO3_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO4_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBB_ChangeO7_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO2",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO3",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO4",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO7",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO2_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO3_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO4_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO7_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO2_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO3_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO4_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtBBChi2_ChangeO7_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );

    h1.addNewTH1( "EvtJ1Q1_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1Chi2_Flavor_j1",       "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "EvtJ1Q1Chi2_MatchedPID_j1q1", "",                        "",                  "Evetns", "", "",  2,   0,   1 );
    h1.addNewTH1( "EvtJ1Q1Chi2_dR_j1q1",         "",                        "",                  "Evetns", "", "", 1000, 0,   10 );
    h1.addNewTH1( "GenJ1Q1_O7",                  "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O7_Mu",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O7_El",               "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O7Asym",              "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O7Asym_Mu",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O7Asym_El",           "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O4",                  "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O4_Mu",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O4_El",               "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O4Asym",              "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O4Asym_Mu",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O4Asym_El",           "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O3",                  "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O3_Mu",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O3_El",               "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O3Asym",              "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O3Asym_Mu",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O3Asym_El",           "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O2",                  "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O2_Mu",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O2_El",               "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1_O2Asym",              "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O2Asym_Mu",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1_O2Asym_El",           "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7",              "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7_Mu",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7_El",           "O7",                      "O_{7}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7Asym",          "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7Asym_Mu",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O7Asym_El",       "A_{O7}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4",              "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4_Mu",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4_El",           "O4",                      "O_{4}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4Asym",          "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4Asym_Mu",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O4Asym_El",       "A_{O4}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3",              "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3_Mu",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3_El",           "O3",                      "O_{3}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3Asym",          "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3Asym_Mu",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O3Asym_El",       "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2",              "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2_Mu",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2_El",           "O2",                      "O_{2}",             "Events", "", "", 100,  -5,   5   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2Asym",          "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2Asym_Mu",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "GenJ1Q1Chi2_O2Asym_El",       "A_{O2}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO2",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO3",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO4",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO7",            "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO2_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO3_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO4_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO7_El",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO2_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO3_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO4_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1_ChangeO7_Mu",         "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO2",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO3",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO4",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO7",        "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO2_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO3_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO4_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO7_El",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO2_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO3_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO4_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "EvtJ1Q1Chi2_ChangeO7_Mu",     "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.CreateTH1( fs );
    h1.Sumw2();

    // Create TH2D
    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.ClearTH2Info();
    h2.addNewTH2("TH2Chi2_dRlb_vs_Mlb",           "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2Chi2_dRlb_vs_Mlb_El",        "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2Chi2_dRlb_vs_Mlb_Mu",        "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2Chi2Matched_dRlb_vs_Mlb",    "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2Chi2Matched_dRlb_vs_Mlb_El", "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2Chi2Matched_dRlb_vs_Mlb_Mu", "",  "", "", "", "",  500,  0,   500, 1000,  0,  10 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Gen_vs_Evt_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Gen_vs_Evt_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_Par_vs_EvtPar_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_Par_vs_EvtPar_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenBB_vs_EvtBB_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenBB_vs_EvtBB_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7",     "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7_El",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.addNewTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7_Mu",  "",  "", "", "", "",  100,  -5,   5, 100,  -5,   5 );
    h2.CreateTH2( fs );
    h2.Sumw2();

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
    chain_  = new TChain(inputTTree_.c_str());

    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
    }

    EvtInfo.Register(chain_);
    GenInfo.Register(chain_);
    //JetInfo.Register(chain_,"PFJetInfo");
    topHadronicInfo.Register( chain_, "topHadronic" );
    IsoLepInfo.Register(      chain_, "isoLepton" );
    BJetInfo.Register(        chain_, "bJet"      );
    BbarJetInfo.Register(     chain_, "bbarJet"   );
    nonBJetColInfo.Register(  chain_, "nonBJetCol");

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
        GenO2 = Obs2( lepton.P3, j1.P3, b_quark.P3, bbar_quark.P3, charge );
        GenO3 = Obs3( lepton.P4, j1.P4, b_quark.P4, bbar_quark.P4, charge );
        GenO4 = Obs4( lepton.P3, j1.P3, b_quark.P3, bbar_quark.P3, charge );
        GenO7 = Obs7( az, b_quark.P3, bbar_quark.P3 );

        TLorentzVector P4_b, P4_bbar;
        P4_b.SetPtEtaPhiE( BJetInfo.Pt, BJetInfo.Eta, BJetInfo.Phi, BJetInfo.Energy );
        P4_bbar.SetPtEtaPhiE( BJetInfo.Pt, BJetInfo.Eta, BJetInfo.Phi, BJetInfo.Energy );

        // * Reconstruct Mlb
        float dRlb=1000000000.;
        int match_BFlavor=0;
        TLorentzVector P4_lb, isoL;
        isoL.SetPtEtaPhiE( IsoLepInfo.Pt, IsoLepInfo.Eta, IsoLepInfo.Phi, IsoLepInfo.Energy );
        if( topHadronicInfo.isAntiB == 1 )
        { 
            if( BJetInfo.GenFlavor == 5 ) match_BFlavor = BJetInfo.GenFlavor; 
            P4_lb = P4_b + isoL;
            dRlb  = P4_b.DeltaR(isoL);
        }
        else if( topHadronicInfo.isAntiB == 0 )
        {
            if( BbarJetInfo.GenFlavor == -5 ) match_BFlavor = BbarJetInfo.GenFlavor; 
            P4_lb = P4_bbar + isoL;
            dRlb  = P4_bbar.DeltaR(isoL);
        } 

        // * Parton ACP
        int ijet1=-1;
        float ipt=-1;
        for( int i=0; i<nonBJetColInfo.Size; i++ )
        {
            if( nonBJetColInfo.Pt[i] > ipt ) 
            {
                ijet1=i;
                ipt = nonBJetColInfo.Pt[i];
            }
        }
        bool isGoodJ1Q1=false;
        bool isGoodBB=false;
        bool isPerfectBB=false;
        double ParO2(0), ParO3(0), ParO4(0), ParO7(0);
        TLorentzVector  bPar, bbarPar, jet1;
        jet1.SetPtEtaPhiE( nonBJetColInfo.Pt[ijet1], nonBJetColInfo.Eta[ijet1], nonBJetColInfo.Phi[ijet1], nonBJetColInfo.Energy[ijet1] );
        bPar.SetPtEtaPhiE(    BJetInfo.GenPt,    BJetInfo.GenEta,     BJetInfo.GenPhi,    BJetInfo.Energy    );
        bbarPar.SetPtEtaPhiE( BbarJetInfo.GenPt, BbarJetInfo.GenEta,  BbarJetInfo.GenPhi, BbarJetInfo.Energy );
        ParO2 = Obs2( isoL.Vect(), jet1.Vect(), bPar.Vect(), bbarPar.Vect(), IsoLepInfo.Charge );
        ParO3 = Obs3( isoL, jet1, bPar, bbarPar, IsoLepInfo.Charge ); 
        ParO4 = Obs4( isoL.Vect(), jet1.Vect(), bPar.Vect(), bbarPar.Vect(), IsoLepInfo.Charge );
        ParO7 = Obs7( az, bPar.Vect(), bbarPar.Vect() );
        h1.GetTH1("Par_PID_b"   )->Fill( BJetInfo.GenFlavor    );
        h1.GetTH1("Par_PID_bbar")->Fill( BbarJetInfo.GenFlavor );

        if( jet1.DeltaR(j1.P4) < 0.5 ) isGoodJ1Q1=true;
        if( BJetInfo.GenFlavor == 5 && BbarJetInfo.GenFlavor == -5 ) isGoodBB=true;
        if( isGoodBB ) isPerfectBB=true; 
        //if( isGoodBB && 
        //        b_quark.P4.DeltaR(P4_b) < 0.5 && 
        //        bbar_quark.P4.DeltaR(P4_bbar) < 0.5 ) isPerfectBB=true;

        if( ParO3 == 0 && isGoodBB ) 
        {
            std::cout<<">> [WARING] Parton O3 == "<<ParO3<<", REOC: "<<O3_<<std::endl;
            isoL.Print();
            jet1.Print();
            bPar.Print();
            bbarPar.Print();
            std::cout<<">>          BJetInfo.GenPt "<<BJetInfo.GenPt<<", BbarJetInfo.GenPt "<<BbarJetInfo.GenPt<<std::endl;
            std::cout<<">>          BJetInfo.GenEta "<<BJetInfo.GenEta<<", BbarJetInfo.GenEta "<<BbarJetInfo.GenEta<<std::endl;
            std::cout<<">>          BJetInfo.GenPhi "<<BJetInfo.GenPhi<<", BbarJetInfo.GenPhi "<<BbarJetInfo.GenPhi<<std::endl;
            std::cout<<">>          BJetInfo.Energy "<<BJetInfo.Energy<<", BbarJetInfo.Energy "<<BbarJetInfo.Energy<<std::endl;
            std::cout<<">>          BJetInfo.GenFlavor "<<BJetInfo.GenFlavor<<", BbarJetInfo.GenFlavor "<<BbarJetInfo.GenFlavor<<std::endl;
        }

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
        h2.GetTH2("TH2_Gen_vs_Evt_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
        h2.GetTH2("TH2_Gen_vs_Evt_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
        h2.GetTH2("TH2_Gen_vs_Evt_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
        h2.GetTH2("TH2_Gen_vs_Evt_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
        if( isGoodBB )
        {
            if( isPerfectBB )
            {
                checkObsChange( h1.GetTH1("EvtBB_ChangeO2"), O2_, GenO2, wrtevt );
                checkObsChange( h1.GetTH1("EvtBB_ChangeO3"), O3_, GenO3, wrtevt );
                checkObsChange( h1.GetTH1("EvtBB_ChangeO4"), O4_, GenO4, wrtevt );
                checkObsChange( h1.GetTH1("EvtBB_ChangeO7"), O7_, GenO7, wrtevt );
                h1.GetTH1("EvtBB_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("EvtBB_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("EvtBB_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("EvtBB_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                h1.GetTH1("GenBB_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("GenBB_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("GenBB_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("GenBB_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtBB_O2Asym"), O2_, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtBB_O3Asym"), O3_, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtBB_O4Asym"), O4_, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtBB_O7Asym"), O7_, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("GenBB_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("GenBB_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("GenBB_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("GenBB_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                h2.GetTH2("TH2_GenBB_vs_EvtBB_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                h2.GetTH2("TH2_GenBB_vs_EvtBB_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                h2.GetTH2("TH2_GenBB_vs_EvtBB_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                h2.GetTH2("TH2_GenBB_vs_EvtBB_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                if( isGoodJ1Q1 )
                {
                    checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO2"), O2_, GenO2, wrtevt );
                    checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO3"), O3_, GenO3, wrtevt );
                    checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO4"), O4_, GenO4, wrtevt );
                    checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO7"), O7_, GenO7, wrtevt );
                    h1.GetTH1("EvtJ1Q1_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("EvtJ1Q1_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("EvtJ1Q1_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("EvtJ1Q1_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                    h1.GetTH1("GenJ1Q1_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("GenJ1Q1_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("GenJ1Q1_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("GenJ1Q1_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("EvtJ1Q1_O2Asym"), O2_, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("EvtJ1Q1_O3Asym"), O3_, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("EvtJ1Q1_O4Asym"), O4_, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("EvtJ1Q1_O7Asym"), O7_, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("GenJ1Q1_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("GenJ1Q1_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("GenJ1Q1_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("GenJ1Q1_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                    h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                }
            }
            if( fabs(ParO3) >= 0 )
            {
                checkObsChange( h1.GetTH1("Par_ChangeO2"), O2_, ParO2, wrtevt );
                checkObsChange( h1.GetTH1("Par_ChangeO3"), O3_, ParO3, wrtevt );
                checkObsChange( h1.GetTH1("Par_ChangeO4"), O4_, ParO4, wrtevt );
                checkObsChange( h1.GetTH1("Par_ChangeO7"), O7_, ParO7, wrtevt );
                h1.GetTH1("EvtPar_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("EvtPar_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("EvtPar_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("EvtPar_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                h1.GetTH1("Par_O2")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                h1.GetTH1("Par_O3")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                h1.GetTH1("Par_O4")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                h1.GetTH1("Par_O7")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("EvtPar_O2Asym"), O2_, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("EvtPar_O3Asym"), O3_, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("EvtPar_O4Asym"), O4_, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("EvtPar_O7Asym"), O7_, wrtevt*acpWrtO7 );
                fillAsym( h1.GetTH1("Par_O2Asym"), ParO2, wrtevt*acpWrtO2 );
                fillAsym( h1.GetTH1("Par_O3Asym"), ParO3, wrtevt*acpWrtO3 );
                fillAsym( h1.GetTH1("Par_O4Asym"), ParO4, wrtevt*acpWrtO4 );
                fillAsym( h1.GetTH1("Par_O7Asym"), ParO7, wrtevt*acpWrtO7 );
                h2.GetTH2("TH2_Par_vs_EvtPar_O2")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                h2.GetTH2("TH2_Par_vs_EvtPar_O3")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                h2.GetTH2("TH2_Par_vs_EvtPar_O4")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                h2.GetTH2("TH2_Par_vs_EvtPar_O7")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                //if( !test  ) 
                //{
                //    std::cout<<">> [WARING] Parton O3 == "<<Obs3( isoL, jet1, bPar, bbarPar, IsoLepInfo.Charge, true )<<", REOC: "<<O3_<<std::endl;
                //    isoL.Print();
                //    jet1.Print();
                //    std::cout<<">>          BJet: "<<endl;
                //    bPar.Print();
                //    std::cout<<">>          BbarJet: "<<endl;
                //    bbarPar.Print();
                //    std::cout<<">>          BJetInfo.GenPt "<<BJetInfo.GenPt<<", BbarJetInfo.GenPt "<<BbarJetInfo.GenPt<<std::endl;
                //    std::cout<<">>          BJetInfo.GenEta "<<BJetInfo.GenEta<<", BbarJetInfo.GenEta "<<BbarJetInfo.GenEta<<std::endl;
                //    std::cout<<">>          BJetInfo.GenPhi "<<BJetInfo.GenPhi<<", BbarJetInfo.GenPhi "<<BbarJetInfo.GenPhi<<std::endl;
                //    std::cout<<">>          BJetInfo.Energy "<<BJetInfo.Energy<<", BbarJetInfo.Energy "<<BbarJetInfo.Energy<<std::endl;
                //    std::cout<<">>          BJetInfo.GenFlavor "<<BJetInfo.GenFlavor<<", BbarJetInfo.GenFlavor "<<BbarJetInfo.GenFlavor<<std::endl;
                //}
            }
        }
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
            h2.GetTH2("TH2_Gen_vs_Evt_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
            h2.GetTH2("TH2_Gen_vs_Evt_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
            h2.GetTH2("TH2_Gen_vs_Evt_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
            h2.GetTH2("TH2_Gen_vs_Evt_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
            if( isGoodBB )
            {
                if( isPerfectBB )
                {
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO2_Mu"), O2_, GenO2, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO3_Mu"), O3_, GenO3, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO4_Mu"), O4_, GenO4, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO7_Mu"), O7_, GenO7, wrtevt );
                    h1.GetTH1("EvtBB_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("EvtBB_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("EvtBB_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("EvtBB_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                    h1.GetTH1("GenBB_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("GenBB_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("GenBB_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("GenBB_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("EvtBB_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("EvtBB_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("EvtBB_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("EvtBB_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("GenBB_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("GenBB_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("GenBB_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("GenBB_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                    if( isGoodJ1Q1 )
                    {
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO2_Mu"), O2_, GenO2, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO3_Mu"), O3_, GenO3, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO4_Mu"), O4_, GenO4, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO7_Mu"), O7_, GenO7, wrtevt );
                        h1.GetTH1("EvtJ1Q1_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                        h1.GetTH1("EvtJ1Q1_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                        h1.GetTH1("EvtJ1Q1_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                        h1.GetTH1("EvtJ1Q1_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                        h1.GetTH1("GenJ1Q1_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                        h1.GetTH1("GenJ1Q1_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                        h1.GetTH1("GenJ1Q1_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                        h1.GetTH1("GenJ1Q1_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                    }
                }
                if( fabs(ParO3) >= 0 )
                {
                    checkObsChange( h1.GetTH1("Par_ChangeO2_Mu"), O2_, ParO2, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO3_Mu"), O3_, ParO3, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO4_Mu"), O4_, ParO4, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO7_Mu"), O7_, ParO7, wrtevt );
                    h1.GetTH1("EvtPar_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("EvtPar_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("EvtPar_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("EvtPar_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                    h1.GetTH1("Par_O2_Mu")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("Par_O3_Mu")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("Par_O4_Mu")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("Par_O7_Mu")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("EvtPar_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("EvtPar_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("EvtPar_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("EvtPar_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("Par_O2Asym_Mu"), ParO2, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("Par_O3Asym_Mu"), ParO3, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("Par_O4Asym_Mu"), ParO4, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("Par_O7Asym_Mu"), ParO7, wrtevt*acpWrtO7 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O2_Mu")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O3_Mu")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O4_Mu")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O7_Mu")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                }
            }

            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Fill( P4_lb.M(), wrtevt );
                h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_Mu")->Fill( P4_lb.M(), wrtevt );
                h2.GetTH2("TH2Chi2_dRlb_vs_Mlb"   )->Fill( P4_lb.M(), dRlb, wrtevt );
                h2.GetTH2("TH2Chi2_dRlb_vs_Mlb_Mu")->Fill( P4_lb.M(), dRlb, wrtevt );
                if( match_BFlavor != 0 )
                {
                    h2.GetTH2("TH2Chi2Matched_dRlb_vs_Mlb"   )->Fill( P4_lb.M(), dRlb, wrtevt );
                    h2.GetTH2("TH2Chi2Matched_dRlb_vs_Mlb_Mu")->Fill( P4_lb.M(), dRlb, wrtevt );
                    h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl"   )->Fill( P4_lb.M(), wrtevt );
                    h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl_Mu")->Fill( P4_lb.M(), wrtevt );
                }
                if( maxMlbCut_ > P4_lb.M() && minMlbCut_ <= P4_lb.M() )
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
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                    if( isGoodBB )
                    {
                        if( isPerfectBB )
                        {
                            h1.GetTH1("EvtBBChi2_dR_j1q1"  )->Fill(jet1.DeltaR(j1.P4));
                            h1.GetTH1("EvtBBChi2_Flavor_j1")->Fill(nonBJetColInfo.GenFlavor[ijet1]);
                            if( nonBJetColInfo.GenFlavor[ijet1] == j1.PdgID ) h1.GetTH1("EvtBBChi2_MatchedPID_j1q1")->Fill(1);
                            else h1.GetTH1("EvtBBChi2_MatchedPID_j1q1")->Fill(0);

                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO2"), O2_, GenO2, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO3"), O3_, GenO3, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO4"), O4_, GenO4, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO7"), O7_, GenO7, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO2_Mu"), O2_, GenO2, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO3_Mu"), O3_, GenO3, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO4_Mu"), O4_, GenO4, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO7_Mu"), O7_, GenO7, wrtevt );
                            h1.GetTH1("EvtBBChi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtBBChi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtBBChi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtBBChi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("EvtBBChi2_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtBBChi2_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtBBChi2_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtBBChi2_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("GenBBChi2_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("GenBBChi2_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("GenBBChi2_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("GenBBChi2_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("GenBBChi2_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("GenBBChi2_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("GenBBChi2_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("GenBBChi2_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("GenBBChi2_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("GenBBChi2_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("GenBBChi2_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("GenBBChi2_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("GenBBChi2_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("GenBBChi2_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("GenBBChi2_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("GenBBChi2_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            if( isGoodJ1Q1 )
                            {
                                h1.GetTH1("EvtJ1Q1Chi2_dR_j1q1"  )->Fill(jet1.DeltaR(j1.P4));
                                h1.GetTH1("EvtJ1Q1Chi2_Flavor_j1")->Fill(nonBJetColInfo.GenFlavor[ijet1]);
                                if( nonBJetColInfo.GenFlavor[ijet1] == j1.PdgID ) h1.GetTH1("EvtJ1Q1Chi2_MatchedPID_j1q1")->Fill(1);
                                else h1.GetTH1("EvtJ1Q1Chi2_MatchedPID_j1q1")->Fill(0);

                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO2"), O2_, GenO2, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO3"), O3_, GenO3, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO4"), O4_, GenO4, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO7"), O7_, GenO7, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO2_Mu"), O2_, GenO2, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO3_Mu"), O3_, GenO3, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO4_Mu"), O4_, GenO4, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO7_Mu"), O7_, GenO7, wrtevt );
                                h1.GetTH1("EvtJ1Q1Chi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("EvtJ1Q1Chi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("EvtJ1Q1Chi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("EvtJ1Q1Chi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("EvtJ1Q1Chi2_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("EvtJ1Q1Chi2_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("EvtJ1Q1Chi2_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("EvtJ1Q1Chi2_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("GenJ1Q1Chi2_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("GenJ1Q1Chi2_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("GenJ1Q1Chi2_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("GenJ1Q1Chi2_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("GenJ1Q1Chi2_O2_Mu")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("GenJ1Q1Chi2_O3_Mu")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("GenJ1Q1Chi2_O4_Mu")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("GenJ1Q1Chi2_O7_Mu")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O2Asym_Mu"), GenO2, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O3Asym_Mu"), GenO3, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O4Asym_Mu"), GenO4, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O7Asym_Mu"), GenO7, wrtevt*acpWrtO7 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2_Mu")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3_Mu")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4_Mu")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7_Mu")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            }
                        }
                        if( fabs(ParO3) >= 0 )
                        {
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO2"), O2_, ParO2, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO3"), O3_, ParO3, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO4"), O4_, ParO4, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO7"), O7_, ParO7, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO2_Mu"), O2_, ParO2, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO3_Mu"), O3_, ParO3, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO4_Mu"), O4_, ParO4, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO7_Mu"), O7_, ParO7, wrtevt );
                            h1.GetTH1("EvtParChi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtParChi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtParChi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtParChi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("EvtParChi2_O2_Mu")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtParChi2_O3_Mu")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtParChi2_O4_Mu")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtParChi2_O7_Mu")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("ParChi2_O2")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("ParChi2_O3")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("ParChi2_O4")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("ParChi2_O7")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("ParChi2_O2_Mu")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("ParChi2_O3_Mu")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("ParChi2_O4_Mu")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("ParChi2_O7_Mu")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtParChi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtParChi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtParChi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtParChi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtParChi2_O2Asym_Mu"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtParChi2_O3Asym_Mu"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtParChi2_O4Asym_Mu"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtParChi2_O7Asym_Mu"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("ParChi2_O2Asym"), ParO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("ParChi2_O3Asym"), ParO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("ParChi2_O4Asym"), ParO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("ParChi2_O7Asym"), ParO7, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("ParChi2_O2Asym_Mu"), ParO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("ParChi2_O3Asym_Mu"), ParO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("ParChi2_O4Asym_Mu"), ParO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("ParChi2_O7Asym_Mu"), ParO7, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O2")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O3")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O4")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O7")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O2_Mu")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O3_Mu")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O4_Mu")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O7_Mu")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                        }
                    }
                }
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
            h2.GetTH2("TH2_Gen_vs_Evt_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
            h2.GetTH2("TH2_Gen_vs_Evt_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
            h2.GetTH2("TH2_Gen_vs_Evt_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
            h2.GetTH2("TH2_Gen_vs_Evt_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
            if( isGoodBB )
            {
                if( isPerfectBB )
                {
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO2_El"), O2_, GenO2, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO3_El"), O3_, GenO3, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO4_El"), O4_, GenO4, wrtevt );
                    checkObsChange( h1.GetTH1("EvtBB_ChangeO7_El"), O7_, GenO7, wrtevt );
                    h1.GetTH1("EvtBB_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("EvtBB_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("EvtBB_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("EvtBB_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                    h1.GetTH1("GenBB_O2_El")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("GenBB_O3_El")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("GenBB_O4_El")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("GenBB_O7_El")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("EvtBB_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("EvtBB_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("EvtBB_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("EvtBB_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("GenBB_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("GenBB_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("GenBB_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("GenBB_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                    h2.GetTH2("TH2_GenBB_vs_EvtBB_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                    if( isGoodJ1Q1 )
                    {
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO2_El"), O2_, GenO2, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO3_El"), O3_, GenO3, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO4_El"), O4_, GenO4, wrtevt );
                        checkObsChange( h1.GetTH1("EvtJ1Q1_ChangeO7_El"), O7_, GenO7, wrtevt );
                        h1.GetTH1("EvtJ1Q1_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                        h1.GetTH1("EvtJ1Q1_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                        h1.GetTH1("EvtJ1Q1_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                        h1.GetTH1("EvtJ1Q1_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                        h1.GetTH1("GenJ1Q1_O2_El")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                        h1.GetTH1("GenJ1Q1_O3_El")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                        h1.GetTH1("GenJ1Q1_O4_El")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                        h1.GetTH1("GenJ1Q1_O7_El")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                        fillAsym( h1.GetTH1("EvtJ1Q1_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
                        fillAsym( h1.GetTH1("GenJ1Q1_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                        h2.GetTH2("TH2_GenJ1Q1_vs_EvtJ1Q1_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                    }
                }
                if( fabs(ParO3) >= 0 )
                {
                    checkObsChange( h1.GetTH1("Par_ChangeO2_El"), O2_, ParO2, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO3_El"), O3_, ParO3, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO4_El"), O4_, ParO4, wrtevt );
                    checkObsChange( h1.GetTH1("Par_ChangeO7_El"), O7_, ParO7, wrtevt );
                    h1.GetTH1("EvtPar_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("EvtPar_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("EvtPar_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("EvtPar_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                    h1.GetTH1("Par_O2_El")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                    h1.GetTH1("Par_O3_El")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                    h1.GetTH1("Par_O4_El")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                    h1.GetTH1("Par_O7_El")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("EvtPar_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("EvtPar_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("EvtPar_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("EvtPar_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                    fillAsym( h1.GetTH1("Par_O2Asym_El"), ParO2, wrtevt*acpWrtO2 );
                    fillAsym( h1.GetTH1("Par_O3Asym_El"), ParO3, wrtevt*acpWrtO3 );
                    fillAsym( h1.GetTH1("Par_O4Asym_El"), ParO4, wrtevt*acpWrtO4 );
                    fillAsym( h1.GetTH1("Par_O7Asym_El"), ParO7, wrtevt*acpWrtO7 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O2_El")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O3_El")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O4_El")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                    h2.GetTH2("TH2_Par_vs_EvtPar_O7_El")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                }
            }
            if( maxChi2Cut_ > minChi2_ && minChi2Cut_ <= minChi2_ )
            {
                h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Fill( P4_lb.M(), wrtevt );
                h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_El")->Fill( P4_lb.M(), wrtevt );
                h2.GetTH2("TH2Chi2_dRlb_vs_Mlb"   )->Fill( P4_lb.M(), dRlb, wrtevt );
                h2.GetTH2("TH2Chi2_dRlb_vs_Mlb_El")->Fill( P4_lb.M(), dRlb, wrtevt );
                if( match_BFlavor != 0 )
                {
                    h2.GetTH2("TH2Chi2Matched_dRlb_vs_Mlb"   )->Fill( P4_lb.M(), dRlb, wrtevt );
                    h2.GetTH2("TH2Chi2Matched_dRlb_vs_Mlb_El")->Fill( P4_lb.M(), dRlb, wrtevt );
                    h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl"   )->Fill( P4_lb.M(), wrtevt );
                    h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl_El")->Fill( P4_lb.M(), wrtevt );
                }
                if( maxMlbCut_ > P4_lb.M() && minMlbCut_ <= P4_lb.M() )
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
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4);
                    h2.GetTH2("TH2Chi2_Gen_vs_Evt_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7);
                    if( isGoodBB )
                    {
                        if( isPerfectBB )
                        {
                            h1.GetTH1("EvtBBChi2_dR_j1q1"  )->Fill(jet1.DeltaR(j1.P4));
                            h1.GetTH1("EvtBBChi2_Flavor_j1")->Fill(nonBJetColInfo.GenFlavor[ijet1]);
                            if( nonBJetColInfo.GenFlavor[ijet1] == j1.PdgID ) h1.GetTH1("EvtBBChi2_MatchedPID_j1q1")->Fill(1);
                            else h1.GetTH1("EvtBBChi2_MatchedPID_j1q1")->Fill(0);

                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO2"), O2_, GenO2, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO3"), O3_, GenO3, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO4"), O4_, GenO4, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO7"), O7_, GenO7, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO2_El"), O2_, GenO2, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO3_El"), O3_, GenO3, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO4_El"), O4_, GenO4, wrtevt );
                            checkObsChange( h1.GetTH1("EvtBBChi2_ChangeO7_El"), O7_, GenO7, wrtevt );
                            h1.GetTH1("EvtBBChi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtBBChi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtBBChi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtBBChi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("EvtBBChi2_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtBBChi2_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtBBChi2_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtBBChi2_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("GenBBChi2_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("GenBBChi2_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("GenBBChi2_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("GenBBChi2_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("GenBBChi2_O2_El")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("GenBBChi2_O3_El")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("GenBBChi2_O4_El")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("GenBBChi2_O7_El")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtBBChi2_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("GenBBChi2_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("GenBBChi2_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("GenBBChi2_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("GenBBChi2_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("GenBBChi2_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("GenBBChi2_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("GenBBChi2_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("GenBBChi2_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_GenBB_vs_EvtBB_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            if( isGoodJ1Q1 )
                            {
                                h1.GetTH1("EvtJ1Q1Chi2_dR_j1q1"  )->Fill(jet1.DeltaR(j1.P4));
                                h1.GetTH1("EvtJ1Q1Chi2_Flavor_j1")->Fill(nonBJetColInfo.GenFlavor[ijet1]);
                                if( nonBJetColInfo.GenFlavor[ijet1] == j1.PdgID ) h1.GetTH1("EvtJ1Q1Chi2_MatchedPID_j1q1")->Fill(1);
                                else h1.GetTH1("EvtJ1Q1Chi2_MatchedPID_j1q1")->Fill(0);

                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO2"), O2_, GenO2, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO3"), O3_, GenO3, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO4"), O4_, GenO4, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO7"), O7_, GenO7, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO2_El"), O2_, GenO2, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO3_El"), O3_, GenO3, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO4_El"), O4_, GenO4, wrtevt );
                                checkObsChange( h1.GetTH1("EvtJ1Q1Chi2_ChangeO7_El"), O7_, GenO7, wrtevt );
                                h1.GetTH1("EvtJ1Q1Chi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("EvtJ1Q1Chi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("EvtJ1Q1Chi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("EvtJ1Q1Chi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("EvtJ1Q1Chi2_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("EvtJ1Q1Chi2_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("EvtJ1Q1Chi2_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("EvtJ1Q1Chi2_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("GenJ1Q1Chi2_O2")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("GenJ1Q1Chi2_O3")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("GenJ1Q1Chi2_O4")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("GenJ1Q1Chi2_O7")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                                h1.GetTH1("GenJ1Q1Chi2_O2_El")->Fill( GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h1.GetTH1("GenJ1Q1Chi2_O3_El")->Fill( GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h1.GetTH1("GenJ1Q1Chi2_O4_El")->Fill( GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h1.GetTH1("GenJ1Q1Chi2_O7_El")->Fill( GenO7/wrtobs, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("EvtJ1Q1Chi2_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O2Asym"), GenO2, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O3Asym"), GenO3, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O4Asym"), GenO4, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O7Asym"), GenO7, wrtevt*acpWrtO7 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O2Asym_El"), GenO2, wrtevt*acpWrtO2 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O3Asym_El"), GenO3, wrtevt*acpWrtO3 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O4Asym_El"), GenO4, wrtevt*acpWrtO4 );
                                fillAsym( h1.GetTH1("GenJ1Q1Chi2_O7Asym_El"), GenO7, wrtevt*acpWrtO7 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O2_El")->Fill(O2_/wrtobs,GenO2/wrtobs, wrtevt*acpWrtO2 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O3_El")->Fill(O3_/wrtobs,GenO3/wrtobs, wrtevt*acpWrtO3 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O4_El")->Fill(O4_/wrtobs,GenO4/wrtobs, wrtevt*acpWrtO4 );
                                h2.GetTH2("TH2Chi2_GenJ1Q1_vs_EvtJ1Q1_O7_El")->Fill(O7_/wrtobs,GenO7/wrtobs, wrtevt*acpWrtO7 );
                            }
                        }
                        if( fabs(ParO3) >= 0 )
                        {
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO2"), O2_, ParO2, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO3"), O3_, ParO3, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO4"), O4_, ParO4, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO7"), O7_, ParO7, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO2_El"), O2_, ParO2, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO3_El"), O3_, ParO3, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO4_El"), O4_, ParO4, wrtevt );
                            checkObsChange( h1.GetTH1("ParChi2_ChangeO7_El"), O7_, ParO7, wrtevt );
                            h1.GetTH1("EvtParChi2_O2")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtParChi2_O3")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtParChi2_O4")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtParChi2_O7")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("EvtParChi2_O2_El")->Fill( O2_/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("EvtParChi2_O3_El")->Fill( O3_/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("EvtParChi2_O4_El")->Fill( O4_/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("EvtParChi2_O7_El")->Fill( O7_/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("ParChi2_O2")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("ParChi2_O3")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("ParChi2_O4")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("ParChi2_O7")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                            h1.GetTH1("ParChi2_O2_El")->Fill( ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h1.GetTH1("ParChi2_O3_El")->Fill( ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h1.GetTH1("ParChi2_O4_El")->Fill( ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h1.GetTH1("ParChi2_O7_El")->Fill( ParO7/wrtobs, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtParChi2_O2Asym"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtParChi2_O3Asym"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtParChi2_O4Asym"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtParChi2_O7Asym"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("EvtParChi2_O2Asym_El"), O2_, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("EvtParChi2_O3Asym_El"), O3_, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("EvtParChi2_O4Asym_El"), O4_, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("EvtParChi2_O7Asym_El"), O7_, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("ParChi2_O2Asym"), ParO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("ParChi2_O3Asym"), ParO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("ParChi2_O4Asym"), ParO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("ParChi2_O7Asym"), ParO7, wrtevt*acpWrtO7 );
                            fillAsym( h1.GetTH1("ParChi2_O2Asym_El"), ParO2, wrtevt*acpWrtO2 );
                            fillAsym( h1.GetTH1("ParChi2_O3Asym_El"), ParO3, wrtevt*acpWrtO3 );
                            fillAsym( h1.GetTH1("ParChi2_O4Asym_El"), ParO4, wrtevt*acpWrtO4 );
                            fillAsym( h1.GetTH1("ParChi2_O7Asym_El"), ParO7, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O2")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O3")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O4")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O7")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O2_El")->Fill(O2_/wrtobs,ParO2/wrtobs, wrtevt*acpWrtO2 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O3_El")->Fill(O3_/wrtobs,ParO3/wrtobs, wrtevt*acpWrtO3 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O4_El")->Fill(O4_/wrtobs,ParO4/wrtobs, wrtevt*acpWrtO4 );
                            h2.GetTH2("TH2Chi2_Par_vs_EvtPar_O7_El")->Fill(O7_/wrtobs,ParO7/wrtobs, wrtevt*acpWrtO7 );
                        }
                    }
                }
            }
        }

    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void SemiLeptanicResultsCheckObs::endJob(){
    TH1D* h_tmp1    = (TH1D*)h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl"   )->Clone(); h_tmp1   ->Rebin(10);
    TH1D* h_tmp1_el = (TH1D*)h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl_El")->Clone(); h_tmp1_el->Rebin(10);
    TH1D* h_tmp1_mu = (TH1D*)h1.GetTH1("EvtChi2Matched_Top_Leptonic_Mbl_Mu")->Clone(); h_tmp1_mu->Rebin(10);
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl"   )->Add(h_tmp1   ); 
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl_El")->Add(h_tmp1_el);
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl_Mu")->Add(h_tmp1_mu);
    TH1D* h_tmp2    = (TH1D*)h1.GetTH1("EvtChi2_Top_Leptonic_Mbl"   )->Clone(); h_tmp2   ->Rebin(10);
    TH1D* h_tmp2_el = (TH1D*)h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_El")->Clone(); h_tmp2_el->Rebin(10);
    TH1D* h_tmp2_mu = (TH1D*)h1.GetTH1("EvtChi2_Top_Leptonic_Mbl_Mu")->Clone(); h_tmp2_mu->Rebin(10);
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl"   )->Divide(h_tmp2   ); 
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl_El")->Divide(h_tmp2_el);
    h1.GetTH1("EffMatched_Top_Leptonic_Mbl_Mu")->Divide(h_tmp2_mu);

    delete h_tmp2; 
    delete h_tmp2_el;
    delete h_tmp2_mu;
    delete h_tmp1;
    delete h_tmp1_el;
    delete h_tmp1_mu;

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
