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

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/GenParticle.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicResultsCheck.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorElectron.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorMuon.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorJet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorVertex.h" 

//
// constructors and destructor
//
SemiLeptanicResultsCheck::SemiLeptanicResultsCheck(const edm::ParameterSet& iConfig) : 
    maxEvents_(             iConfig.getParameter<int>("MaxEvents")), 
    reportEvery_(           iConfig.getParameter<int>("ReportEvery")),
    inputTTree_(            iConfig.getParameter<std::string>("InputTTree")),
    inputFiles_(            iConfig.getParameter<std::vector<std::string> >("InputFiles")),
    HLT_MuChannel_(         iConfig.getParameter<std::vector<int>>("HLT_MuChannel")),
    HLT_ElChannel_(         iConfig.getParameter<std::vector<int>>("HLT_ElChannel")),
    selPars_Vertex_(        iConfig.getParameter<edm::ParameterSet>("SelPars_Vertex")),
    selPars_Jet_(           iConfig.getParameter<edm::ParameterSet>("SelPars_Jet")),
    selPars_BJet_(          iConfig.getParameter<edm::ParameterSet>("SelPars_BJet")),
    selPars_NonBJet_(       iConfig.getParameter<edm::ParameterSet>("SelPars_NonBJet")),
    selPars_LooseLepton_(   iConfig.getParameter<edm::ParameterSet>("SelPars_LooseLepton")),
    selPars_TightMuon_(     iConfig.getParameter<edm::ParameterSet>("SelPars_TightMuon")),
    selPars_TightElectron_( iConfig.getParameter<edm::ParameterSet>("SelPars_TightElectron")),
    dR_IsoLeptonFromJets_(  iConfig.getParameter<double>("dR_IsoLeptonFromJets")),
    Owrt_(                  iConfig.getParameter<double>("Owrt")),
    //isSignal_(              iConfig.getParameter<bool>("isSingal"))
    Debug_(                 iConfig.getParameter<bool>("Debug"))
{}

SemiLeptanicResultsCheck::~SemiLeptanicResultsCheck()
{ 
    delete chain_;
}

// ------------ Other function -------------
std::string SemiLeptanicResultsCheck::int2str( int i )
{
    std::string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}
template<class TH1>
void SemiLeptanicResultsCheck::setObservableHist( TH1* h, string ob ){
    h->GetXaxis()->SetBinLabel(1,(ob+"<0").c_str());
    h->GetXaxis()->SetBinLabel(2,(ob+">0").c_str());
}
template<class TH1>
void SemiLeptanicResultsCheck::fillObservableHist( TH1* h, double ob, string obs, double wrt ){
    if( ob > 0 ){
        h->Fill((obs+">0").c_str(), wrt);
    }else{
        h->Fill((obs+"<0").c_str(), wrt);
    }
}
template <class Object>
bool SemiLeptanicResultsCheck::getHighPtSelectMo( vector<Object> col, Object &obj, int mo )
{
    bool matched=false;
    int o1(-1);
    double pt(0);
    const int size=col.size();
    for( int i=0; i<size; i++)
    {
        if( mo == 0 || abs(col[i].Mo1PdgID) == mo || abs(col[i].Mo2PdgID) == mo ){
            if( pt < col[i].Pt )
            {
                pt=col[i].Pt;
                o1=i;
                matched=true;
            }
        }
    }
    obj=col[o1];
    return matched;
}
    template <class Object>
void SemiLeptanicResultsCheck::getHighPtObject( vector<Object> col, Object &obj )
{
    int o1(-1);
    double pt(0);
    const int size=col.size();
    for( int i=0; i<size; i++)
    {
        if( pt < col[i].Pt )
        {
            pt=col[i].Pt;
            o1=i;
        }
    }
    obj=col[o1];
}
template <class Object>
void SemiLeptanicResultsCheck::get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 )
{
    int o1, o2;
    double pt1, pt2;
    o1=o2=-1;
    pt1=pt2=0;
    const int size=col.size();
    for( int i=0; i<size; i++){
        if( pt1 < col[i].Pt ){
            pt2=pt1;
            pt1=col[i].Pt;
            o2=o1;
            o1=i;
        }else if( pt2 < col[i].Pt ){
            pt2=col[i].Pt;
            o2=i;
        }
    }
    obj1=col[o1];
    obj2=col[o2];
}

bool SemiLeptanicResultsCheck::isIsoLeptonFromJets( Lepton lepton, vector<Jet> jetCol, double dR )
{
    bool isIsoLepFromJets = true;
    int const jetcolSize = jetCol.size();
    for( int i=0; i<jetcolSize; i++ )
    {
        double dr = lepton.P4.DeltaR( jetCol[i].P4 );
        if( dr < dR )
        { 
            isIsoLepFromJets = false;
            break;
        }
    }
    return isIsoLepFromJets;
}
double SemiLeptanicResultsCheck::Obs2( TVector3 isoLep, TVector3 hardJet, TVector3 bjet1, TVector3 bjet2 )
{
    TVector3 O2_1v =  bjet1 + bjet2;
    TVector3 O2_2v = isoLep.Cross( hardJet );
    double O2 = O2_1v.Dot( O2_2v );
    return O2;
}
double SemiLeptanicResultsCheck::Obs7( TVector3 beam, TVector3 bjet1, TVector3 bjet2 )
{
    double O7_1z = beam.Dot( bjet1 - bjet2 );
    double O7_2z = beam.Dot( bjet1.Cross( bjet2 ));
    double O7 = O7_1z * O7_2z;
    return O7;
}
// ------------ method called once each job just before starting event loop  ------------
void SemiLeptanicResultsCheck::beginJob()
{
    // Create TH1D
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.ClearTH1Info();
    h1.addNewTH1( "Evt_O7",            "O7",                        "O_{7}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O7_Mu",         "O7",                        "O_{7}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O7_El",         "O7",                        "O_{7}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O7Asym",        "A_{O7}",                    "",                   "Events", "",    "",  2,   0,   2   );
    h1.addNewTH1( "Evt_O7Asym_Mu",     "A_{O7}",                    "",                   "Events", "",    "",  2,   0,   2   );
    h1.addNewTH1( "Evt_O7Asym_El",     "A_{O7}",                    "",                   "Events", "",    "",  2,   0,   2   );
    h1.addNewTH1( "Evt_O2",            "O2",                        "O_{2}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O2_Mu",         "O2",                        "O_{2}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O2_El",         "O2",                        "O_{2}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O2Asym",        "A_{O2}",                    "",                   "Events", "",    "",  2,   0,   2   );
    h1.addNewTH1( "Evt_O2Asym_Mu",     "A_{O2}",                    "",                   "Events", "",    "",  2,   0,   2   );
    h1.addNewTH1( "Evt_O2Asym_El",     "A_{O2}",                    "",                   "Events", "",    "",  2,   0,   2   );

    h1.addNewTH1( "Evt_isoMu_Pt",      "pT of isoMuon",             "p_{T}(#mu)",         "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoMu_Et",      "ET of isoMuon",             "E_{T}(#mu)",         "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoMu_E",       "Energy of isoMuon",         "Energy(#mu)",        "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoMu_Eta",     "Eta of isoMuon",            "#eta(#mu)",          "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_isoMu_Phi",     "Phi of isoMuon",            "#phi(#mu)",          "Yields", "",    "", 64,  -3.2, 3.2 );

    h1.addNewTH1( "Evt_isoEl_Pt",      "pT of isEle",               "p_{T}(e)",           "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoEl_Et",      "ET of isEle",               "E_{T}(e)",           "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoEl_E",       "Energy of isEle",           "Energy(e)",          "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_isoEl_Eta",     "Eta of isEle",              "#eta(e)",            "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_isoEl_Phi",     "Phi of isEle",              "#phi(e)",            "Yields", "",    "", 64,  -3.2, 3.2 );

    h1.addNewTH1( "Evt_HardJet_Pt",    "pT of HardJet",             "p_{T}(HardJet)",     "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_M",     "Mass of HardJet",           "Mass(HardJet)",      "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_E",     "Energy of HardJet",         "Energy(HardJet)",    "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_HardJet_Eta",   "Eta of HardJet",            "#eta(HardJet)",      "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_HardJet_Phi",   "Phi of HardJet",            "#phi(HardJet)",      "Yields", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_HardJet_BTag",  "HardJet b-tagged",          "bTag",               "Yields", "",    "", 100,  0,   1   );
    
    h1.addNewTH1( "Evt_Events",        "",                          "",                   "Events", "",    "",  1,   1,   2   );
    h1.addNewTH1( "Evt_NJets",         "Num. of jets",              "N(j)",               "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NSelJets",      "Num. of selected jets",     "N(selected j)",      "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NbJets",        "Num. of b-jets",            "N(B-tagged j)",      "Events", "",    "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NLeptons",      "Num. of leptons",           "N(lep)",             "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NSelLeptons",   "Num. of selected leptons",  "N(selected lep)",    "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NMuons",        "Num. of muon",              "N(#mu)",             "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NSelMuons",     "Num. of selected muon",     "N(selected #mu)",    "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseMuIsoMu", "Num. of loose muon",        "N(loose #mu)",       "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseElIsoMu", "Num. of loose electron",    "N(loose e)",         "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NElectrons",    "Num. of electron",          "N(e)",               "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NSelElectrons", "Num. of selected electron", "N(selected e)",      "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseMuIsoEl", "Num. of loose muon",        "N(loose #mu)",       "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseElIsoEl", "Num. of loose electron",    "N(loose e)",         "Events", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_Channel",       "",                          "",                   "Evetns", "",    "", 1,    0,   1   );

    h1.addNewTH1( "Gen_PID",           "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo1",      "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo2",      "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo1",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo2",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo1",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo2",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo1",    "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo2",    "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa1",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa2",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa1",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa2",     "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_NbQ",           "",                          "",                   "Evetns", "",    "", 10,   0,   10 );
    h1.addNewTH1( "Gen_NlightQ",       "",                          "",                   "Evetns", "",    "", 10,   0,   10 );
    h1.addNewTH1( "Gen_NchargeLep",    "",                          "",                   "Evetns", "",    "", 10,   0,   10 );
    h1.addNewTH1( "Gen_Nneutrinos",    "",                          "",                   "Evetns", "",    "", 10,   0,   10 );
    h1.addNewTH1( "Gen_O2_highPt",     "O2",                        "",                   "Events", "",    "", 40,  -2,   2  );
    h1.addNewTH1( "Gen_O2Asym_highPt", "A_{O2}",                    "",                   "Events", "",    "",  2,   0,   2  );
    h1.addNewTH1( "Gen_highPtLQ_PDG",  "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtLQ_Mo1",  "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtLQ_Mo2",  "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtLep_PDG", "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtLep_Mo1", "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtLep_Mo2", "",                          "",                   "Evetns", "",    "", 60, -30,   30 );
    h1.addNewTH1( "Gen_highPtBQ_Pair", "",                          "",                   "Events", "",    "",  2,   0,   2  );
    h1.addNewTH1( "Gen_highPt_HasO2",  "No b bbar and mo select",   "",                   "Evetns", "",    "",  2,   0,   2  );
    h1.addNewTH1( "Gen_O2_ideal",      "O2",                        "",                   "Events", "",    "", 40,  -2,   2  );
    h1.addNewTH1( "Gen_O2Asym_ideal",  "A_{O2}",                    "",                   "Events", "",    "",  2,   0,   2  );
    h1.addNewTH1( "Gen_ideal_HasO2",   "O2 in theory",              "",                   "Evetns", "",    "",  2,   0,   2  );

    h1.CreateTH1( fs );
    h1.Sumw2();

    setObservableHist( h1.GetTH1("Evt_O7Asym"),    "O_{7}");
    setObservableHist( h1.GetTH1("Evt_O7Asym_Mu"), "O_{7}");
    setObservableHist( h1.GetTH1("Evt_O7Asym_El"), "O_{7}");
    setObservableHist( h1.GetTH1("Evt_O2Asym"),    "O_{2}");
    setObservableHist( h1.GetTH1("Evt_O2Asym_Mu"), "O_{2}");
    setObservableHist( h1.GetTH1("Evt_O2Asym_El"), "O_{2}");
    setObservableHist( h1.GetTH1("Gen_O2Asym_highPt"), "O_{2}");
    setObservableHist( h1.GetTH1("Gen_O2Asym_ideal"),  "O_{2}");

    // Create TH2D
    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.ClearTH2Info();
    h2.addNewTH2("TH2_O2_vs_LepCharge",     "", "", "O_{2}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepCharge_Mu",  "", "", "O_{2}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepCharge_El",  "", "", "O_{2}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPz",         "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPz_Mu",      "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPz_El",      "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPt",         "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPt_Mu",      "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPt_El",      "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepEta",        "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepEta_Mu",     "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepEta_El",     "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPhi",        "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPhi_Mu",     "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_LepPhi_El",     "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPz",     "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPz_Mu",  "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPz_El",  "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPt",     "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPt_Mu",  "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPt_El",  "", "", "O_{2}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetEta",    "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetEta_Mu", "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetEta_El", "", "", "O_{2}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPhi",    "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPhi_Mu", "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O2_vs_HardJetPhi_El", "", "", "O_{2}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepCharge",     "", "", "O_{7}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepCharge_Mu",  "", "", "O_{7}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepCharge_El",  "", "", "O_{7}", "", "", 3,   -1,   2,   40, -2, 2);
    h2.CreateTH2( fs );
    h2.Sumw2();

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
    chain_  = new TChain(inputTTree_.c_str());

    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
    }

    EvtInfo.Register(chain_);
    GenInfo.Register(chain_);
    VtxInfo.Register(chain_);
    JetInfo.Register(chain_,"PFJetInfo");
    LepInfo.Register(chain_,"PFLepInfo");

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void SemiLeptanicResultsCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
    using namespace edm;
    using namespace std;

    if(  chain_ == 0 ) return;

    cout<<">> [INFO] Starting analysis loop with "<<maxEvents_<<" events..."<<endl;

    TVector3 ax, ay, az;
    ax.SetXYZ(1, 0, 0);
    ay.SetXYZ(0, 1, 0);
    az.SetXYZ(0, 0, 1);

    // Create all selectors

    //SelectorVertex  VtxSelection( selPars_Vertex_, Debug_ );

    SelectorJet     JetSelection(     selPars_Jet_, Debug_ );
    SelectorJet    BJetSelection(    selPars_BJet_, Debug_ );
    SelectorJet NonBJetSelection( selPars_NonBJet_, Debug_ );

    SelectorMuon MuonSelectionTight(   selPars_TightMuon_, Debug_ );
    SelectorMuon MuonSelectionLoose( selPars_LooseLepton_, Debug_ );

    SelectorElectron ElectronSelectionLoose(   selPars_LooseLepton_, Debug_ );
    SelectorElectron ElectronSelectionTight( selPars_TightElectron_, Debug_ );

    // Loop evetns
    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);

        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        h1.GetTH1("Evt_Events")->Fill(1);

        // ---- * Checking generator info
        vector<GenParticle> particles; 
        vector<GenParticle>   bquarks,   lquarks; 
        vector<GenParticle> chargeLep, neutrinos; // Only keep electron and muon
        for( int idx=0; idx < GenInfo.Size; idx++ )
        {
            GenParticle particle( GenInfo, idx );
            if( particle.Status == 3 )
            { 
                particles.push_back( particle );
                h1.GetTH1("Gen_PID")->Fill( particle.PdgID );  
                if( particle.PdgID == -24 )
                {
                    h1.GetTH1("Gen_PID_wpDa1")->Fill( particle.Da1PdgID );   
                    h1.GetTH1("Gen_PID_wpDa2")->Fill( particle.Da2PdgID );   
                }
                if( particle.PdgID == 24 )
                {
                    h1.GetTH1("Gen_PID_wmDa1")->Fill( particle.Da1PdgID );   
                    h1.GetTH1("Gen_PID_wmDa2")->Fill( particle.Da2PdgID );   
                }
                if( abs(particle.PdgID) == 5 )
                { 
                    bquarks.push_back( particle ); 
                    h1.GetTH1("Gen_PID_bMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_bMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID)  < 5 )
                { 
                    lquarks.push_back( particle );
                    h1.GetTH1("Gen_PID_lqMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lqMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) == 11 || abs(particle.PdgID) == 13 )
                {
                    chargeLep.push_back( particle );
                    h1.GetTH1("Gen_PID_lepMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lepMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) == 12 || abs(particle.PdgID) == 14 )
                { 
                    neutrinos.push_back( particle );
                    h1.GetTH1("Gen_PID_nuMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_nuMo2")->Fill( particle.Mo2PdgID );   
                }
            }
        }
        h1.GetTH1("Gen_NbQ")       ->Fill( bquarks.size()   ); 
        h1.GetTH1("Gen_NlightQ")   ->Fill( lquarks.size()   ); 
        h1.GetTH1("Gen_NchargeLep")->Fill( chargeLep.size() ); 
        h1.GetTH1("Gen_Nneutrinos")->Fill( neutrinos.size() );
        
        // Gen level O2
        int bqsize = bquarks.size(); 
        int lqsize = lquarks.size(); 
        int lpsize = chargeLep.size(); 
        if( bqsize >= 2 && lqsize >= 1 && lpsize >=1 ) 
        {
            vector<GenParticle> bjets, bbarjets;
            // blined: ignore b or bbar
            for( int bq1=0; bq1<bqsize; bq1++)
            {
                if( bquarks[bq1].PdgID ==  5 )   bjets.push_back( bquarks[bq1] ); 
                if( bquarks[bq1].PdgID == -5 ) bbarjets.push_back( bquarks[bq1] ); 
            }

            // Choose high pT
            GenParticle b1, b2, lq, lep;
            get2HighPtObject( bquarks, b1, b2 );
            getHighPtObject( lquarks,   lq  );
            getHighPtObject( chargeLep, lep );
            double O2 = Obs2( lep.P3, lq.P3, b1.P3, b2.P3 );

            fillObservableHist( h1.GetTH1("Gen_O2Asym_highPt"), O2, "O_{2}");
            h1.GetTH1("Gen_O2_highPt")    ->Fill(  O2/Owrt_    );
            h1.GetTH1("Gen_highPtLQ_PDG") ->Fill(  lq.PdgID    );
            h1.GetTH1("Gen_highPtLQ_Mo1") ->Fill(  lq.Mo1PdgID );
            h1.GetTH1("Gen_highPtLQ_Mo2") ->Fill(  lq.Mo2PdgID );
            h1.GetTH1("Gen_highPtLep_PDG")->Fill( lep.PdgID    );
            h1.GetTH1("Gen_highPtLep_Mo1")->Fill( lep.Mo1PdgID );
            h1.GetTH1("Gen_highPtLep_Mo2")->Fill( lep.Mo2PdgID );
            if( b1.PdgID/b2.PdgID == -1 )
            { 
                h1.GetTH1("Gen_highPtBQ_Pair")->Fill(1);
            }
            else 
            {    
                h1.GetTH1("Gen_highPtBQ_Pair")->Fill(0);
            }
            h1.GetTH1("Gen_highPt_HasO2")->Fill(1);
            
            // ideal
            if( bjets.size() >= 1 && bbarjets.size() >=1 )
            {
                GenParticle bjet, bbarjet, isolep, hardjet;
                if( getHighPtSelectMo(     bjets,    bjet, 6  ) && // b from top
                    getHighPtSelectMo(  bbarjets, bbarjet, 6  ) && // b from top
                    getHighPtSelectMo(   lquarks, hardjet, 24 ) && // non-b jet from W
                    getHighPtSelectMo( chargeLep,  isolep, 24 ))   // lepton from W
                {
                    O2 = Obs2( isolep.P3, hardjet.P3, bbarjet.P3, bjet.P3 );
                    h1.GetTH1("Gen_O2_ideal")->Fill( O2/Owrt_ );
                    h1.GetTH1("Gen_ideal_HasO2")->Fill(1);
                    fillObservableHist( h1.GetTH1("Gen_O2Asym_ideal"), O2, "O_{2}");
                }
            }  
            else
            {
                h1.GetTH1("Gen_ideal_HasO2")->Fill(0);
            }  
        } 
        else
        {
            h1.GetTH1("Gen_ideal_HasO2") ->Fill(0);
            h1.GetTH1("Gen_highPt_HasO2")->Fill(0);
        }   

        // ---- * Reconstruction info
        //* Jet selection
        vector<Jet> JetColSelected, BJetCol, nonBJetCol;
        for( int idx=0; idx < JetInfo.Size; idx++)
        {
            Jet jet( JetInfo, idx );
            if( JetSelection.isPass(jet)    ) JetColSelected.push_back(jet);
            if( BJetSelection.isPass(jet)   ) BJetCol.push_back(jet);
            if( NonBJetSelection.isPass(jet)) nonBJetCol.push_back(jet);
        }
        h1.GetTH1("Evt_NJets")->Fill(JetInfo.Size);
        h1.GetTH1("Evt_NSelJets")->Fill(JetColSelected.size());
        h1.GetTH1("Evt_NbJets")->Fill(BJetCol.size());
        if( JetColSelected.size() != ( nonBJetCol.size()+BJetCol.size()) ) std::cout<<">> [WARNING] JetColSelected.size() "<<JetColSelected.size()<<" != ( nonBJetCol.size() "<<nonBJetCol.size()<<" + BJetCol.size() "<<BJetCol.size()<<" )"<<std::endl;

        //* Lepton selection
        vector<Lepton> MuColTight, MuColLoose_MuChannel, ElColLoose_MuChannel;
        vector<Lepton> ElColTight, MuColLoose_ElChannel, ElColLoose_ElChannel;
        for( int idx=0; idx < LepInfo.Size; idx++)
        {
            Lepton lepton( LepInfo, idx );
            // Electron selections
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
            // Muon selections
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

        //* Event selection
        //Jet hardJet, bjet1, bjet2;
        //Lepton isoMu, isoEl;
        //bool isGoodMuonEvt(false), isGoodElectronEvt(false);
    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void SemiLeptanicResultsCheck::endJob(){
    std::cout<<">> [INFO] End of Job!"<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SemiLeptanicResultsCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SemiLeptanicResultsCheck);
