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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TopCandidate.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicAnalysis.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorElectron.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorMuon.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorJet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorVertex.h" 

//
// constructors and destructor
//
SemiLeptanicAnalysis::SemiLeptanicAnalysis(const edm::ParameterSet& iConfig) : 
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
    dR_IsoLeptonFromJets_( iConfig.getParameter<double>("dR_IsoLeptonFromJets")),
    Owrt_(                  iConfig.getParameter<double>("Owrt")),
    NJets_(                 iConfig.getParameter<int>("NJets")),
    Debug_(                 iConfig.getParameter<bool>("Debug")),
    isSkim_(                iConfig.getParameter<bool>("IsSkim")),
    doSaveTree_(            iConfig.getParameter<bool>("DoSaveTree"))
{
    if( inputTTree_.compare("Skim/root") == 0 ){ isSkim_=true; }
}

SemiLeptanicAnalysis::~SemiLeptanicAnalysis()
{ 
    delete chain_;
}

// ------------ Other function -------------
std::string SemiLeptanicAnalysis::int2str( int i )
{
    std::string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}
    template<class TH1>
void SemiLeptanicAnalysis::setCutFlow( TH1* h )
{
    if( isSkim_ )
    {
        h->GetXaxis()->SetBinLabel(1,"All");
        h->GetXaxis()->SetBinLabel(2,"PreSelect");
        h->GetXaxis()->SetBinLabel(3,"#geq1 goodVtx");
        h->GetXaxis()->SetBinLabel(4,"HLT");
        h->GetXaxis()->SetBinLabel(5,"#geq1 Lep");
        h->GetXaxis()->SetBinLabel(6,"1 isoLep");
        h->GetXaxis()->SetBinLabel(7,"veto(Loose #mu)");
        h->GetXaxis()->SetBinLabel(8,"veto(Loose e)");
        h->GetXaxis()->SetBinLabel(9,("#geq"+int2str(NJets_)+" Jets").c_str());
        h->GetXaxis()->SetBinLabel(10,"=2 bjets");
    }else{
        h->GetXaxis()->SetBinLabel(1,"All");
        h->GetXaxis()->SetBinLabel(2,"#geq1 goodVtx");
        h->GetXaxis()->SetBinLabel(3,"HLT");
        h->GetXaxis()->SetBinLabel(4,"#geq1 Lep");
        h->GetXaxis()->SetBinLabel(5,"1 isoLep");
        h->GetXaxis()->SetBinLabel(6,"veto(Loose #mu)");
        h->GetXaxis()->SetBinLabel(7,"veto(Loose e)");
        h->GetXaxis()->SetBinLabel(8,("#geq"+int2str(NJets_)+" Jets").c_str());
        h->GetXaxis()->SetBinLabel(9,"=2 bjets");
    }
    return ;
}

template<class TH1>
void SemiLeptanicAnalysis::setObservableHist(TH1* h, string ob ){
    h->GetXaxis()->SetBinLabel(1,(ob+"<0").c_str());
    h->GetXaxis()->SetBinLabel(2,(ob+">0").c_str());
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

    template <class Object>
void SemiLeptanicAnalysis::get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 )
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

bool SemiLeptanicAnalysis::isIsoLeptonFromJets( Lepton lepton, vector<Jet> jetCol, double dR )
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

float SemiLeptanicAnalysis::getChi2( Jet jet1, Jet jet2, Jet bjet, float M_top, float Wth_top, float M_W, float Wth_W )
{
    TLorentzVector qq_v  = jet1.P4 + jet2.P4;
    TLorentzVector bqq_v = bjet.P4 + qq_v;
    float iTop = ( bqq_v.M() - M_top )/Wth_top;
    float iW   = (  qq_v.M() - M_W   )/Wth_W;
    return iTop*iTop + iW*iW;
}

double SemiLeptanicAnalysis::Obs2( Lepton isoLep, Jet hardJet, Jet bjet1, Jet bjet2 )
{
    TVector3 O2_1v =  bjet1.P3 + bjet2.P3;
    TVector3 O2_2v = isoLep.P3.Cross( hardJet.P3 );
    double O2 = O2_1v.Dot( O2_2v );
    return O2;
}

double SemiLeptanicAnalysis::Obs7( TVector3 beam, Jet bjet1, Jet bjet2 )
{
    double O7_1z = beam.Dot( bjet1.P3 - bjet2.P3 );
    double O7_2z = beam.Dot( bjet1.P3.Cross( bjet2.P3 ));
    double O7 = O7_1z * O7_2z;
    return O7;
}
// ------------ method called once each job just before starting event loop  ------------
void SemiLeptanicAnalysis::beginJob()
{
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.addNewTH1( "Evt_O7_Mu",         "O7",                        "O_{7}",              "Events", "",    "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O7_El",         "O7",                        "O_{7}",              "Events", "",    "", 40,  -2,   2   );
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

    h1.addNewTH1( "Evt_bJet1_Pt",      "pT of b-Jet",               "p_{T}(B-tagged j)",  "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet1_M",       "Mass of b-Jet",             "Mass(B-tagged j)",   "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet1_E",       "Energy of b-Jet",           "Energy(B-tagged j)", "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet1_Eta",     "Eta of b-Jet",              "#eta(B-tagged j)",   "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bJet1_Phi",     "Phi of b-Jet",              "#phi(B-tagged j)",   "Yields", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bJet1_BTag",    "b-Jet b-tagged",            "bTag",               "Yields", "",    "", 100,  0,   1   );

    h1.addNewTH1( "Evt_bJet2_Pt",      "pT of b-Jet",               "p_{T}(B-tagged j)",  "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet2_M",       "Mass of b-Jet",             "Mass(B-tagged j)",   "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet2_E",       "Energy of b-Jet",           "Energy(B-tagged j)", "Yields", "GeV", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_bJet2_Eta",     "Eta of b-Jet",              "#eta(B-tagged j)",   "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_bJet2_Phi",     "Phi of b-Jet",              "#phi(B-tagged j)",   "Yields", "",    "", 64,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_bJet2_BTag",    "b-Jet b-tagged",            "bTag",               "Yields", "",    "", 100,  0,   1   );

    h1.addNewTH1( "Evt_Top_Hadronic_Chi2", "",                      "#Chi^{2}",           "Yields", "",    "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass", "",                      "Mass",               "Yields", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt",   "",                      "Pt",                 "Yields", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta",  "",                      "Eta",                "Yields", "",    "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi",  "",                      "Phi",                "Yields", "",    "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt",   "",                      "Mass",               "Yields", "",    "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi",  "",                      "Phi",                "Yields", "",    "", 65,  -3.2, 3.2 );

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
    h1.addNewTH1( "Evt_CutFlow_Mu",    "",                          "",                   "Evetns", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_CutFlow_El",    "",                          "",                   "Evetns", "",    "", 10,   0,   10  );
    h1.addNewTH1( "Evt_MuCut",         "isoMu:looseMu:looseEl",     "",                   "Evetns", "",    "", 7,    0,   7   );
    h1.addNewTH1( "Evt_ElCut",         "isoEl:looseMu:looseEl",     "",                   "Evetns", "",    "", 7,    0,   7   );
    h1.addNewTH1( "Evt_SameChannel",   "",                          "",                   "Evetns", "",    "", 1,    0,   1   );
    h1.addNewTH1( "Evt_Triger",        "",                          "",                   "",       "",    "", 5900, 0,   5900);

    h1.CreateTH1( fs );
    h1.Sumw2();

    setCutFlow( h1.GetTH1("Evt_CutFlow") );
    setCutFlow( h1.GetTH1("Evt_CutFlow_Mu") );
    setCutFlow( h1.GetTH1("Evt_CutFlow_El") );
    setObservableHist(h1.GetTH1("Evt_O7Asym"),    "O_{7}");
    setObservableHist(h1.GetTH1("Evt_O7Asym_Mu"), "O_{7}");
    setObservableHist(h1.GetTH1("Evt_O7Asym_El"), "O_{7}");
    setObservableHist(h1.GetTH1("Evt_O2Asym"),    "O_{2}");
    setObservableHist(h1.GetTH1("Evt_O2Asym_Mu"), "O_{2}");
    setObservableHist(h1.GetTH1("Evt_O2Asym_El"), "O_{2}");
    setLeptonSelHist(h1.GetTH1("Evt_MuCut"));
    setLeptonSelHist(h1.GetTH1("Evt_ElCut"));

    h2 = TH2InfoClass<TH2D>(Debug_);
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
    h2.addNewTH2("TH2_O7_vs_LepPz",         "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPz_Mu",      "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPz_El",      "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPt",         "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPt_Mu",      "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPt_El",      "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepEta",        "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepEta_Mu",     "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepEta_El",     "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPhi",        "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPhi_Mu",     "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_LepPhi_El",     "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPz",     "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPz_Mu",  "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPz_El",  "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPt",     "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPt_Mu",  "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPt_El",  "", "", "O_{7}", "", "", 500,  0,   500, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetEta",    "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetEta_Mu", "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetEta_El", "", "", "O_{7}", "", "", 100, -5,   5,   40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPhi",    "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPhi_Mu", "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);
    h2.addNewTH2("TH2_O7_vs_HardJetPhi_El", "", "", "O_{7}", "", "", 64,  -3.2, 3.2, 40, -2, 2);

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
    if( doSaveTree_ )
    {
        fs->cd();
        newtree_ = chain_->CloneTree(0);
    }else{
        newtree_ = fs->make<TTree>("tree", "");
        newtree_->Branch("EvtInfo.RunNo",   &RunNo_,    "EvtInfo.RunNo/I"  );
        newtree_->Branch("EvtInfo.EvtNo",   &EvtNo_,    "EvtInfo.EvtNo/L"  );
        newtree_->Branch("EvtInfo.BxNo",    &BxNo_,     "EvtInfo.BxNo/I"   );
        newtree_->Branch("EvtInfo.LumiNo",  &LumiNo_,   "EvtInfo.LumiNo/I" );
    }
    newtree_->Branch("EvtInfo.isMuonEvt",   &isMuonEvt_, "EvtInfo.isMuonEvt/I" );
    newtree_->Branch("EvtInfo.isEleEvt",    &isEleEvt_,  "EvtInfo.isEleEvt/I"  );

    if( isSkim_ )
    { 
        h1.GetTH1("Evt_CutFlow")->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow")->SetBinError(1, sqrt(allEvents));
        h1.GetTH1("Evt_CutFlow_Mu")->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow_Mu")->SetBinError(1, sqrt(allEvents));
        h1.GetTH1("Evt_CutFlow_El")->Fill("All", allEvents);
        h1.GetTH1("Evt_CutFlow_El")->SetBinError(1, sqrt(allEvents));
    }

    EvtInfo.Register(chain_);
    VtxInfo.Register(chain_);
    JetInfo.Register(chain_,"PFJetInfo");
    LepInfo.Register(chain_,"PFLepInfo");

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void SemiLeptanicAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    SelectorVertex  VtxSelection( selPars_Vertex_, Debug_ );

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

        // HLT selection
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

        // Vertex selection
        vector<Vertex> VxtColSelected;
        for( int idx=0; idx < VtxInfo.Size; idx++)
        {
            Vertex vtx( VtxInfo, idx );
            if( VtxSelection.isPass(vtx) ) VxtColSelected.push_back(vtx);
        }

        //* Jet selection
        vector<Jet> JetColSelected, BJetCol, nonBJetCol;
        for( int idx=0; idx < JetInfo.Size; idx++)
        {
            Jet jet( JetInfo, idx );
            h1.GetTH1("Jet_Pt")->Fill( jet.Pt );
            h1.GetTH1("Jet_Px")->Fill( jet.Px );
            h1.GetTH1("Jet_Py")->Fill( jet.Py );
            h1.GetTH1("Jet_Pz")->Fill( jet.Pz );
            h1.GetTH1("Jet_M")->Fill(  jet.Mass );
            h1.GetTH1("Jet_E")->Fill(  jet.Energy );
            h1.GetTH1("Jet_Eta")->Fill(  jet.Eta );
            h1.GetTH1("Jet_Phi")->Fill(  jet.Phi );
            h1.GetTH1("Jet_BTag")->Fill( jet.CombinedSVBJetTags );

            if( JetSelection.isPass(jet) )
            { 
                JetColSelected.push_back(jet);
                h1.GetTH1("SelJet_Pt")->Fill( jet.Pt );
                h1.GetTH1("SelJet_Px")->Fill( jet.Px );
                h1.GetTH1("SelJet_Py")->Fill( jet.Py );
                h1.GetTH1("SelJet_Pz")->Fill( jet.Pz );
                h1.GetTH1("SelJet_M")->Fill(  jet.Mass );
                h1.GetTH1("SelJet_E")->Fill(  jet.Energy );
                h1.GetTH1("SelJet_Eta")->Fill(  jet.Eta );
                h1.GetTH1("SelJet_Phi")->Fill(  jet.Phi );
                h1.GetTH1("SelJet_BTag")->Fill( jet.CombinedSVBJetTags );
            }
            if( BJetSelection.isPass(jet) )
            { 
                BJetCol.push_back(jet);
                h1.GetTH1("bJet_Pt")->Fill( jet.Pt );
                h1.GetTH1("bJet_Px")->Fill( jet.Px );
                h1.GetTH1("bJet_Py")->Fill( jet.Py );
                h1.GetTH1("bJet_Pz")->Fill( jet.Pz );
                h1.GetTH1("bJet_M")->Fill(  jet.Mass );
                h1.GetTH1("bJet_E")->Fill(  jet.Energy );
                h1.GetTH1("bJet_Eta")->Fill(  jet.Eta );
                h1.GetTH1("bJet_Phi")->Fill(  jet.Phi );
                h1.GetTH1("bJet_BTag")->Fill( jet.CombinedSVBJetTags );
            }
            if( NonBJetSelection.isPass(jet) ) nonBJetCol.push_back(jet);
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
        Jet hardJet, bjet1, bjet2;
        Lepton isoMu, isoEl;
        bool isGoodMuonEvt(false), isGoodElectronEvt(false);

        if( isSkim_ )
        {
            h1.GetTH1("Evt_CutFlow")->Fill("PreSelect", 1);
            h1.GetTH1("Evt_CutFlow_El")->Fill("PreSelect", 1);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("PreSelect", 1);
        }else{
            h1.GetTH1("Evt_CutFlow")->Fill("All", 1);
            h1.GetTH1("Evt_CutFlow_El")->Fill("All", 1);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("All", 1);
        }

        if( VxtColSelected.size() )
        {
            h1.GetTH1("Evt_CutFlow")->Fill("#geq1 goodVtx", 1);
            h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 goodVtx", 1);
            h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 goodVtx", 1);

            bool passElectronSel(false), passMuonSel(false);
            // Muon channel
            if( passMuonHLT )
            {
                h1.GetTH1("Evt_CutFlow_Mu")->Fill("HLT", 1);
                if( MuColTight.size() > 0 )
                { 
                    h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 Lep", 1);
                    if( MuColTight.size() == 1 )
                    {
                        isoMu = MuColTight[0];
                        h1.GetTH1("Evt_CutFlow_Mu")->Fill("1 isoLep", 1);

                        if( MuColLoose_MuChannel.size() == 0 )
                        {
                            h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose #mu)", 1);

                            if( ElColLoose_MuChannel.size() == 0 )
                            {
                                h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose e)", 1);
                                passMuonSel=true;
                            }
                        }
                    }
                }
            } // [END] Muon channel

            // Electron channel
            if( passElectronHLT )
            {
                h1.GetTH1("Evt_CutFlow_El")->Fill("HLT", 1);
                if( ElColTight.size() > 0 )
                { 
                    h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 Lep", 1);
                    if( ElColTight.size() == 1 ) 
                    {
                        isoEl=ElColTight[0];
                        h1.GetTH1("Evt_CutFlow_El")->Fill("1 isoLep", 1);

                        if( MuColLoose_ElChannel.size() == 0 )
                        {
                            h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose #mu)", 1);

                            if( ElColLoose_ElChannel.size() == 0 )
                            {
                                h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose e)", 1);
                                passElectronSel=true;
                            }
                        }
                    }
                }
            } // [END] Electron channel

            // jets and b-jet cut flow
            if( passElectronSel || passMuonSel ) 
            {
                //if( JetColSelected.size() >= NJets_ )
                int unsigned selectedJetsSize = BJetCol.size() + nonBJetCol.size();
                if( selectedJetsSize >= NJets_ )
                { 
                    if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill(("#geq"+int2str(NJets_)+" Jets").c_str(), 1);
                    if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill(("#geq"+int2str(NJets_)+" Jets").c_str(), 1);

                    if( BJetCol.size() == 2 )
                    {
                        // Lable the hardest non_bjet 
                        int j1=-1;
                        double pt1=0;
                        const int sizeNonBJetCol = nonBJetCol.size();
                        for( int i=0; i < sizeNonBJetCol; i++)
                        {
                            if( pt1 < nonBJetCol[i].Pt ){
                                j1=i;
                                pt1=nonBJetCol[i].Pt;
                            }
                        }
                        if( j1 == -1 ){ std::cout<<">>[WARNING] "<<entry<<"Doesn't find hard jet!"<<endl; }
                        hardJet = nonBJetCol[j1];

                        // Lable bjet by Pt
                        get2HighPtObject( BJetCol, bjet1, bjet2 );
                        h1.GetTH1("bJet1_Pt")->Fill(bjet1.Pt);
                        h1.GetTH1("bJet2_Pt")->Fill(bjet2.Pt);

                        // Distinguish hadronic-top and leptonic-top's b-jet by chi^2
                        TopCandidate top_hadronic, top_leptonic;         
                        float chi2 = +1E10;
                        int topjet1(-1), topjet2(-1), topbjet1(-1), topbjet2(-1); 
                        for( int ij1=1; ij1<sizeNonBJetCol; ij1++){
                            for( int ij2=ij1-1; ij2<ij1; ij2++){
                                for( int bj=0; bj<2; bj++)
                                {
                                    float chi2_ = getChi2( nonBJetCol[ij1], nonBJetCol[ij2], BJetCol[bj] );
                                    if( chi2_ < chi2 )
                                    { 
                                        topjet1  = ij1;
                                        topjet2  = ij2;
                                        topbjet1 = bj;
                                        chi2 = chi2_;
                                    }
                                }
                            }
                        }
                        topbjet2 = (topbjet1==0)? 1:0;
                        top_hadronic.Fill( BJetCol[topbjet1], nonBJetCol[topjet1], nonBJetCol[topjet2] );
                        
                        
                        // Fill cutflow hist to each channel
                        if( passMuonSel )
                        {     
                            isGoodMuonEvt=true;
                            top_leptonic.Fill( BJetCol[topbjet2], isoMu, EvtInfo.PFMET, EvtInfo.PFMETPhi );
                            h1.GetTH1("Evt_CutFlow_Mu")->Fill("=2 bjets", 1);
                        }
                        if( passElectronSel )
                        { 
                            isGoodElectronEvt=true;
                            top_leptonic.Fill( BJetCol[topbjet2], isoEl, EvtInfo.PFMET, EvtInfo.PFMETPhi );
                            h1.GetTH1("Evt_CutFlow_El")->Fill("=2 bjets", 1);
                        }
                        h1.GetTH1("Evt_Top_Hadronic_Chi2")->Fill( chi2               );
                        h1.GetTH1("Evt_Top_Hadronic_Mass")->Fill( top_hadronic.Mass  );
                        h1.GetTH1("Evt_Top_Hadronic_Pt"  )->Fill( top_hadronic.Pt    );
                        h1.GetTH1("Evt_Top_Hadronic_Eta" )->Fill( top_hadronic.Eta   );
                        h1.GetTH1("Evt_Top_Hadronic_Phi" )->Fill( top_hadronic.Phi   );
                        h1.GetTH1("Evt_Top_Leptonic_Mt"  )->Fill( top_leptonic.MassT );
                        h1.GetTH1("Evt_Top_Leptonic_Phi" )->Fill( top_leptonic.Phi   );
                    }
                }
            }//[END] Jet and bjet cutflow

            // Check events
            if( passElectronHLT || passMuonHLT ){ 
                h1.GetTH1("Evt_CutFlow")->Fill("HLT", 1);
                if( ( MuColTight.size() + ElColTight.size()) > 0 ) h1.GetTH1("Evt_CutFlow")->Fill("#geq1 Lep", 1);
                if( ( MuColTight.size() + ElColTight.size()) == 1 ){ 
                    h1.GetTH1("Evt_CutFlow")->Fill("1 isoLep", 1);
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
        } // [END] Vxt selection

        //* Fill other events plots 
        if( isGoodMuonEvt && isGoodElectronEvt ) h1.GetTH1("Evt_SameChannel")->Fill(0);
        if( isGoodMuonEvt || isGoodElectronEvt )
        {
            if( !doSaveTree_ )
            {
                RunNo_ = EvtInfo.RunNo;
                EvtNo_ = EvtInfo.EvtNo;
                BxNo_  = EvtInfo.BxNo;
                LumiNo_= EvtInfo.LumiNo;
            }
            if( isGoodMuonEvt )
            {
                isMuonEvt_ = 1;
                isEleEvt_  = 0;
            }else{
                isMuonEvt_ = 0;
                isEleEvt_  = 1;
            }
            newtree_->Fill();

            h1.GetTH1("Evt_bJet1_Pt")->Fill(bjet1.Pt);
            h1.GetTH1("Evt_bJet1_Eta")->Fill(bjet1.Eta);
            h1.GetTH1("Evt_bJet1_Phi")->Fill(bjet1.Phi);
            h1.GetTH1("Evt_bJet1_E")->Fill(bjet1.Energy);
            h1.GetTH1("Evt_bJet1_M")->Fill(bjet1.Mass);
            h1.GetTH1("Evt_bJet1_BTag")->Fill(bjet1.CombinedSVBJetTags);
            h1.GetTH1("Evt_bJet2_Pt")->Fill(bjet2.Pt);
            h1.GetTH1("Evt_bJet2_Eta")->Fill(bjet2.Eta);
            h1.GetTH1("Evt_bJet2_Phi")->Fill(bjet2.Phi);
            h1.GetTH1("Evt_bJet2_E")->Fill(bjet2.Energy);
            h1.GetTH1("Evt_bJet2_M")->Fill(bjet2.Mass);
            h1.GetTH1("Evt_bJet2_BTag")->Fill(bjet2.CombinedSVBJetTags);
            h1.GetTH1("Evt_HardJet_Pt")->Fill(hardJet.Pt);
            h1.GetTH1("Evt_HardJet_Eta")->Fill(hardJet.Eta);
            h1.GetTH1("Evt_HardJet_Phi")->Fill(hardJet.Phi);
            h1.GetTH1("Evt_HardJet_E")->Fill(hardJet.Energy);
            h1.GetTH1("Evt_HardJet_M")->Fill(hardJet.Mass);
            h1.GetTH1("Evt_HardJet_BTag")->Fill(hardJet.CombinedSVBJetTags);
        }
        //* Fill observables O7 and O2
        // -- Muon channel
        if( isGoodMuonEvt )
        {
            h1.GetTH1("Evt_isoMu_Pt")->Fill(isoMu.Pt);
            h1.GetTH1("Evt_isoMu_Et")->Fill(isoMu.Et);
            h1.GetTH1("Evt_isoMu_Eta")->Fill(isoMu.Eta);
            h1.GetTH1("Evt_isoMu_Phi")->Fill(isoMu.Phi);
            h1.GetTH1("Evt_isoMu_E")->Fill(isoMu.Energy);

            double O2 = Obs2( isoMu, hardJet, bjet1, bjet2 );
            h1.GetTH1("Evt_O2")->Fill( O2/Owrt_ );
            h1.GetTH1("Evt_O2_Mu")->Fill( O2/Owrt_ );
            if( O2 > 0 ){
                h1.GetTH1("Evt_O2Asym")->Fill("O_{2}>0",1);
                h1.GetTH1("Evt_O2Asym_Mu")->Fill("O_{2}>0",1);
            }else{
                h1.GetTH1("Evt_O2Asym")->Fill("O_{2}<0",1);
                h1.GetTH1("Evt_O2Asym_Mu")->Fill("O_{2}<0",1);
            }

            double O7 = Obs7( az, bjet1, bjet2 );
            h1.GetTH1("Evt_O7")->Fill( O7/Owrt_ );
            h1.GetTH1("Evt_O7_Mu")->Fill( O7/Owrt_ );
            if( O7 > 0 ){
                h1.GetTH1("Evt_O7Asym")->Fill("O_{7}>0",1);
                h1.GetTH1("Evt_O7Asym_Mu")->Fill("O_{7}>0",1);
            }else{
                h1.GetTH1("Evt_O7Asym")->Fill("O_{7}<0",1);
                h1.GetTH1("Evt_O7Asym_Mu")->Fill("O_{7}<0",1);
            }
            h2.GetTH2("TH2_O2_vs_LepCharge")    ->Fill( isoMu.Charge, O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepCharge_Mu") ->Fill( isoMu.Charge, O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPz")        ->Fill( isoMu.Pz,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPz_Mu")     ->Fill( isoMu.Pz,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPt")        ->Fill( isoMu.Pt,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPt_Mu")     ->Fill( isoMu.Pt,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepEta")       ->Fill( isoMu.Eta,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepEta_Mu")    ->Fill( isoMu.Eta,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPhi")       ->Fill( isoMu.Phi,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPhi_Mu")    ->Fill( isoMu.Phi,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPz")    ->Fill( hardJet.Pz,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPz_Mu") ->Fill( hardJet.Pz,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPt")    ->Fill( hardJet.Pt,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPt_Mu") ->Fill( hardJet.Pt,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetEta")   ->Fill( hardJet.Eta,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetEta_Mu")->Fill( hardJet.Eta,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPhi")   ->Fill( hardJet.Phi,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPhi_Mu")->Fill( hardJet.Phi,  O2/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepCharge")    ->Fill( isoMu.Charge, O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepCharge_Mu") ->Fill( isoMu.Charge, O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPz")        ->Fill( isoMu.Pz,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPz_Mu")     ->Fill( isoMu.Pz,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPt")        ->Fill( isoMu.Pt,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPt_Mu")     ->Fill( isoMu.Pt,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepEta")       ->Fill( isoMu.Eta,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepEta_Mu")    ->Fill( isoMu.Eta,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPhi")       ->Fill( isoMu.Phi,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPhi_Mu")    ->Fill( isoMu.Phi,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPz")    ->Fill( hardJet.Pz,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPz_Mu") ->Fill( hardJet.Pz,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPt")    ->Fill( hardJet.Pt,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPt_Mu") ->Fill( hardJet.Pt,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetEta")   ->Fill( hardJet.Eta,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetEta_Mu")->Fill( hardJet.Eta,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPhi")   ->Fill( hardJet.Phi,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPhi_Mu")->Fill( hardJet.Phi,  O7/Owrt_ );
        }
        // -- Electron channel
        if( isGoodElectronEvt )
        {
            h1.GetTH1("Evt_isoEl_Pt")->Fill(isoEl.Pt);
            h1.GetTH1("Evt_isoEl_Et")->Fill(isoEl.Et);
            h1.GetTH1("Evt_isoEl_Eta")->Fill(isoEl.Eta);
            h1.GetTH1("Evt_isoEl_Phi")->Fill(isoEl.Phi);
            h1.GetTH1("Evt_isoEl_E")->Fill(isoEl.Energy);

            double O2 = Obs2( isoEl, hardJet, bjet1, bjet2 );
            h1.GetTH1("Evt_O2")->Fill(O2/Owrt_);
            h1.GetTH1("Evt_O2_El")->Fill(O2/Owrt_);
            if( O2 > 0 ){
                h1.GetTH1("Evt_O2Asym")->Fill("O_{2}>0",1);
                h1.GetTH1("Evt_O2Asym_El")->Fill("O_{2}>0",1);
            }else{
                h1.GetTH1("Evt_O2Asym")->Fill("O_{2}<0",1);
                h1.GetTH1("Evt_O2Asym_El")->Fill("O_{2}<0",1);
            }

            double O7 = Obs7( az, bjet1, bjet2 );
            h1.GetTH1("Evt_O7")->Fill(O7/Owrt_);
            h1.GetTH1("Evt_O7_El")->Fill(O7/Owrt_);
            if( O7 > 0 ){
                h1.GetTH1("Evt_O7Asym")->Fill("O_{7}>0", 1);
                h1.GetTH1("Evt_O7Asym_El")->Fill("O_{7}>0",1);
            }else{
                h1.GetTH1("Evt_O7Asym")->Fill("O_{7}<0",1);
                h1.GetTH1("Evt_O7Asym_El")->Fill("O_{7}<0",1);
            }
            h2.GetTH2("TH2_O2_vs_LepCharge")    ->Fill( isoEl.Charge, O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepCharge_El") ->Fill( isoEl.Charge, O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPz")        ->Fill( isoEl.Pz,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPz_El")     ->Fill( isoEl.Pz,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPt")        ->Fill( isoEl.Pt,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPt_El")     ->Fill( isoEl.Pt,     O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepEta")       ->Fill( isoEl.Eta,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepEta_El")    ->Fill( isoEl.Eta,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPhi")       ->Fill( isoEl.Phi,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_LepPhi_El")    ->Fill( isoEl.Phi,    O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPz")    ->Fill( hardJet.Pz,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPz_El") ->Fill( hardJet.Pz,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPt")    ->Fill( hardJet.Pt,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPt_El") ->Fill( hardJet.Pt,   O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetEta")   ->Fill( hardJet.Eta,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetEta_El")->Fill( hardJet.Eta,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPhi")   ->Fill( hardJet.Phi,  O2/Owrt_ );
            h2.GetTH2("TH2_O2_vs_HardJetPhi_El")->Fill( hardJet.Phi,  O2/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepCharge")    ->Fill( isoEl.Charge, O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepCharge_El") ->Fill( isoEl.Charge, O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPz")        ->Fill( isoEl.Pz,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPz_El")     ->Fill( isoEl.Pz,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPt")        ->Fill( isoEl.Pt,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPt_El")     ->Fill( isoEl.Pt,     O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepEta")       ->Fill( isoEl.Eta,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepEta_El")    ->Fill( isoEl.Eta,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPhi")       ->Fill( isoEl.Phi,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_LepPhi_El")    ->Fill( isoEl.Phi,    O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPz")    ->Fill( hardJet.Pz,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPz_El") ->Fill( hardJet.Pz,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPt")    ->Fill( hardJet.Pt,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPt_El") ->Fill( hardJet.Pt,   O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetEta")   ->Fill( hardJet.Eta,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetEta_El")->Fill( hardJet.Eta,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPhi")   ->Fill( hardJet.Phi,  O7/Owrt_ );
            h2.GetTH2("TH2_O7_vs_HardJetPhi_El")->Fill( hardJet.Phi,  O7/Owrt_ );
        }
    }//// [END] entry loop 
}

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
