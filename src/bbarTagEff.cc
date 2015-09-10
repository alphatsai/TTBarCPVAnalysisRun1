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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TopCandidate.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/GenParticle.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/bbarTagEff.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorElectron.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorMuon.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorJet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SelectorVertex.h" 

//
// constructors and destructor
//
bbarTagEff::bbarTagEff(const edm::ParameterSet& iConfig) : 
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
    dR_Matching_(           iConfig.getParameter<double>("dR_Matching")),
    dR_IsoLeptonFromJets_(  iConfig.getParameter<double>("dR_IsoLeptonFromJets")),
    Owrt_(                  iConfig.getParameter<double>("Owrt")),
    //isSignal_(              iConfig.getParameter<bool>("isSingal"))
    Debug_(                 iConfig.getParameter<bool>("Debug"))
{}

bbarTagEff::~bbarTagEff()
{ 
    delete chain_;
}

// ------------ Other function -------------
    template <class TH1>
void bbarTagEff::integralFromLowerBins( TH1* h_in, TH1* h_out )
{
    if( h_in->GetNbinsX() != h_out->GetNbinsX() ) printf(">> [ERROR] Deffirent bin size between h_in(%d) and h_out(%d)\n", h_in->GetNbinsX(), h_out->GetNbinsX() ); 
    int minBin = 1;
    int maxBin = h_in->GetNbinsX();
    double sum   = 0;
    double sumw2 = 0;
    for( int b=minBin; b<=maxBin; b++ )
    {
        sum   += h_in->GetBinContent(b);
        sumw2 += h_in->GetBinError(b)*h_in->GetBinError(b);
        h_out->SetBinContent(b, sum);
        h_out->SetBinError(b, sqrt(sumw2));
    }
}

    template <class Object>
bool bbarTagEff::getHighPtSelectMo( vector<Object> col, Object &obj, int mo )
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
    if( matched ) obj=col[o1];
    else std::cout<<"[WARING] Can't find matched particle in bbarTagEff::getHighPtSelectMo"<<std::endl;
    return matched;
}

// ------------ method called once each job just before starting event loop  ------------
void bbarTagEff::beginJob()
{
    // Create TH1D
    h1 = TH1InfoClass<TH1D>(Debug_);
    h1.ClearTH1Info();
    h1.addNewTH1( "Evt_O3",                "O3",                      "O_{3}",             "Events", "", "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O3_Mu",             "O3",                      "O_{3}",             "Events", "", "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O3_El",             "O3",                      "O_{3}",             "Events", "", "", 40,  -2,   2   );
    h1.addNewTH1( "Evt_O3Asym",            "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_Mu",         "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );
    h1.addNewTH1( "Evt_O3Asym_El",         "A_{O3}",                  "",                  "Events", "", "",  2,  -1,   1   );

    h1.addNewTH1( "Evt_Events",            "",                        "",                  "Events", "", "",  1,   1,   2   );
    h1.addNewTH1( "Evt_Channel",           "",                        "",                  "Events", "", "",  4,   0,   4   );
    h1.addNewTH1( "Evt_NJets",             "Num. of jets",            "N(j)",              "Events", "", "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NSelJets",          "Num. of selected jets",   "N(selected j)",     "Events", "", "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NbJets",            "Num. of b-jets",          "N(B-tagged j)",     "Events", "", "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NnonbJets",         "Num. of non b-jets",      "N(non B-tagged j)", "Events", "", "", 20,   0,   20  );
    h1.addNewTH1( "Evt_NTightMuons",       "Num. of tight muon",      "N(tight #mu)",      "Events", "", "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseMuIsoMu",     "Num. of loose muon",      "N(loose #mu)",      "Events", "", "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseMuIsoEl",     "Num. of loose muon",      "N(loose #mu)",      "Events", "", "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NTightElectrons",   "Num. of tight electron",  "N(tight e)",        "Events", "", "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseElIsoMu",     "Num. of loose electron",  "N(loose e)",        "Events", "", "", 10,   0,   10  );
    h1.addNewTH1( "Evt_NLooseElIsoEl",     "Num. of loose electron",  "N(loose e)",        "Events", "", "", 10,   0,   10  );

    h1.addNewTH1( "Evt_Top_Hadronic_Chi2", "",                        "#chi^{2}",          "Events", "", "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass", "",                        "Mass",              "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt",   "",                        "Pt",                "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta",  "",                        "Eta",               "Events", "", "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi",  "",                        "Phi",               "Events", "", "", 65,  -3.2, 3.2 );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt",   "",                        "Mass",              "Events", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi",  "",                        "Phi",               "Events", "", "", 65,  -3.2, 3.2 );

    h1.addNewTH1( "Evt_tMatched_Chi2",      "",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2",      "",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_lbbar","",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_bl",   "",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_ll",   "",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_bbarb","",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_bbarl","",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_lb",   "",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2_other","",                        "#chi^{2}",          "Events", "", "", 200,  0, 200 );

    h1.addNewTH1( "Evt_bJet_PdgID",        "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Evt_bbarJet_PdgID",     "",                        "",                  "Evetns", "", "", 60, -30,   30  );

    h1.addNewTH1( "Evt_Top_Hadronic_Chi2CutInt", "",                 "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Chi2CutEff", "",                 "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2CutInt",     "",                 "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_bMatched_Chi2CutEff",     "",                 "#chi^{2}",          "Events", "", "", 200,  0, 200 );
    h1.addNewTH1( "Evt_Chi2CutEff",              "",                 "#chi^{2}",          "Events", "", "", 200,  0, 200 );

    h1.addNewTH1( "Gen_PID",               "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_bMo1",          "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_bMo2",          "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_lqMo1",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_lqMo2",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_nuMo1",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_nuMo2",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_lepMo1",        "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_lepMo2",        "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_wpDa1",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_wpDa2",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_wmDa1",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_PID_wmDa2",         "",                        "",                  "Evetns", "", "", 60, -30,   30  );
    h1.addNewTH1( "Gen_NbQ",               "",                        "",                  "Evetns", "", "", 10,   0,   10  );
    h1.addNewTH1( "Gen_NlightQ",           "",                        "",                  "Evetns", "", "", 10,   0,   10  );
    h1.addNewTH1( "Gen_NchargeLeps",       "",                        "",                  "Evetns", "", "", 10,   0,   10  );
    h1.addNewTH1( "Gen_Nneutrinos",        "",                        "",                  "Evetns", "", "", 10,   0,   10  );
    h1.CreateTH1( fs );
    h1.Sumw2();

    // Create TH2D
    h2 = TH2InfoClass<TH2D>(Debug_);
    h2.ClearTH2Info();
    h2.addNewTH2( "TH2_Chi2_vs_dR_GenPar_TopCan", "", "", "#chi^{2}", "dR(Gen_t,t_reco)", "", 200, 0, 200, 300,  0, 3 );
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

    chain_->SetBranchAddress("EvtInfo.RunNo",     &RunNo_         );
    chain_->SetBranchAddress("EvtInfo.EvtNo",     &EvtNo_         );
    chain_->SetBranchAddress("EvtInfo.BxNo",      &BxNo_          );
    chain_->SetBranchAddress("EvtInfo.LumiNo",    &LumiNo_        );
    chain_->SetBranchAddress("EvtInfo.isMuonEvt", &isMuonEvt_     );    
    chain_->SetBranchAddress("EvtInfo.isEleEvt",  &isElectronEvt_ );    

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void bbarTagEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
    using namespace edm;
    using namespace std;

    if(  chain_ == 0 ) return;

    cout<<">> [INFO] Starting analysis loop with "<<maxEvents_<<" events..."<<endl;

    TVector3 ax, ay, az;
    ax.SetXYZ(1, 0, 0);
    ay.SetXYZ(0, 1, 0);
    az.SetXYZ(0, 0, 1);

    //// ----- * Create all selectors
    SelectorJet     JetSelection(     selPars_Jet_, Debug_ );
    SelectorJet    BJetSelection(    selPars_BJet_, Debug_ );
    SelectorJet NonBJetSelection( selPars_NonBJet_, Debug_ );

    SelectorMuon MuonSelectionTight(   selPars_TightMuon_, Debug_ );
    SelectorMuon MuonSelectionLoose( selPars_LooseLepton_, Debug_ );

    SelectorElectron ElectronSelectionLoose(   selPars_LooseLepton_, Debug_ );
    SelectorElectron ElectronSelectionTight( selPars_TightElectron_, Debug_ );

    //// ----- * Loop evetns
    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);

        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        h1.GetTH1("Evt_Events")->Fill(1);

             if(  isMuonEvt_ &&  isElectronEvt_ ) h1.GetTH1("Evt_Channel")->Fill(3);
        else if(  isMuonEvt_ && !isElectronEvt_ ) h1.GetTH1("Evt_Channel")->Fill(2);
        else if( !isMuonEvt_ &&  isElectronEvt_ ) h1.GetTH1("Evt_Channel")->Fill(1);
        else if( !isMuonEvt_ && !isElectronEvt_ ) h1.GetTH1("Evt_Channel")->Fill(0);

        // ---- * Checking generator info
        vector<GenParticle>  particles; 
        vector<GenParticle>    bquarks,   lquarks, quarks1, quarks2, tops; 
        vector<GenParticle> chargeLeps, neutrinos; // Only keep electron and muon
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
                if( abs(particle.PdgID) == 6 )
                {
                    tops.push_back( particle ); 
                }
                if( abs(particle.PdgID) == 5 )
                { 
                    quarks1.push_back( particle ); 
                    quarks2.push_back( particle ); 
                    bquarks.push_back( particle ); 
                    h1.GetTH1("Gen_PID_bMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_bMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) <= 4 || particle.PdgID == 21 )
                {
                    if( !checkMo( particle, 5 ) && !checkMo( particle, -5 )) quarks2.push_back( particle ); 
                    quarks1.push_back( particle );
                    lquarks.push_back( particle );
                    h1.GetTH1("Gen_PID_lqMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lqMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) == 11 || abs(particle.PdgID) == 13 || abs(particle.PdgID) == 15)
                {
                    chargeLeps.push_back( particle );
                    h1.GetTH1("Gen_PID_lepMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lepMo2")->Fill( particle.Mo2PdgID );  
                    //printf(" ID: %d, Mo: %d %d\n", particle.PdgID, particle.Mo1PdgID, particle.Mo2PdgID );
                    //int mo1 = particle.Mo1; 
                    //int mo2 = particle.Mo2; 
                    //printf("    Mo1: %3d:%3d, status: %3d, Da: %3d %3d\n", particle.Mo1PdgID, GenInfo.PdgID[mo1], GenInfo.Status[mo1], GenInfo.Da1PdgID[mo1], GenInfo.Da2PdgID[mo1]); 
                    //printf("                               Da: %3d %3d (%d %d)\n", GenInfo.PdgID[GenInfo.Da1[mo1]], GenInfo.PdgID[GenInfo.Da2[mo1]], GenInfo.Status[GenInfo.Da1[mo1]], GenInfo.Status[GenInfo.Da2[mo1]]); 
                    //printf("    Mo2: %3d:%3d, status: %3d, Da: %3d %3d\n", particle.Mo2PdgID, GenInfo.PdgID[mo2], GenInfo.Status[mo2], GenInfo.Da1PdgID[mo2], GenInfo.Da2PdgID[mo2]); 
                    //printf("                               Da: %3d %3d (%d %d)\n", GenInfo.PdgID[GenInfo.Da1[mo2]], GenInfo.PdgID[GenInfo.Da2[mo2]], GenInfo.Status[GenInfo.Da1[mo2]], GenInfo.Status[GenInfo.Da2[mo2]]); 
                }
                if( abs(particle.PdgID) == 12 || abs(particle.PdgID) == 14 || abs(particle.PdgID) == 16 )
                { 
                    neutrinos.push_back( particle );
                    h1.GetTH1("Gen_PID_nuMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_nuMo2")->Fill( particle.Mo2PdgID );   
                }
            }
        }
        h1.GetTH1("Gen_NbQ"        )->Fill( bquarks.size()    ); 
        h1.GetTH1("Gen_NlightQ"    )->Fill( lquarks.size()    ); 
        h1.GetTH1("Gen_NchargeLeps")->Fill( chargeLeps.size() ); 
        h1.GetTH1("Gen_Nneutrinos" )->Fill( neutrinos.size()  );
        
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

        h1.GetTH1("Evt_NJets"          )->Fill( JetInfo.Size               );
        h1.GetTH1("Evt_NSelJets"       )->Fill( JetColSelected.size()      );
        h1.GetTH1("Evt_NbJets"         )->Fill( BJetCol.size()             );
        h1.GetTH1("Evt_NnonbJets"      )->Fill( nonBJetCol.size()          );
        h1.GetTH1("Evt_NTightMuons"    )->Fill( MuColTight.size()          );
        h1.GetTH1("Evt_NLooseMuIsoMu"  )->Fill( MuColLoose_MuChannel.size());
        h1.GetTH1("Evt_NLooseMuIsoEl"  )->Fill( MuColLoose_ElChannel.size());
        h1.GetTH1("Evt_NTightElectrons")->Fill( ElColTight.size()          );
        h1.GetTH1("Evt_NLooseElIsoEl"  )->Fill( ElColLoose_ElChannel.size());
        h1.GetTH1("Evt_NLooseElIsoMu"  )->Fill( ElColLoose_MuChannel.size());

        if( isMuonEvt_     && MuColTight.size() != 1 ) printf(">> [WARNING] Size of tight Mu is not 1 in Muon channle: %ld\n", MuColTight.size());
        if( isElectronEvt_ && ElColTight.size() != 1 ) printf(">> [WARNING] Size of tight El is not 1 in Electron channle: %ld\n", ElColTight.size());
        if( BJetCol.size() != 2                      ) printf(">> [WARNING] Size of b-jets are not 2: %ld\n",BJetCol.size());
        if( BJetCol.size() + nonBJetCol.size() < 4   ) printf(">> [WARNING] Size of all selected jets are above 3: %ld\n",BJetCol.size()+nonBJetCol.size());

        Jet hardJet1, hardJet2;
        Jet b_jet;
        Jet bbar_jet;
        Lepton isoLep; 
        // Lable isolated lepton
             if(  isMuonEvt_ && !isElectronEvt_ ) isoLep = MuColTight[0];
        else if( !isMuonEvt_ &&  isElectronEvt_ ) isoLep = ElColTight[0];
        else printf(">> [WARNING] Ambiguous isolated lepton, either eletron or muon\n");

        // Lable the hardest non_bjet
        get2HighPtObject( nonBJetCol, hardJet1, hardJet2 );

        // Distinguish hadronic-top and leptonic-top's b-jet by chi^2
        int iHadronic(0), iLeptonic(1);
        vector<TopCandidate> topPair;         
        const int sizeNonBJetCol = nonBJetCol.size();
        float chi2 = +1E10;
        int topjet1(-1), topjet2(-1), hadronicTopbjet(-1), leptonicTopbjet(-1); 
        for( int ij1=1; ij1<sizeNonBJetCol; ij1++){
            for( int ij2=ij1-1; ij2<ij1; ij2++){
                for( int bj=0; bj<2; bj++)
                {
                    float chi2_ = getChi2( nonBJetCol[ij1], nonBJetCol[ij2], BJetCol[bj] );
                    if( chi2_ < chi2 )
                    { 
                        topjet1  = ij1;
                        topjet2  = ij2;
                        hadronicTopbjet = bj;
                        chi2 = chi2_;
                    }
                }
            }
        }       

        leptonicTopbjet = ( hadronicTopbjet==0 )? 1:0;
        for( int i=0; i<2; i++){
            TopCandidate top_;
            topPair.push_back(top_);
        }
        topPair[iHadronic].Fill( BJetCol[hadronicTopbjet], nonBJetCol[topjet1], nonBJetCol[topjet2] );
        topPair[iLeptonic].Fill( BJetCol[leptonicTopbjet], isoLep, EvtInfo.PFMET, EvtInfo.PFMETPhi  );

        h1.GetTH1("Evt_Top_Hadronic_Chi2")->Fill( chi2               );
        h1.GetTH1("Evt_Top_Hadronic_Mass")->Fill( topPair[iHadronic].Mass  );
        h1.GetTH1("Evt_Top_Hadronic_Pt"  )->Fill( topPair[iHadronic].Pt    );
        h1.GetTH1("Evt_Top_Hadronic_Eta" )->Fill( topPair[iHadronic].Eta   );
        h1.GetTH1("Evt_Top_Hadronic_Phi" )->Fill( topPair[iHadronic].Phi   );
        h1.GetTH1("Evt_Top_Leptonic_Mt"  )->Fill( topPair[iLeptonic].MassT );
        h1.GetTH1("Evt_Top_Leptonic_Phi" )->Fill( topPair[iLeptonic].Phi   );

        if( isoLep.Charge < 0 ) // tbar->bbar+w-, w- -> l- v
        {
            b_jet    = BJetCol[hadronicTopbjet];
            bbar_jet = BJetCol[leptonicTopbjet];
        }
        else if( isoLep.Charge > 0 ) //t->b+w+, w+ -> l+ v
        {
            b_jet    = BJetCol[leptonicTopbjet];
            bbar_jet = BJetCol[hadronicTopbjet];
        }
        else
        { std::cout<<">> [ERROR] There an nuetral lepton!? "<<std::endl; }

        vector<Jet> jetCandidates;
        jetCandidates.push_back(    b_jet ); //0
        jetCandidates.push_back( bbar_jet ); //1
        jetCandidates.push_back( hardJet1 ); //2
        jetCandidates.push_back( hardJet2 ); //3

        // 1* Gen matching to top
        if( tops.size() != 0 )
        {
            GenParticle GenParticle_matchedHadronicTop;
            matchObject( topPair[iHadronic], GenParticle_matchedHadronicTop, tops, 1E+8);

            double dR = topPair[iHadronic].P4.DeltaR(GenParticle_matchedHadronicTop.P4); 
            h2.GetTH2("TH2_Chi2_vs_dR_GenPar_TopCan")->Fill( chi2, dR );           

            if( dR < 0.5 )
            {
                // Check whether the matched top is hadronic decay
                int iW = (abs(GenParticle_matchedHadronicTop.Da1PdgID) == 24 )? GenParticle_matchedHadronicTop.Da1:GenParticle_matchedHadronicTop.Da2;
                GenParticle W( GenInfo, iW );
                if( abs(W.Da1PdgID) <=  5 || abs(W.Da2PdgID) <= 5 ) 
                {
                    h1.GetTH1("Evt_tMatched_Chi2")->Fill( chi2 );            
                }
            }
        }

        // 2* Gen matching to b and bbar jets 
        vector<GenParticle> GenParticle_matchedBJet;
        matchMultiObject( jetCandidates, quarks2, GenParticle_matchedBJet );

        h1.GetTH1("Evt_bJet_PdgID"   )->Fill( GenParticle_matchedBJet[0].PdgID ); 
        h1.GetTH1("Evt_bbarJet_PdgID")->Fill( GenParticle_matchedBJet[1].PdgID );
 
        if( GenParticle_matchedBJet[0].PdgID == 5 && GenParticle_matchedBJet[1].PdgID == -5 )
            h1.GetTH1("Evt_bMatched_Chi2")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID !=  5 && GenParticle_matchedBJet[1].PdgID == -5 )
            h1.GetTH1("Evt_bMatched_Chi2_lbbar")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID ==  5 && GenParticle_matchedBJet[1].PdgID != -5 )
            h1.GetTH1("Evt_bMatched_Chi2_bl")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID == -5 && GenParticle_matchedBJet[1].PdgID == 5 )
            h1.GetTH1("Evt_bMatched_Chi2_bbarb")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID == -5 && GenParticle_matchedBJet[1].PdgID != 5 )
            h1.GetTH1("Evt_bMatched_Chi2_bbarl")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID != -5 && GenParticle_matchedBJet[1].PdgID == 5 )
            h1.GetTH1("Evt_bMatched_Chi2_lb")->Fill( chi2 );
        else if( GenParticle_matchedBJet[0].PdgID !=  5 && GenParticle_matchedBJet[1].PdgID != -5 )
            h1.GetTH1("Evt_bMatched_Chi2_ll")->Fill( chi2 );
        else
            h1.GetTH1("Evt_bMatched_Chi2_other")->Fill( chi2 );
       
    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void bbarTagEff::endJob(){
    integralFromLowerBins( h1.GetTH1("Evt_bMatched_Chi2"),     h1.GetTH1("Evt_bMatched_Chi2CutInt")    );
    integralFromLowerBins( h1.GetTH1("Evt_Top_Hadronic_Chi2"), h1.GetTH1("Evt_Top_Hadronic_Chi2CutInt"));
    h1.GetTH1("Evt_bMatched_Chi2CutEff"    )->Add( h1.GetTH1("Evt_bMatched_Chi2CutInt")    );
    h1.GetTH1("Evt_Top_Hadronic_Chi2CutEff")->Add( h1.GetTH1("Evt_Top_Hadronic_Chi2CutInt"));

    double scale = 1/h1.GetTH1("Evt_Top_Hadronic_Chi2")->Integral();
    h1.GetTH1("Evt_Top_Hadronic_Chi2CutEff")->Scale(scale);
    h1.GetTH1("Evt_bMatched_Chi2CutEff")->Divide( h1.GetTH1("Evt_Top_Hadronic_Chi2CutInt") );

    h1.GetTH1("Evt_Chi2CutEff")->Add( h1.GetTH1("Evt_Top_Hadronic_Chi2CutEff") );
    h1.GetTH1("Evt_Chi2CutEff")->Multiply( h1.GetTH1("Evt_bMatched_Chi2CutEff"));
    
    std::cout<<">> [INFO] End of Job!"<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void bbarTagEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbarTagEff);
