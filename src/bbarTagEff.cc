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

    h1.addNewTH1( "Evt_Top_Hadronic_Chi2", "",                        "#chi^{2}",          "Yields", "", "", 200,  0,   200 );
    h1.addNewTH1( "Evt_Top_Hadronic_Mass", "",                        "Mass",              "Yields", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Pt",   "",                        "Pt",                "Yields", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Hadronic_Eta",  "",                        "Eta",               "Yields", "", "", 100, -5,   5   );
    h1.addNewTH1( "Evt_Top_Hadronic_Phi",  "",                        "Phi",               "Yields", "", "", 65,  -3.2, 3.2   );
    h1.addNewTH1( "Evt_Top_Leptonic_Mt",   "",                        "Mass",              "Yields", "", "", 500,  0,   500 );
    h1.addNewTH1( "Evt_Top_Leptonic_Phi",  "",                        "Phi",               "Yields", "", "", 65,  -3.2, 3.2 );

    h1.addNewTH1( "Gen_PID",               "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo1",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_bMo2",          "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo1",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lqMo2",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo1",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_nuMo2",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo1",        "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_lepMo2",        "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa1",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wpDa2",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa1",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_PID_wmDa2",         "",                        "",                  "Evetns", "", "", 60, -30,   30 );
    h1.addNewTH1( "Gen_NbQ",               "",                        "",                  "Evetns", "", "", 10,   0,   10 );
    h1.addNewTH1( "Gen_NlightQ",           "",                        "",                  "Evetns", "", "", 10,   0,   10 );
    h1.addNewTH1( "Gen_NchargeLeps",       "",                        "",                  "Evetns", "", "", 10,   0,   10 );
    h1.addNewTH1( "Gen_Nneutrinos",        "",                        "",                  "Evetns", "", "", 10,   0,   10 );
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
        vector<GenParticle>    bquarks,   lquarks, quarks; 
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
                if( abs(particle.PdgID) == 5 )
                { 
                     quarks.push_back( particle ); 
                    bquarks.push_back( particle ); 
                    h1.GetTH1("Gen_PID_bMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_bMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID)  < 5 || particle.PdgID == 21 )
                { 
                     quarks.push_back( particle );
                    lquarks.push_back( particle );
                    h1.GetTH1("Gen_PID_lqMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lqMo2")->Fill( particle.Mo2PdgID );   
                }
                if( abs(particle.PdgID) == 11 || abs(particle.PdgID) == 13 || abs(particle.PdgID) == 15)
                {
                    chargeLeps.push_back( particle );
                    h1.GetTH1("Gen_PID_lepMo1")->Fill( particle.Mo1PdgID );   
                    h1.GetTH1("Gen_PID_lepMo2")->Fill( particle.Mo2PdgID );   
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

        Jet hardJet;
        Jet b_jet;
        Jet bbar_jet;
        Lepton isoLep; 
        // Lable isolated lepton
             if(  isMuonEvt_ && !isElectronEvt_ ) isoLep = MuColTight[0];
        else if( !isMuonEvt_ &&  isElectronEvt_ ) isoLep = ElColTight[0];
        else printf(">> [WARNING] Ambiguous isolated lepton, either eletron or muon\n");

        // Lable the hardest non_bjet
        getHighPtObject( nonBJetCol, hardJet );

        // Distinguish hadronic-top and leptonic-top's b-jet by chi^2
        TopCandidate top_hadronic, top_leptonic;         
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
        top_hadronic.Fill( BJetCol[hadronicTopbjet], nonBJetCol[topjet1], nonBJetCol[topjet2] );
        top_leptonic.Fill( BJetCol[leptonicTopbjet], isoLep, EvtInfo.PFMET, EvtInfo.PFMETPhi  );

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

        h1.GetTH1("Evt_Top_Hadronic_Chi2")->Fill( chi2               );
        h1.GetTH1("Evt_Top_Hadronic_Mass")->Fill( top_hadronic.Mass  );
        h1.GetTH1("Evt_Top_Hadronic_Pt"  )->Fill( top_hadronic.Pt    );
        h1.GetTH1("Evt_Top_Hadronic_Eta" )->Fill( top_hadronic.Eta   );
        h1.GetTH1("Evt_Top_Hadronic_Phi" )->Fill( top_hadronic.Phi   );
        h1.GetTH1("Evt_Top_Leptonic_Mt"  )->Fill( top_leptonic.MassT );
        h1.GetTH1("Evt_Top_Leptonic_Phi" )->Fill( top_leptonic.Phi   );

    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void bbarTagEff::endJob(){
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
