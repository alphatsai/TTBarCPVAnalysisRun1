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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/checkEvtTool.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h" 

class LeptonJetFilter: public edm::EDAnalyzer
{
    public:
        explicit LeptonJetFilter(const edm::ParameterSet&);
        ~LeptonJetFilter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        std::string int2str(int i); 

        int maxEvents_;
        const int reportEvery_;
        std::vector<std::string> inputFiles_;
        std::vector<std::string> inputJsons_;
        const std::string inputTTree_; 
        std::vector<int> MuonHLT_; 
        std::vector<int> ElectronHLT_; 
        const double LepPt_;
        const double JetPt_;
        const int NJet_;
        const int NLep_;
        const bool Debug_;

        edm::Service<TFileService> fs;
        TChain* chain_;
        TTree* newtree;
        TH1D* h_CutFlow;
        TH1D* h_Events;

        EvtInfoBranches    EvtInfo;
        JetInfoBranches    JetInfo;
        LepInfoBranches    LepInfo;

};
//
// constructors and destructor
//
LeptonJetFilter::LeptonJetFilter(const edm::ParameterSet& iConfig) : 
    maxEvents_(     iConfig.getParameter<int>("MaxEvents")), 
    reportEvery_(   iConfig.getParameter<int>("ReportEvery")),
    inputFiles_(    iConfig.getParameter<std::vector<std::string> >("InputFiles")),
    inputJsons_(    iConfig.getParameter<std::vector<std::string> >("InputJsons")),
    inputTTree_(    iConfig.getParameter<std::string>("InputTTree")),
    MuonHLT_(       iConfig.getParameter<std::vector<int>>("MuonHLT")),
    ElectronHLT_(   iConfig.getParameter<std::vector<int>>("ElectronHLT")),
    LepPt_(         iConfig.getParameter<double>("LepPt")),
    JetPt_(         iConfig.getParameter<double>("JetPt")),
    NJet_(          iConfig.getParameter<int>("NJet")),
    NLep_(          iConfig.getParameter<int>("NLep")),
    Debug_(         iConfig.getParameter<bool>("Debug"))
{}

LeptonJetFilter::~LeptonJetFilter()
{ 
    delete chain_;
}

// ------------ Other function -------------
std::string LeptonJetFilter::int2str( int i )
{
    std::string s;
    std::stringstream ss(s);
    ss << i;
    return ss.str(); 
}
// ------------ method called once each job just before starting event loop  ------------
void LeptonJetFilter::beginJob()
{ 
    h_CutFlow = fs->make<TH1D>("Evt_CutFlow", "", 4, 0, 4);
    h_Events  = fs->make<TH1D>("Evt_Events",  "", 1, 0, 1);

    h_CutFlow->GetXaxis()->SetBinLabel(1, "Original");
    h_CutFlow->GetXaxis()->SetBinLabel(2, "JSON");
    h_CutFlow->GetXaxis()->SetBinLabel(3, "HLT");
    h_CutFlow->GetXaxis()->SetBinLabel(4, ("#geq "+int2str(NLep_)+" Leptons").c_str());
    h_CutFlow->GetXaxis()->SetBinLabel(5, ("#geq "+int2str(NJet_)+" Jets").c_str());

    std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<std::endl;
    chain_  = new TChain(inputTTree_.c_str());

    for(unsigned i=0; i<inputFiles_.size(); ++i){
        chain_->Add(inputFiles_.at(i).c_str());
        TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
        f->Close();
    }
    //newtree = fs->make<TTree>("tree", "") ;
    fs->cd();
    newtree = chain_->CloneTree(0);

    EvtInfo.Register(chain_);
    JetInfo.Register(chain_,"PFJetInfo");
    LepInfo.Register(chain_,"PFLepInfo");

    if( maxEvents_<0 || maxEvents_>chain_->GetEntries() ) maxEvents_ = chain_->GetEntries();

    return;  
}

// ------------ method called for each event  ------------
void LeptonJetFilter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    for(int entry=0; entry<maxEvents_; entry++)
    {
        chain_->GetEntry(entry);

        if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

        h_CutFlow->Fill("Original", 1);
        h_Events->Fill(0);

        // Check evt by json 
        bool isdata(0);
        isdata  = EvtInfo.McFlag ? 0 : 1; 
        if( isdata && !checkEvt.isGoodEvt( EvtInfo.RunNo, EvtInfo.LumiNo ) ) continue;
        h_CutFlow->Fill("JSON", 1);

        // HLT selection
        bool passHLT=false;
        for( std::vector<int>::const_iterator ihlt = ElectronHLT_.begin(); ihlt != ElectronHLT_.end(); ++ihlt )
        { 
            if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
                passHLT=true;
                break;
            }
        }
        if( !passHLT )
        { 
            for( std::vector<int>::const_iterator ihlt = MuonHLT_.begin(); ihlt != MuonHLT_.end(); ++ihlt )
            {
                if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
                    passHLT=true;
                    break;
                }
            }
        }

        if( !passHLT ) continue;
        h_CutFlow->Fill("HLT", 1);

        //* Lepton selection
        int numPassLep(0);
        for( int idx=0; idx < LepInfo.Size; idx++)
        {
            if( LepInfo.Pt[idx] > LepPt_ ) numPassLep++;
        }
        if( numPassLep < NLep_ ) continue;
        h_CutFlow->Fill(("#geq "+int2str(NLep_)+" Leptons").c_str(), 1);

        //* Jet selection
        int numPassJets(0);
        for( int idx=0; idx < JetInfo.Size; idx++)
        {
            if( JetInfo.Pt[idx] > JetPt_ ) numPassJets++;
        }
        if( numPassJets < NJet_ ) continue;
        h_CutFlow->Fill(("#geq "+int2str(NJet_)+" Jets").c_str(), 1);

        newtree->Fill();

    }//// [END] entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void LeptonJetFilter::endJob(){
    std::cout<<">> [INFO] End of Job!"<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void LeptonJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonJetFilter);
