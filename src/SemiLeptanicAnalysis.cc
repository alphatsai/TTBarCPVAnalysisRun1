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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1InfoClass.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/SemiLeptanicAnalysis.h" 

//
// constructors and destructor
//
SemiLeptanicAnalysis::SemiLeptanicAnalysis(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	reportEvery_(iConfig.getParameter<int>("ReportEvery")),
	inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
	inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
	MuonHLT_(iConfig.getParameter<std::vector<int>>("MuonHLT")),
	ElectronHLT_(iConfig.getParameter<std::vector<int>>("ElectronHLT")),
	NJets_(iConfig.getParameter<double>("NJets")),
	IsoEleEt_(iConfig.getParameter<double>("IsoEleEt")),
	IsoMuonPt_(iConfig.getParameter<double>("IsoMuonPt")),
	NonBjetCSVThr_(iConfig.getParameter<double>("NonBjetCSVThr")),
	Owrt_(iConfig.getParameter<double>("Owrt")),
	Debug_(iConfig.getParameter<bool>("Debug")),
	isSkim_(false)
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
	h1.addNewTH1( "Evt_O7_Mu",	 "O7",	  	"O_{7}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt_O7_El",	 "O7",	  	"O_{7}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt_O7Asym_Mu",	 "A_{O7}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt_O7Asym_El",	 "A_{O7}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt_O2",		 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt_O2_Mu",	 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt_O2_El",	 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt_O2Asym",	 "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt_O2Asym_Mu",	 "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt_O2Asym_El",	 "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	
	h1.addNewTH1( "Evt_isoMu_Pt",   "pT of isoMuon",	  "p_{T}(#mu)", "Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoMu_Et",   "ET of isoMuon",	  "E_{T}(#mu)", "Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoMu_E",	"Energy of isoMuon",	  "Energy(#mu)","Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoMu_Eta",	"Eta of isoMuon",	  "#eta(#mu)", 	"Yields", 	"", 	"",	100, -5, 5 );
	h1.addNewTH1( "Evt_isoMu_Phi",	"Phi of isoMuon",	  "#phi(#mu)", 	"Yields", 	"", 	"",	64, -3.2,   3.2 ) ;

	h1.addNewTH1( "Evt_isoEl_Pt",   "pT of isEle",	  "p_{T}(e)", "Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoEl_Et",   "ET of isEle",	  "E_{T}(e)", "Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoEl_E",	"Energy of isEle","Energy(e)","Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_isoEl_Eta",	"Eta of isEle",	  "#eta(e)",  "Yields", 	"", 	"",	100, -5, 5 );
	h1.addNewTH1( "Evt_isoEl_Phi",	"Phi of isEle",	  "#phi(e)",  "Yields", 	"", 	"",	64, -3.2,   3.2 ) ;

	h1.addNewTH1( "Evt_HardJet_Pt",    "pT of HardJet",	 "p_{T}(HardJet)", 	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_HardJet_M",	   "Mass of HardJet",    "Mass(HardJet)", 	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_HardJet_E",	   "Energy of HardJet",  "Energy(HardJet)",	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_HardJet_Eta",   "Eta of HardJet",	 "#eta(HardJet)", 	"Yields", 	"", 	"",	100, -5, 5 );
	h1.addNewTH1( "Evt_HardJet_Phi",   "Phi of HardJet",	 "#phi(HardJet)", 	"Yields", 	"", 	"",	64, -3.2,   3.2 ) ;
	h1.addNewTH1( "Evt_HardJet_BTag",  "HardJet b-tagged",   "bTag", 		"Yields", 	"", 	"",	100, 0,   1 ) ;

	h1.addNewTH1( "Evt_bJet1_Pt",    "pT of b-Jet",	   "p_{T}(B-tagged j)",	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet1_M",	 "Mass of b-Jet",  "Mass(B-tagged j)", 	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet1_E",	 "Energy of b-Jet","Energy(B-tagged j)","Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet1_Eta",	 "Eta of b-Jet",   "#eta(B-tagged j)", 	"Yields", 	"", 	"",	100, -5, 5 );
	h1.addNewTH1( "Evt_bJet1_Phi",	 "Phi of b-Jet",   "#phi(B-tagged j)", 	"Yields", 	"", 	"",	64, -3.2,   3.2 ) ;
	h1.addNewTH1( "Evt_bJet1_BTag",	 "b-Jet b-tagged", "bTag", 		"Yields", 	"", 	"",	100, 0,   1 ) ;

	h1.addNewTH1( "Evt_bJet2_Pt",    "pT of b-Jet",	   "p_{T}(B-tagged j)",	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet2_M",	 "Mass of b-Jet",  "Mass(B-tagged j)", 	"Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet2_E",	 "Energy of b-Jet","Energy(B-tagged j)","Yields", 	"GeV", 	"",	500, 0,   500 ) ;
	h1.addNewTH1( "Evt_bJet2_Eta",	 "Eta of b-Jet",   "#eta(B-tagged j)", 	"Yields", 	"", 	"",	100, -5, 5 );
	h1.addNewTH1( "Evt_bJet2_Phi",	 "Phi of b-Jet",   "#phi(B-tagged j)", 	"Yields", 	"", 	"",	64, -3.2,   3.2 ) ;
	h1.addNewTH1( "Evt_bJet2_BTag",	 "b-Jet b-tagged", "bTag", 		"Yields", 	"", 	"",	100, 0,   1 ) ;

	h1.addNewTH1( "Evt_NLeptons",	   "Num. of leptons",	  	"N(lep)", 	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NSelLeptons",   "Num. of selected leptons", 	"N(selected lep)","Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NMuons",	   "Num. of muon",	     	"N(#mu)", 	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NSelMuons",	   "Num. of selected muon",    	"N(selected #mu)","Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NLooseMuIsoMu", "Num. of loose muon", 	"N(loose #mu)",	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NLooseElIsoMu", "Num. of loose electron",   	"N(loose e)",	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NElectrons",	   "Num. of electron",	  	"N(e)", 	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NSelElectrons", "Num. of selected electron",	"N(selected e)",  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NLooseMuIsoEl", "Num. of loose muon", 	"N(loose #mu)",	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_NLooseElIsoEl", "Num. of loose electron", 	"N(loose e)",	  "Events", 	"", 	"",		10, 0,   10) ;
	h1.addNewTH1( "Evt_CutFlow_Mu",    "",         	 		"",     "Evetns", "", "",    10, 0, 10 );
	h1.addNewTH1( "Evt_CutFlow_El",    "",         	  		"",     "Evetns", "", "",    10, 0, 10 );
	h1.addNewTH1( "Evt_MuCut",     	   "isoMu:looseMu:looseEl",    	"",     "Evetns", "", "",    7, 0, 7 );
	h1.addNewTH1( "Evt_ElCut",     	   "isoEl:looseMu:looseEl",     "",     "Evetns", "", "",    7, 0, 7 );
	h1.addNewTH1( "Evt_SameChannel",   "",     "",     "Evetns", "", "",    1, 0, 1 );
	h1.addNewTH1( "Evt_Triger", "",     "",     "", "", "",    5900, 0, 5900 );

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
		for( std::vector<int>::const_iterator ihlt = ElectronHLT_.begin(); ihlt != ElectronHLT_.end(); ++ihlt )
		{	
			if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
				passElectronHLT=true;
				break;
			}
		}	
		for( std::vector<int>::const_iterator ihlt = MuonHLT_.begin(); ihlt != MuonHLT_.end(); ++ihlt )
		{
			if( int(EvtInfo.TrgBook[*ihlt]) == 1 ){
				passMuonHLT=true;
				break;
			}
		}

		// Vertex selection
		vector<Vertex> selVertex;
		for( int idx=0; idx < VtxInfo.Size; idx++)
		{
			Vertex vtx( VtxInfo, idx );
			if( vtx.Type != 1 ) continue;
			if( vtx.isFake ) continue;
			if( vtx.Ndof < 4 ) continue;
			if( abs(vtx.z) > 24 ) continue;
			if( vtx.Rho > 2 ) continue;
			selVertex.push_back(vtx);
		}

		//* Jet selection
		vector<Jet> seljetCol, bjetCol, nonbjetCol;
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
			
			// ID tight
			if( jet.CHF == 0 ) continue;
			if( jet.NEF > 0.99 ) continue;
			if( jet.NHF > 0.99 ) continue;
			if( jet.CEF > 0.99 ) continue;
			if( jet.NCH == 0 ) continue;
			if( jet.NConstituents <= 1 ) continue;
			if( abs( jet.Eta ) >= 2.4 ) continue;
			if( jet.Pt > 30. )
			{ 
				seljetCol.push_back(jet);
				h1.GetTH1("SelJet_Pt")->Fill( jet.Pt );			
				h1.GetTH1("SelJet_Px")->Fill( jet.Px );			
				h1.GetTH1("SelJet_Py")->Fill( jet.Py );			
				h1.GetTH1("SelJet_Pz")->Fill( jet.Pz );			
				h1.GetTH1("SelJet_M")->Fill(  jet.Mass );			
				h1.GetTH1("SelJet_E")->Fill(  jet.Energy );			
				h1.GetTH1("SelJet_Eta")->Fill(  jet.Eta );			
				h1.GetTH1("SelJet_Phi")->Fill(  jet.Phi );			
				h1.GetTH1("SelJet_BTag")->Fill( jet.CombinedSVBJetTags );	
	
				if( jet.CombinedSVBJetTags < NonBjetCSVThr_ ) nonbjetCol.push_back(jet);	
				if( jet.CombinedSVBJetTags > 0.679 )
				{ 
					bjetCol.push_back(jet);
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
			}
		}
		h1.GetTH1("Evt_NJets")->Fill(JetInfo.Size);	
		h1.GetTH1("Evt_NSelJets")->Fill(seljetCol.size());	
		h1.GetTH1("Evt_NbJets")->Fill(bjetCol.size());

		//* Lepton selection
		vector<Lepton> selMuCol, looseMuCol_isoMu, looseElCol_isoMu;
		vector<Lepton> selElCol, looseMuCol_isoEl, looseElCol_isoEl;
		for( int idx=0; idx < LepInfo.Size; idx++)
		{
			Lepton lepton( LepInfo, idx );
			// Electron selections
			if( lepton.LeptonType == 11 )
			{
				float relIso1=( lepton.ChargedHadronIsoR03 + lepton.NeutralHadronIsoR03 + lepton.PhotonIsoR03)/fabs(lepton.Pt);
				if( abs(lepton.Eta) > 2.5 ) continue;
				if( abs(lepton.Eta) < 1.566 && abs(lepton.Eta) > 1.4442) continue;
				if( lepton.Et < 20 ) continue;
				if( relIso1 > 0.15 ) continue;
	
				if( lepton.EgammaMVATrig > 0 && lepton.EgammaMVATrig < 1. )
				{
					looseElCol_isoMu.push_back(lepton);	
					if( lepton.Et < IsoEleEt_ )
					{
						looseElCol_isoEl.push_back(lepton);
					}
				}
				if( relIso1 < 0.1 &&
				    lepton.Et > IsoEleEt_ &&
				    lepton.EgammaMVATrig > 0.5 && 
				    lepton.NumberOfExpectedInnerHits <= 0 && 
				    abs(lepton.Eta) < 2.1 && 
				    abs(lepton.ElTrackDxy_PV)<0.02 
				  )
				{
					selElCol.push_back(lepton);
				}	
			}
			// Muon selections
			if( lepton.LeptonType == 13 )
			{
				float relIso1=( lepton.ChargedHadronIsoR04 + lepton.NeutralHadronIsoR04 + lepton.PhotonIsoR04)/fabs(lepton.Pt);
				if( relIso1 > 0.2 ) continue;	
				if( abs(lepton.Eta) > 2.5 ) continue;
				if( lepton.Pt < 10) continue;
				if( (lepton.MuType&0x02) == 0 && (lepton.MuType&0x04) == 0 ) continue; 

				looseMuCol_isoEl.push_back(lepton);	

				if( lepton.Pt < IsoMuonPt_ )
				{
					looseMuCol_isoMu.push_back(lepton);	
				}
				else if( relIso1 < 0.12 &&
					 lepton.Pt >= IsoMuonPt_ && 
					 lepton.MuNMuonhits > 0  &&
					 lepton.MuNMatchedStations > 1  &&
					 lepton.MuGlobalNormalizedChi2 < 10 &&
					 lepton.MuNTrackLayersWMeasurement > 5  &&
					 (lepton.MuType&0x02) != 0 && // Global muon 
					 abs(lepton.MuInnerTrackDxy_PV) < 0.2 &&
					 abs(lepton.Eta) < 2.1  
					 //lepton.MuNPixelLayers > 0  &&
					 //lepton.MuNTrackerHits > 10  &&
					 //lepton.MuInnerTrackNHits > 10  &&
					 //abs(lepton.MuInnerTrackDxy_BS) < 0.02 &&
					) // tight muon
				{
					selMuCol.push_back(lepton);
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

		if( selVertex.size() )
		{
			h1.GetTH1("Evt_CutFlow")->Fill("#geq1 goodVtx", 1);
			h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 goodVtx", 1);
			h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 goodVtx", 1);

			bool passElectronSel(false), passMuonSel(false);
			// Muon channel
			if( passMuonHLT )
			{
				h1.GetTH1("Evt_CutFlow_Mu")->Fill("HLT", 1);
				if( selMuCol.size() > 0 )
				{ 
					h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 Lep", 1);
					if( selMuCol.size() == 1 )
					{
						isoMu = selMuCol[0];
						h1.GetTH1("Evt_CutFlow_Mu")->Fill("1 isoLep", 1);

						if( looseMuCol_isoMu.size() == 0 )
						{
							h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose #mu)", 1);

							if( looseElCol_isoMu.size() == 0 )
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
				if( selElCol.size() > 0 )
				{ 
					h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 Lep", 1);
					if( selElCol.size() == 1 ) 
					{
						isoEl=selElCol[0];
						h1.GetTH1("Evt_CutFlow_El")->Fill("1 isoLep", 1);

						if( looseMuCol_isoEl.size() == 0 )
						{
							h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose #mu)", 1);

							if( looseElCol_isoEl.size() == 0 )
							{
								h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose e)", 1);
								passElectronSel=true;
							}	
						}
					}				
				}
			} // [END] Electron channel

			if( passElectronSel || passMuonSel )
			{
				if( seljetCol.size() >= NJets_ )
				{ 	
					if( passMuonSel )     h1.GetTH1("Evt_CutFlow_Mu")->Fill(("#geq"+int2str(NJets_)+" Jets").c_str(), 1);
					if( passElectronSel ) h1.GetTH1("Evt_CutFlow_El")->Fill(("#geq"+int2str(NJets_)+" Jets").c_str(), 1);

					if( bjetCol.size() == 2 )
					{
						// Lable the hardest non_bjet 
						int j1=-1;
						double pt1=0;
						const int size_seljetCol = nonbjetCol.size();	
						for( int i=0; i < size_seljetCol; i++)
						{
							if( pt1 < nonbjetCol[i].Pt ){
								j1=i;
								pt1=nonbjetCol[i].Pt;
							}
						}
						hardJet = nonbjetCol[j1];
						if( j1 == -1 ){ std::cout<<">>[WARNING] "<<entry<<"Doesn't find hard jet!"<<endl; }

						//Lable bjet by Pt
						get2HighPtObject( bjetCol, bjet1, bjet2 );
						h1.GetTH1("bJet12_Px")->Fill(bjet1.Px);
						h1.GetTH1("bJet12_Px")->Fill(bjet2.Px);
						h1.GetTH1("bJet12_Py")->Fill(bjet1.Py);
						h1.GetTH1("bJet12_Py")->Fill(bjet2.Py);
						h1.GetTH1("bJet12_Pz")->Fill(bjet1.Pz);
						h1.GetTH1("bJet12_Pz")->Fill(bjet2.Pz);

						if( passMuonSel )
						{     
							isGoodMuonEvt=true;	
							h1.GetTH1("Evt_CutFlow_Mu")->Fill("=2 bjets", 1);
						}
						if( passElectronSel )
						{ 
							isGoodElectronEvt=true;
							h1.GetTH1("Evt_CutFlow_El")->Fill("=2 bjets", 1);
						}
					}
				}
			}//[END] Jet and bjet cutflow

			// Check events
			if( passElectronHLT || passMuonHLT ){ 
				h1.GetTH1("Evt_CutFlow")->Fill("HLT", 1);
				if( ( selMuCol.size() + selElCol.size()) > 0 ) h1.GetTH1("Evt_CutFlow")->Fill("#geq1 Lep", 1);
				if( ( selMuCol.size() + selElCol.size()) == 1 ){ 
					h1.GetTH1("Evt_CutFlow")->Fill("1 isoLep", 1);
					if( selMuCol.size() == 1 && looseMuCol_isoMu.size() == 0 && looseElCol_isoMu.size() == 0 )
						h1.GetTH1("Evt_MuCut")->Fill("1:0:0", 1);
					if( selMuCol.size() == 0 && looseMuCol_isoMu.size() == 1 && looseElCol_isoMu.size() == 0 )
						h1.GetTH1("Evt_MuCut")->Fill("0:1:0", 1);
					if( selMuCol.size() == 0 && looseMuCol_isoMu.size() == 0 && looseElCol_isoMu.size() == 1 )
						h1.GetTH1("Evt_MuCut")->Fill("0:0:1", 1);
					if( selMuCol.size() == 1 && looseMuCol_isoMu.size() == 1 && looseElCol_isoMu.size() == 0 )
						h1.GetTH1("Evt_MuCut")->Fill("1:1:0", 1);
					if( selMuCol.size() == 1 && looseMuCol_isoMu.size() == 0 && looseElCol_isoMu.size() == 1 )
						h1.GetTH1("Evt_MuCut")->Fill("1:0:1", 1);
					if( selMuCol.size() == 0 && looseMuCol_isoMu.size() == 1 && looseElCol_isoMu.size() == 1 )
						h1.GetTH1("Evt_MuCut")->Fill("0:1:1", 1);
					if( selMuCol.size() == 1 && looseMuCol_isoMu.size() == 1 && looseElCol_isoMu.size() == 1 )
						h1.GetTH1("Evt_MuCut")->Fill("1:1:1", 1);
					if( selElCol.size() == 1 && looseMuCol_isoEl.size() == 0 && looseElCol_isoEl.size() == 0 )
						h1.GetTH1("Evt_ElCut")->Fill("1:0:0", 1);
					if( selElCol.size() == 0 && looseMuCol_isoEl.size() == 1 && looseElCol_isoEl.size() == 0 )
						h1.GetTH1("Evt_ElCut")->Fill("0:1:0", 1);
					if( selElCol.size() == 0 && looseMuCol_isoEl.size() == 0 && looseElCol_isoEl.size() == 1 )
						h1.GetTH1("Evt_ElCut")->Fill("0:0:1", 1);
					if( selElCol.size() == 1 && looseMuCol_isoEl.size() == 1 && looseElCol_isoEl.size() == 0 )
						h1.GetTH1("Evt_ElCut")->Fill("1:1:0", 1);
					if( selElCol.size() == 1 && looseMuCol_isoEl.size() == 0 && looseElCol_isoEl.size() == 1 )
						h1.GetTH1("Evt_ElCut")->Fill("1:0:1", 1);
					if( selElCol.size() == 0 && looseMuCol_isoEl.size() == 1 && looseElCol_isoEl.size() == 1 )
						h1.GetTH1("Evt_ElCut")->Fill("0:1:1", 1);
					if( selElCol.size() == 1 && looseMuCol_isoEl.size() == 1 && looseElCol_isoEl.size() == 1 )
						h1.GetTH1("Evt_ElCut")->Fill("1:1:1", 1);
				}
			}
		} // [END] Vxt selection

		//* Fill other events plots 
		if( isGoodMuonEvt && isGoodElectronEvt ) h1.GetTH1("Evt_SameChannel")->Fill(0);
		if( isGoodMuonEvt || isGoodElectronEvt )
		{
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
		}
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
