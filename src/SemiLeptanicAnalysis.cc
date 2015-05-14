// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
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
	Owrt_(iConfig.getParameter<double>("Owrt")),
	Debug_(iConfig.getParameter<bool>("Debug"))
{}

SemiLeptanicAnalysis::~SemiLeptanicAnalysis()
{ 
	delete chain_;
}

// ------------ Other function -------------
template<class TH1>
void SemiLeptanicAnalysis::setCutFlow( TH1* h, std::string channel )
{
	if( channel.compare("lj") == 0 || channel.compare("ljm") == 0 || channel.compare("lje") == 0){
		h->GetXaxis()->SetBinLabel(1,"All");
		h->GetXaxis()->SetBinLabel(2,"#geq1 goodVtx");
		if( channel.compare("lj") == 0  ) h->GetXaxis()->SetBinLabel(3,"#geq1 Lep");
		if( channel.compare("ljm") == 0 ) h->GetXaxis()->SetBinLabel(3,"#geq1 Mu");
		if( channel.compare("lje") == 0 ) h->GetXaxis()->SetBinLabel(3,"#geq1 El");
		if( channel.compare("lj") == 0  ) h->GetXaxis()->SetBinLabel(4,"1 isoLep");
		if( channel.compare("ljm") == 0 ) h->GetXaxis()->SetBinLabel(4,"1 isoMu");
		if( channel.compare("lje") == 0 ) h->GetXaxis()->SetBinLabel(4,"1 isoEl");
		h->GetXaxis()->SetBinLabel(5,"veto(Loose #mu)");
		h->GetXaxis()->SetBinLabel(6,"veto(Loose e)");
		h->GetXaxis()->SetBinLabel(7,"#geq3 Jets");
		h->GetXaxis()->SetBinLabel(8,"#geq2 bjets");
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
	h1.addNewTH1( "Evt2b_O7",	 "O7",	  	"O_{7}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O7Asym",	 "A_{O7}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt2b_O7_Mu",	 "O7",	  	"O_{7}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O7_El",	 "O7",	  	"O_{7}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O7Asym_Mu", "A_{O7}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt2b_O7Asym_El", "A_{O7}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt2b_O2",	 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O2_Mu",	 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O2_El",	 "O2",	  	"O_{2}", "Events", 	"", 	"", 40, -2,   2) ;
	h1.addNewTH1( "Evt2b_O2Asym",	 "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt2b_O2Asym_Mu", "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	h1.addNewTH1( "Evt2b_O2Asym_El", "A_{O2}",	"", 	 "Events", 	"", 	"",  2, 0,   2) ;
	
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
	h1.addNewTH1( "Evt_CutFlow_Mu",    "",         	 		"",     "Evetns", "", "",    9, 0, 9 );
	h1.addNewTH1( "Evt_CutFlow_El",    "",         	  		"",     "Evetns", "", "",    9, 0, 9 );
	h1.addNewTH1( "Evt_MuCut",     	   "isoMu:looseMu:looseEl",    	"",     "Evetns", "", "",    7, 0, 7 );
	h1.addNewTH1( "Evt_ElCut",     	   "isoEl:looseMu:looseEl",     "",     "Evetns", "", "",    7, 0, 7 );
	h1.addNewTH1( "Evt_SameChannel",   "",     "",     "Evetns", "", "",    1, 0, 1 );
	h1.addNewTH1( "Evt2b_SameChannel", "",     "",     "Evetns", "", "",    1, 0, 1 );

	h1.CreateTH1( fs );
	h1.Sumw2();

	setCutFlow( h1.GetTH1("Evt_CutFlow"), "lj" );
	setCutFlow( h1.GetTH1("Evt_CutFlow_Mu"), "ljm" );
	setCutFlow( h1.GetTH1("Evt_CutFlow_El"), "lje" );
	setObservableHist(h1.GetTH1("Evt_O7Asym"),    "O_{7}");
	setObservableHist(h1.GetTH1("Evt_O7Asym_Mu"), "O_{7}");
	setObservableHist(h1.GetTH1("Evt_O7Asym_El"), "O_{7}");
	setObservableHist(h1.GetTH1("Evt_O2Asym"),    "O_{2}");
	setObservableHist(h1.GetTH1("Evt_O2Asym_Mu"), "O_{2}");
	setObservableHist(h1.GetTH1("Evt_O2Asym_El"), "O_{2}");
	setObservableHist(h1.GetTH1("Evt2b_O7Asym"),    "O_{7}");
	setObservableHist(h1.GetTH1("Evt2b_O7Asym_Mu"), "O_{7}");
	setObservableHist(h1.GetTH1("Evt2b_O7Asym_El"), "O_{7}");
	setObservableHist(h1.GetTH1("Evt2b_O2Asym"),    "O_{2}");
	setObservableHist(h1.GetTH1("Evt2b_O2Asym_Mu"), "O_{2}");
	setObservableHist(h1.GetTH1("Evt2b_O2Asym_El"), "O_{2}");
	setLeptonSelHist(h1.GetTH1("Evt_MuCut"));
	setLeptonSelHist(h1.GetTH1("Evt_ElCut"));

	std::cout<<">> [INFO] Chaining "<<inputFiles_.size()<<" files..."<<endl;
	chain_  = new TChain(inputTTree_.c_str());

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		f->Close();
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

		// Vertex selection
		vector<Vertex> selVertex;
		for( int idx=0; idx < VtxInfo.Size; idx++)
		{
			Vertex vtx( VtxInfo, idx );
			if( vtx.isFake ) continue;
			if( vtx.Ndof < 4 ) continue;
			if( abs(vtx.z) > 24 ) continue;
			if( sqrt((vtx.x*vtx.x)+(vtx.y*vtx.y)) > 2 ) continue;
			selVertex.push_back(vtx);
		}

		//* Jet selection
		vector<Jet> seljetCol, bjetCol;
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
			}

			if(  jet.Pt > 30. && jet.CombinedSVBJetTags > 0.679 )
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
		h1.GetTH1("Evt_NJets")->Fill(JetInfo.Size);	
		h1.GetTH1("Evt_NSelJets")->Fill(seljetCol.size());	
		h1.GetTH1("Evt_NbJets")->Fill(bjetCol.size());

		//* Lepton selection
		vector<Lepton> selMuCol, looseMuCol_isoMu, looseElCol_isoMu;
		vector<Lepton> selElCol, looseMuCol_isoEl, looseElCol_isoEl;
		for( int idx=0; idx < LepInfo.Size; idx++)
		{
			Lepton lepton( LepInfo, idx );
			float relIso1=( lepton.ChargedHadronIso + lepton.NeutralHadronIso + lepton.PhotonIso)/fabs(lepton.Pt);
			// Electron selections
			if( lepton.LeptonType == 11 )
			{
				if( abs(lepton.Eta) > 2.5 ) continue;
				if( abs(lepton.Eta) < 1.5 && abs(lepton.Eta) > 1.4442) continue;
				if( lepton.Pt < 15) continue;
				if( relIso1 < 0.2 )
					looseElCol_isoMu.push_back(lepton);	
				if( lepton.Pt < 20) continue;
				if( lepton.Pt < 45 &&
				    relIso1 < 1.
  				  )
					looseElCol_isoEl.push_back(lepton);
				else if( relIso1 < 0.1 &&
					 lepton.Pt > 45	&&
					 lepton.EgammaCutBasedEleIdTRIGGERWP70 
					)
					selElCol.push_back(lepton);	
			}
			// Muon selections
			if( lepton.LeptonType == 13 )
			{
				if( relIso1 > 0.2 ) continue;	
				if( abs(lepton.Eta) > 2.5 ) continue;
				if( lepton.Pt < 10) continue;
				if( (lepton.MuType&0x02) != 0 ) 
					looseMuCol_isoEl.push_back(lepton);	
				if( lepton.Pt < 35)
					looseMuCol_isoMu.push_back(lepton);	
				else if( relIso1 < 0.125 &&
					 (lepton.MuType&0x02) != 0 && // Global muon 
					 (lepton.MuType&0x04) != 0 && // Tracker muon 
					 lepton.Pt >= 35 && 
					 lepton.MuNMuonhits > 0  &&
					 lepton.MuNPixelLayers > 0  &&
					 lepton.MuNTrackerHits > 10  &&
					 lepton.MuInnerTrackNHits > 10  &&
					 lepton.MuNMatchedStations > 1  &&
					 lepton.MuInnerTrackDxy_BS < 0.02 &&
					 lepton.MuGlobalNormalizedChi2 < 10 &&
					 abs(lepton.Eta) < 2.1  
					) // tight muon
					selMuCol.push_back(lepton);
			}
		}

		//* Event selection
		Jet jet1, bjet1, bjet2;
		Lepton isoMu, isoEl;
		bool isMuonChannel(false), isMuonChannel2b(false), isEleChannel(false), isEleChannel2b(false);

		h1.GetTH1("Evt_CutFlow")->Fill("All", 1);
		h1.GetTH1("Evt_CutFlow_El")->Fill("All", 1);
		h1.GetTH1("Evt_CutFlow_Mu")->Fill("All", 1);

		if( selVertex.size() )
		{
			h1.GetTH1("Evt_CutFlow")->Fill("#geq1 goodVtx", 1);
			h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 goodVtx", 1);
			h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 goodVtx", 1);
			//if( ( selMuCol.size() + selElCol.size()) == 1  )
			if( ( selMuCol.size() + selElCol.size()) > 0 )
			{
				h1.GetTH1("Evt_CutFlow")->Fill("#geq1 Lep", 1);
				// Muon channel
				if( selMuCol.size() > 0 )
				{ 
					h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq1 Mu", 1);
					if( selMuCol.size() == 1 )
					{
						isoMu = selMuCol[0];
						h1.GetTH1("Evt_CutFlow")->Fill("1 isoLep", 1);
						h1.GetTH1("Evt_CutFlow_Mu")->Fill("1 isoMu", 1);

						if( looseMuCol_isoMu.size() == 0 )
						{
							h1.GetTH1("Evt_CutFlow")->Fill("veto(Loose #mu)", 1);
							h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose #mu)", 1);

							if( looseElCol_isoMu.size() == 0 )
							{
								h1.GetTH1("Evt_CutFlow")->Fill("veto(Loose e)", 1);
								h1.GetTH1("Evt_CutFlow_Mu")->Fill("veto(Loose e)", 1);
								if( seljetCol.size() >= 3 )
								{ 	
									h1.GetTH1("Evt_CutFlow")->Fill("#geq3 Jets", 1);
									h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq3 Jets", 1);
									// Lable the hardest non_bjet 
									int j1=-1;
									double pt1=0;
									const int size_seljetCol = seljetCol.size();	
									for( int i=0; i < size_seljetCol; i++)
									{
										if( seljetCol[i].CombinedSVBJetTags < 0.679 && pt1 < seljetCol[i].Pt ){
											j1=i;
											pt1=seljetCol[i].Pt;
										}
									}
									jet1 = seljetCol[j1];
								}
								if( bjetCol.size() >= 2 )
								{
									isMuonChannel=true;	
									h1.GetTH1("Evt_CutFlow")->Fill("#geq2 bjets", 1);
									h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq2 bjets", 1);
									//Lable bjet by Pt
									get2HighPtObject( bjetCol, bjet1, bjet2 );
									h1.GetTH1("bJet12_Px")->Fill(bjet1.Px);
									h1.GetTH1("bJet12_Px")->Fill(bjet2.Px);
									h1.GetTH1("bJet12_Py")->Fill(bjet1.Py);
									h1.GetTH1("bJet12_Py")->Fill(bjet2.Py);
									h1.GetTH1("bJet12_Pz")->Fill(bjet1.Pz);
									h1.GetTH1("bJet12_Pz")->Fill(bjet2.Pz);
								}
								if( bjetCol.size() == 2 )
								{	
									isMuonChannel2b=true;	
									h1.GetTH1("Evt_CutFlow")->Fill("=2 bjets", 1);
									h1.GetTH1("Evt_CutFlow_Mu")->Fill("=2 bjets", 1);
								}
							}
						}
					}
				} // [END] Muon channel
				//else if( selElCol.size() == 1 ) // Electron channel	
				if( selElCol.size() > 0 )
				{ 
					h1.GetTH1("Evt_CutFlow_El")->Fill("#geq1 El", 1);
					if( selElCol.size() == 1 ) // Electron channel	
					{
						isoEl=selElCol[0];
						h1.GetTH1("Evt_CutFlow")->Fill("1 isoLep", 1);
						h1.GetTH1("Evt_CutFlow_El")->Fill("1 isoEl", 1);

						if( looseMuCol_isoEl.size() == 0 )
						{
							h1.GetTH1("Evt_CutFlow")->Fill("veto(Loose #mu)", 1);
							h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose #mu)", 1);

							if( looseElCol_isoEl.size() == 0 )
							{
								h1.GetTH1("Evt_CutFlow")->Fill("veto(Loose e)", 1);
								h1.GetTH1("Evt_CutFlow_El")->Fill("veto(Loose e)", 1);

								if( seljetCol.size() >= 3 )
								{ 	
									h1.GetTH1("Evt_CutFlow")->Fill("#geq3 Jets", 1);
									h1.GetTH1("Evt_CutFlow_El")->Fill("#geq3 Jets", 1);
									// Lable the hardest non_bjet 
									int j1=-1;
									double pt1=0;	
									const int size_seljetCol = seljetCol.size();	
									for( int i=0; i<size_seljetCol; i++)
									{
										if( seljetCol[i].CombinedSVBJetTags < 0.679 && pt1 < seljetCol[i].Pt ){
											j1=i;
											pt1=seljetCol[i].Pt;
										}
									}
									jet1=seljetCol[j1];

									if( bjetCol.size() >= 2 )
									{ 
										isEleChannel=true;
										h1.GetTH1("Evt_CutFlow")->Fill("#geq2 bjets", 1);
										h1.GetTH1("Evt_CutFlow_El")->Fill("#geq2 bjets", 1);
										// Lable bjet by Pt
										get2HighPtObject( bjetCol, bjet1, bjet2 );
										h1.GetTH1("bJet12_Px")->Fill(bjet1.Px);
										h1.GetTH1("bJet12_Px")->Fill(bjet2.Px);
										h1.GetTH1("bJet12_Py")->Fill(bjet1.Py);
										h1.GetTH1("bJet12_Py")->Fill(bjet2.Py);
										h1.GetTH1("bJet12_Pz")->Fill(bjet1.Pz);
										h1.GetTH1("bJet12_Pz")->Fill(bjet2.Pz);
									}
									if( bjetCol.size() == 2 )
									{
										isEleChannel2b=true;	
										h1.GetTH1("Evt_CutFlow")->Fill("=2 bjets", 1);
										h1.GetTH1("Evt_CutFlow_El")->Fill("=2 bjets", 1);
									}
								}
							}	
						}
					}				
				} // [END] Electron channel
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
			} // [END] one lepton selection
		} // [END] Vxt selection

		//* Fill observables O7 and O2
		//* bJets >= 2
		Lepton isoLep;
		if( isMuonChannel && !isEleChannel ) isoLep = isoMu;
		else if( !isMuonChannel && isEleChannel ) isoLep = isoEl;

		TVector3 O2_1v =  bjet1.P3 + bjet2.P3;
		TVector3 O2_2v = isoLep.P3.Cross( jet1.P3 );
		double O2 = O2_1v.Dot( O2_2v );

		double O7_1z = az.Dot( bjet1.P3 - bjet2.P3 );
		double O7_2z = az.Dot( bjet1.P3.Cross( bjet2.P3 ));
		double O7 = O7_1z * O7_2z;

		//if( isMuonChannel && !isEleChannel )
		if( isMuonChannel   && isEleChannel   ) h1.GetTH1("Evt_SameChannel")->Fill(0);
		if( isMuonChannel2b && isEleChannel2b ) h1.GetTH1("Evt2b_SameChannel")->Fill(0);
		if( isMuonChannel )
		{
			h1.GetTH1("Evt_O2")->Fill( O2/Owrt_ );	
			h1.GetTH1("Evt_O2_Mu")->Fill( O2/Owrt_ );
			if( O2 > 0 ){
				h1.GetTH1("Evt_O2Asym")->Fill("O_{2}>0",1);
				h1.GetTH1("Evt_O2Asym_Mu")->Fill("O_{2}>0",1);
			}else{
				h1.GetTH1("Evt_O2Asym")->Fill("O_{2}<0",1);
				h1.GetTH1("Evt_O2Asym_Mu")->Fill("O_{2}<0",1);
			}

			h1.GetTH1("Evt_O7")->Fill( O7/Owrt_ );	
			h1.GetTH1("Evt_O7_Mu")->Fill( O7/Owrt_ );
			h1.GetTH1("Evt_O7_term1")->Fill( O7_1z );
			h1.GetTH1("Evt_O7_term2")->Fill( O7_2z );
			if( O7 > 0 ){
				h1.GetTH1("Evt_O7Asym")->Fill("O_{7}>0",1);
				h1.GetTH1("Evt_O7Asym_Mu")->Fill("O_{7}>0",1);
			}else{
				h1.GetTH1("Evt_O7Asym")->Fill("O_{7}<0",1);
				h1.GetTH1("Evt_O7Asym_Mu")->Fill("O_{7}<0",1);
			}
		}
		//else if( !isMuonChannel && isEleChannel )
		if( isEleChannel )
		{
			h1.GetTH1("Evt_O2")->Fill(O2/Owrt_);	
			h1.GetTH1("Evt_O2_El")->Fill(O2/Owrt_);
			if( O2 > 0 ){
				h1.GetTH1("Evt_O2Asym")->Fill("O_{2}>0",1);
				h1.GetTH1("Evt_O2Asym_El")->Fill("O_{2}>0",1);
			}else{
				h1.GetTH1("Evt_O2Asym")->Fill("O_{2}<0",1);
				h1.GetTH1("Evt_O2Asym_El")->Fill("O_{2}<0",1);
			}

			h1.GetTH1("Evt_O7")->Fill(O7/Owrt_);	
			h1.GetTH1("Evt_O7_El")->Fill(O7/Owrt_);	
			h1.GetTH1("Evt_O7_term1")->Fill(O7_1z);
			h1.GetTH1("Evt_O7_term2")->Fill(O7_2z);
			if( O7 > 0 ){
				h1.GetTH1("Evt_O7Asym")->Fill("O_{7}>0", 1);
				h1.GetTH1("Evt_O7Asym_El")->Fill("O_{7}>0",1);
			}else{
				h1.GetTH1("Evt_O7Asym")->Fill("O_{7}<0",1);
				h1.GetTH1("Evt_O7Asym_El")->Fill("O_{7}<0",1);
			}
		}
		//* bJets == 2
		//if( isMuonChannel2b && !isEleChannel2b )
		if( isMuonChannel2b )
		{
			h1.GetTH1("Evt2b_O2")->Fill(O2/Owrt_);	
			h1.GetTH1("Evt2b_O2_Mu")->Fill(O2/Owrt_);
			if( O2 > 0 ){
				h1.GetTH1("Evt2b_O2Asym")->Fill("O_{2}>0",1);
				h1.GetTH1("Evt2b_O2Asym_Mu")->Fill("O_{2}>0",1);
			}else{
				h1.GetTH1("Evt2b_O2Asym")->Fill("O_{2}<0",1);
				h1.GetTH1("Evt2b_O2Asym_Mu")->Fill("O_{2}<0",1);
			}

			h1.GetTH1("Evt2b_O7")->Fill(O7/Owrt_);	
			h1.GetTH1("Evt2b_O7_Mu")->Fill(O7/Owrt_);
			if( O7 > 0 ){
				h1.GetTH1("Evt2b_O7Asym")->Fill("O_{7}>0",1);
				h1.GetTH1("Evt2b_O7Asym_Mu")->Fill("O_{7}>0",1);
			}else{
				h1.GetTH1("Evt2b_O7Asym")->Fill("O_{7}<0",1);
				h1.GetTH1("Evt2b_O7Asym_Mu")->Fill("O_{7}<0",1);
			}
		}
		//else if( !isMuonChannel2b && isEleChannel2b )
		if( isEleChannel2b )
		{
			h1.GetTH1("Evt2b_O2")->Fill(O2/Owrt_);	
			h1.GetTH1("Evt2b_O2_El")->Fill(O2/Owrt_);
			if( O2 > 0 ){
				h1.GetTH1("Evt2b_O2Asym")->Fill("O_{2}>0",1);
				h1.GetTH1("Evt2b_O2Asym_El")->Fill("O_{2}>0",1);
			}else{
				h1.GetTH1("Evt2b_O2Asym")->Fill("O_{2}<0",1);
				h1.GetTH1("Evt2b_O2Asym_El")->Fill("O_{2}<0",1);
			}

			h1.GetTH1("Evt2b_O7")->Fill(O7/Owrt_);	
			h1.GetTH1("Evt2b_O7_El")->Fill(O7/Owrt_);	
			if( O7 > 0 ){
				h1.GetTH1("Evt2b_O7Asym")->Fill("O_{7}>0", 1);
				h1.GetTH1("Evt2b_O7Asym_El")->Fill("O_{7}>0",1);
			}else{
				h1.GetTH1("Evt2b_O7Asym")->Fill("O_{7}<0",1);
				h1.GetTH1("Evt2b_O7Asym_El")->Fill("O_{7}<0",1);
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
