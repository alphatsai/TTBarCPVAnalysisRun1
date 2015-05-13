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
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Vertex.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
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
	Debug_(iConfig.getParameter<bool>("Debug"))
{}

SemiLeptanicAnalysis::~SemiLeptanicAnalysis()
{ 
	delete chain_;
}

// ------------ Other function -------------
template<class TH1>
void SemiLeptanicAnalysis::setCutFlow( TH1* h, std::string channel ){
	if( channel.compare("lj") == 0 || channel.compare("ljm") == 0 || channel.compare("lje") == 0){
		h->GetXaxis()->SetBinLabel(1,"All");
		h->GetXaxis()->SetBinLabel(2,"#geq1 goodVtx");
		if( channel.compare("lj") == 0  ) h->GetXaxis()->SetBinLabel(3,"1 isoLep");
		if( channel.compare("ljm") == 0 ) h->GetXaxis()->SetBinLabel(3,"1 isoMu");
		if( channel.compare("lje") == 0 ) h->GetXaxis()->SetBinLabel(3,"1 isoEl");
		h->GetXaxis()->SetBinLabel(4,"veto(Loose #mu)");
		h->GetXaxis()->SetBinLabel(5,"veto(Loose e)");
		h->GetXaxis()->SetBinLabel(6,"#geq3 Jets");
		h->GetXaxis()->SetBinLabel(7,"#geq2 bjets");
		h->GetXaxis()->SetBinLabel(8,"=2 bjets");
	}
	return ;
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
	h1.addNewTH1( "Evt_CutFlow_Mu",    "",         	 		"",     "Evetns", "", "",    8, 0, 8 );
	h1.addNewTH1( "Evt_CutFlow_El",    "",         	  		"",     "Evetns", "", "",    8, 0, 8 );
	h1.addNewTH1( "Evt_MuCut",     	   "isoMu:looseMu:looseEl",    	"",     "Evetns", "", "",    7, 0, 7 );
	h1.addNewTH1( "Evt_ElCut",     	   "isoEl:looseMu:looseEl",     "",     "Evetns", "", "",    7, 0, 7 );

	h1.CreateTH1( fs );
	h1.Sumw2();

	setCutFlow( h1.GetTH1("Evt_CutFlow"), "lj" );
	setCutFlow( h1.GetTH1("Evt_CutFlow_Mu"), "ljm" );
	setCutFlow( h1.GetTH1("Evt_CutFlow_El"), "lje" );

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

	vector<double> ax, ay, az;
	ax.push_back(1); ax.push_back(0); ax.push_back(0);
	ay.push_back(0); ay.push_back(1); ay.push_back(0);
	az.push_back(0); az.push_back(0); az.push_back(1);

	for(int entry=0; entry<maxEvents_; entry++)
	{
		chain_->GetEntry(entry);

		if( Debug_ ){ if( ( entry%reportEvery_) == 0 ) cout<<">> [DEBUG] "<<entry<<" of "<< maxEvents_<<endl; }

		h1.GetTH1("Evt_Events")->Fill(1);
		h1.GetTH1("Evt_CutFlow")->Fill("All", 1);

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
				looseElCol_isoMu.push_back(lepton);	
				if( lepton.Pt < 20) continue;
				if( lepton.Pt < 45)
					looseElCol_isoEl.push_back(lepton);
				else
					selElCol.push_back(lepton);	
			}
			// Muon selections
			if( lepton.LeptonType == 13 )
			{
				if( relIso1 > 0.2 ) continue;	
				if( abs(lepton.Eta) > 2.5 ) continue;
				if( lepton.Pt < 10) continue; 
					looseMuCol_isoEl.push_back(lepton);	
				if( lepton.Pt < 35)
					looseMuCol_isoMu.push_back(lepton);	
				else if( relIso1 < 0.125 && 
					 lepton.Pt >= 35 && 
					 lepton.MuType == 14 && 
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
		Jet jet1;
		if( selVertex.size() ){
			h1.GetTH1("Evt_CutFlow")->Fill("#geq1 goodVtx", 1);
			if( seljetCol.size() >= 3 ){ 	
				h1.GetTH1("Evt_CutFlow")->Fill("#geq3 Jets", 1);
				h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq3 Jets", 1);
			}
			// Lable the hardest non_bjet 
			int j1=-1;
			double pt1=0;
			const int size_seljetCol = seljetCol.size();	
			for( int i=0; i < size_seljetCol; i++){
				if( seljetCol[i].CombinedSVBJetTags < 0.679 && pt1 < seljetCol[i].Pt ){
					j1=i;
					pt1=seljetCol[i].Pt;
				}
			}
			jet1 = seljetCol[j1];
			if( bjetCol.size() >= 2 ){
				//isMuCh=true;	
				h1.GetTH1("Evt_CutFlow")->Fill("#geq2 bjets", 1);
				//h1.GetTH1("Evt_CutFlow_Mu")->Fill("#geq2 bjets", 1);
				// Lable bjet by Pt
				//sort2HighPt( bjetCol, bjet1, bjet2 );

				//h1.GetTH1("bJet12_Px")->Fill(bjet1.P4().Px());
				//h1.GetTH1("bJet12_Px")->Fill(bjet2.P4().Px());
				//h1.GetTH1("bJet12_Py")->Fill(bjet1.P4().Py());
				//h1.GetTH1("bJet12_Py")->Fill(bjet2.P4().Py());
				//h1.GetTH1("bJet12_Pz")->Fill(bjet1.P4().Pz());
				//h1.GetTH1("bJet12_Pz")->Fill(bjet2.P4().Pz());
			}

			if( bjetCol.size() == 2 ){	
				//is2bMuCh=true;	
				h1.GetTH1("Evt_CutFlow")->Fill("=2 bjets", 1);
				//h1.GetTH1("Evt_CutFlow_Mu")->Fill("=2 bjets", 1);
			}
		}	

  	} //// entry loop 
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
