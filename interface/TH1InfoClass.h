#ifndef TH1INFOCLASS_H
#define TH1INFOCLASS_H

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH1Info.h"

using namespace std;

template<typename TH1> 
class TH1InfoClass
{
	public:
		TH1InfoClass( bool deBug=false );
		void addNewTH1( std::string name, std::string title, std::string xtitle, std::string ytitle, std::string xunit, std::string yunit, int bin, double min, double max);
		//void initTH1Info(); // Be writen in custumer way
		void defaultTH1Info(); 
		void CreateTH1();
		void CreateTH1( TFile* f, std::string dir_name=""); 
		void CreateTH1( edm::Service<TFileService> f); 
		void SetTitles(); 
		void Sumw2();
		TH1* GetTH1(std::string name);
		TH1Info GetInfo(std::string name);
		TH1Info GetInfo(int index);
		
		bool debug;
		int size;

	private:
		map<std::string, TH1*> mapTH1;
		map<std::string, int> indexTH1;
		vector<TH1Info> Info; 
};

template<typename TH1>
TH1InfoClass<TH1>::TH1InfoClass( bool deBug )
{
	debug=deBug;	
	defaultTH1Info();
	size=Info.size();
}

template<typename TH1>
void TH1InfoClass<TH1>::addNewTH1(std::string name, std::string title, std::string xtitle, std::string ytitle, std::string xunit, std::string yunit, int bin, double min, double max)
{
	Info.push_back( TH1Info( name, title, xtitle, ytitle, xunit, yunit, bin, min, max ) );
	size=Info.size();
	if( debug ) printf(">> [DEBUG] Add new TH1 %-24s bin/min/max:%d/%.2f/%.2f, new size is %d\n", name.c_str(), bin, min, max, size);
}
template<typename TH1>
void TH1InfoClass<TH1>::defaultTH1Info()
{
	// Info.push_back( TH1Info( Name,		Title,			xTitle, 	yTitle, 	xUnit,	yUnit, 	Bin, Min, Max ) );
	Info.push_back( TH1Info( "Jet_Pt",		"pT of Jet",			  "p_{T}(j)", 		"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "Jet_Px",		"px of Jet",			  "p_{x}(j)", 		"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "Jet_Py",		"py of Jet",			  "p_{y}(j)", 		"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "Jet_Pz",		"pz of Jet",			  "p_{z}(j)", 		"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "Jet_M",		"Mass of Jet",			  "Mass(j)", 		"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "Jet_E",		"Energy of Jet",		  "Energy(j)",		"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "Jet_Eta",		"Eta of Jet",	 		  "#eta(j)", 		"Yields", 	"", 	"",		100, -5, 5 ) );
	Info.push_back( TH1Info( "Jet_Phi",		"Phi of Jet",	 		  "#phi(j)", 		"Yields", 	"", 	"",		64, -3.2,   3.2 ) );
	Info.push_back( TH1Info( "Jet_BTag",	"Jet b-tagged",	 		  "bTag", 		    "Yields", 	"", 	"",		100, 0,  1  ) );
	Info.push_back( TH1Info( "SelJet_Pt",	"pT of selected Jet",	  "p_{T}(selected j)", 	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "SelJet_Px",	"px of selected Jet",	  "p_{x}(selected j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "SelJet_Py",	"py of selected Jet",	  "p_{y}(selected j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "SelJet_Pz",	"pz of selected Jet",	  "p_{z}(selected j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "SelJet_M",	"Mass of selected Jet",	  "Mass(selected j)", 	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "SelJet_E",	"Energy of selected Jet", "Energy(selected j)",	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "SelJet_Eta",	"Eta of selected Jet",	  "#eta(selected j)", 	"Yields", 	"", 	"",		100, -5, 5 ) );
	Info.push_back( TH1Info( "SelJet_Phi",	"Phi of selected Jet",	  "#phi(selected j)", 	"Yields", 	"", 	"",		64, -3.2,   3.2 ) );
	Info.push_back( TH1Info( "SelJet_BTag",	"Selected Jet b-tagged",  "bTag", 		    "Yields", 	"", 	"",		100, 0,   1 ) );
	Info.push_back( TH1Info( "bJet_Pt",     "pT of b-Jet",	  "p_{T}(B-tagged j)", 	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "bJet_Px",     "px of b-Jet",	  "p_{x}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "bJet_Py",     "py of b-Jet",	  "p_{y}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "bJet_Pz",     "pz of b-Jet",	  "p_{z}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "bJet_M",	"Mass of b-Jet",  "Mass(B-tagged j)", 	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "bJet_E",	"Energy of b-Jet","Energy(B-tagged j)",	"Yields", 	"GeV", 	"",		500, 0,   500 ) );
	Info.push_back( TH1Info( "bJet_Eta",	"Eta of b-Jet",	  "#eta(B-tagged j)", 	"Yields", 	"", 	"",		100, -5, 5 ) );
	Info.push_back( TH1Info( "bJet_Phi",	"Phi of b-Jet",	  "#phi(B-tagged j)", 	"Yields", 	"", 	"",		64, -3.2,   3.2 ) );
	Info.push_back( TH1Info( "bJet_BTag",	"b-Jet b-tagged",  "bTag", 		    "Yields", 	"", 	"",	100, 0,   1 ) );
	Info.push_back( TH1Info( "bJet12_Px",     "px of b-Jet",  "p_{x}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "bJet12_Py",     "py of b-Jet",  "p_{y}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "bJet12_Pz",     "pz of b-Jet",  "p_{z}(B-tagged j)", 	"Yields", 	"GeV", 	"",		1000, -500,   500 ) );
	Info.push_back( TH1Info( "Evt_Events",	"",	  	  "", 				"Events", 	"", 	"",		1, 1,   2) );
	Info.push_back( TH1Info( "Evt_Channel",	"",	  	  "", 				"", 		"", 	"",		2, 0,   2) );
	Info.push_back( TH1Info( "Evt_CutFlow", "",        	  "",     			"Evetns", 	"", 	"",     8, 0, 8 ));
	Info.push_back( TH1Info( "Evt_NJets",	"Num. of jets",	  		  "N(j)", 			"Events", 	"", 	"",		20, 0,   20) );
	Info.push_back( TH1Info( "Evt_NSelJets","Num. of selected jets",  "N(selected j).", "Events", 	"", 	"",		20, 0,   20) );
	Info.push_back( TH1Info( "Evt_NbJets",	"Num. of b-jets",	  	  "N(B-tagged j)", 	"Events", 	"", 	"",		20, 0,   20) );
	Info.push_back( TH1Info( "Evt_O7",	"O7",	  			"O_{7}", 	"Events", 	"", 	"",		40, -2,   2) );
	Info.push_back( TH1Info( "Evt_O7Asym",	"O7 Asymetric",	  	"", 	"Events", 	"", 	"",		2, 0,   2) );
	Info.push_back( TH1Info( "Evt_O7_term1","Pz(pbj1-pbj2)",	"GeV", 	"Events", 	"", 	"",		1000, -500,   500) );
	Info.push_back( TH1Info( "Evt_O7_term2","Pz(pbj1xpbj2)",	"GeV", 	"Events", 	"", 	"",		1000, -500,   500) );
}

//* Create Histogram
template<typename TH1> 
void TH1InfoClass<TH1>::CreateTH1()
{
	if( debug ) printf(">> [DEBUG] %d TH1s will be created...\n", size);
	for(int i=0; i<size; i++){
		if( debug ) printf(">> [DEBUG] #%3d %-24s created. Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].Bin, Info[i].Min, Info[i].Max);
		indexTH1[Info[i].Name] = i;
		mapTH1[Info[i].Name] = new TH1(Info[i].Name.c_str(),"",Info[i].Bin, Info[i].Min, Info[i].Max);
	}
}
template<typename TH1> 
void TH1InfoClass<TH1>::CreateTH1( edm::Service<TFileService> f )
{
	if( debug ) printf(">> [DEBUG] %d TH1s will be created...\n", size);
	for(int i=0; i<size; i++){
		if( debug ) printf(">> [DEBUG] #%3d %-24s created. Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].Bin, Info[i].Min, Info[i].Max);
		indexTH1[Info[i].Name] = i;
		mapTH1[Info[i].Name] = f->make<TH1>(Info[i].Name.c_str(),"",Info[i].Bin, Info[i].Min, Info[i].Max);
	}
}
template<typename TH1> 
void TH1InfoClass<TH1>::CreateTH1( TFile* f, std::string dir_name)
{
	printf(">> %d TH1s will be got...\n", size);
	for(int i=0; i<size; i++){ 
		if( debug ) printf(">> [DEBUG] #%3d %-24s created. Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].Bin, Info[i].Min, Info[i].Max);
		indexTH1[Info[i].Name] = i;
		mapTH1[Info[i].Name] =(TH1*)f->Get( (dir_name+Info[i].Name).c_str() );
	}
}

//* Set some option for Histogram
template<typename TH1> 
void TH1InfoClass<TH1>::SetTitles(){
	for(int i=0; i<size; i++){ 
		//mapTH1[Info[i].Name]->SetTile(Info[i]._title.c_str());
		string xt, yt;
		if( Info[i].xUnit.size() == 0 ) xt = Info[i].xTitle;
		else xt = Info[i].xTitle+" ["+Info[i].xUnit+"]";	
		if( Info[i].yUnit.size() == 0 ) yt = Info[i].yTitle;
		else yt = Info[i].yTitle+" ["+Info[i].yUnit+"]";	
		mapTH1[Info[i].Name]->SetXtitle( xt.c_str() );
		mapTH1[Info[i].Name]->SetYtitle( yt.c_str() );
	}
}

template<typename TH1> 
void TH1InfoClass<TH1>::Sumw2(){
	for(int i=0; i<size; i++){ 
		mapTH1.find(Info[i].Name)->second->Sumw2();
	}
}

//* Get Histogram
template<typename TH1> 
TH1* TH1InfoClass<TH1>::GetTH1(std::string name){
	if( mapTH1.find(name) == mapTH1.end() ){
		printf(">> [ERROR] %s is not found in TH1InfoClass::GetTH1(std::string)\n", name.c_str());
		return NULL;
	}else{
		return mapTH1.find(name)->second;
	}
}

//* Get Info variables
template<typename TH1> 
TH1Info TH1InfoClass<TH1>::GetInfo(std::string name){
	TH1Info info;
	if( mapTH1.find(name) == mapTH1.end() ){
		printf(">> [ERROR] %s is not found in TH1InfoClass::GetInfo(std::string)\n", name.c_str());
	}else{	
		info=Info[indexTH1.find(name)->second];
	}
	return info;
}
template<typename TH1> 
TH1Info TH1InfoClass<TH1>::GetInfo(int index){
	TH1Info info;
	if( index >= size ){
		printf(">> [ERROR] %d is over the size(%d) in TH1InfoClass::GetInfo(int)\n", index, size);
	}else{
		info=Info[index];
	}
	return info;
}

#endif
//
