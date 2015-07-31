#ifndef SELECTORMUON_H 
#define SELECTORMUON_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Lepton.h"

class SelectorMuon{
	public:
		// constructor
		SelectorMuon()
		{
			min=0;
			max=1;
			hasCuts=false;
		}
		SelectorMuon( const edm::ParameterSet& iconfig, bool debug=false )
		{
			min=0;
			max=1;
			hasCuts=false;
			setCuts(iconfig, debug);
		}
		~SelectorMuon(){};
 
		// functions
		void setCuts( const edm::ParameterSet& iconfig, bool debug=false )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;

			iConfig = iconfig;
	
			leptype	         = iConfig.getParameter<std::string>("lepType");
			checkGlobalMuon  = iConfig.getParameter<bool>("CheckGlobalMuon");
			checkTrackerMuon = iConfig.getParameter<bool>("CheckTrackerMuon");
			setPars( "lepPt"                      );
			setPars( "lepAbsEta"                  );
			setPars( "lepRelIsoR04"               );
			setPars( "MuAbsInnerTrackDxyPV"       );
			setPars( "MuGlobalNormalizedChi2"     );
			setPars( "MuNMuonhits"                );
			setPars( "MuNMatchedStations"         );
			setPars( "MuNTrackLayersWMeasurement" );

			if( debug )
			{
				printf(">> [DEBUG] List current cuts in selector: %s\n", leptype.c_str());
				printf("%30s %10s %10s\n", "Selection", "Min", "Max");
				printf("%30s %10d \n",     "checkGlobalMuon",  checkGlobalMuon );
				printf("%30s %10d \n",     "checkTrackerMuon", checkTrackerMuon );
				printCuts("pT",                         "lepPt"                         );
				printCuts("|Eta|",                      "lepAbsEta"                     );
				printCuts("relIsoR04",                  "lepRelIsoR04"                  );
				printCuts("|MuInnerTrackDxy_PV|",       "MuAbsInnerTrackDxyPV"          );
				printCuts("MuGlobalNormalizedChi2",     "MuGlobalNormalizedChi2"        );
				printCuts("MuNMuonhits",                "MuNMuonhits"                   );
				printCuts("MuNMatchedStations",         "MuNMatchedStations"            );
				printCuts("MuNTrackLayersWMeasurement", "MuNTrackLayersWMeasurement"    );
			}
		}
		void setPars( std::string parName )
		{
			setPars( parName+"Min", parName+"Max");
		}
		void setPars( std::string parMin, std::string parMax )
		{
			mapPars[parMin]= iConfig.getParameter<double>(parMin.c_str());
			mapPars[parMax]= iConfig.getParameter<double>(parMax.c_str());
		}


		bool isPass( Lepton lepton )
		{
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorMuon::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return false;
			}

			// Muon selections
			if( !pass( lepton.Pt,				"lepPt"                       )) return false;
			if( !pass( getRelIsoR04(lepton),		"lepRelIsoR40"                )) return false;
			if( !pass( fabs(lepton.Eta),			"lepAbsEta"                   )) return false;
			if( !pass( fabs(lepton.MuInnerTrackDxy_PV), 	"MuAbsInnerTrackDxyPV"        )) return false;
			if( !pass( lepton.MuNMuonhits,			"MuNMuonhits"                 )) return false;
			if( !pass( lepton.MuNMatchedStations,		"MuNMatchedStations"          )) return false;
			if( !pass( lepton.MuGlobalNormalizedChi2,	"MuGlobalNormalizedChi2"      )) return false;
			if( !pass( lepton.MuNTrackLayersWMeasurement,	"MuNTrackLayersWMeasurement"  )) return false;
			if( checkGlobalMuon && !checkTrackerMuon ){
				 if( (lepton.MuType&0x02) == 0 ) return false; // Only it's global muon
			}else if( !checkGlobalMuon && checkTrackerMuon ){
				 if( (lepton.MuType&0x04) == 0 ) return false; // Only it's tracker muon
			}else if( checkGlobalMuon && checkTrackerMuon ){
				 if( (lepton.MuType&0x02) == 0 && (lepton.MuType&0x04) == 0 ) return false; // Either tracker or global muon
			}

			return true;
		}

		// Helper
		template<typename parType>
		bool pass( parType value, double min, double max )
		{
			if( value > min && value < max ) return true;
			return false;
		}
		template<typename parType>
		bool pass( parType value, std::string parName )
		{ 
			return pass( value, getCut(parName+"Min"), getCut(parName+"Max")); 
		}

		double getCut( std::string parName ){ return mapPars.find(parName)->second; }

		float getRelIsoR04( Lepton lepton )
		{
			float a = lepton.ChargedHadronIsoR04 + lepton.NeutralHadronIsoR04 + lepton.PhotonIsoR04;
			float b = fabs(lepton.Pt);
			float reliso = a/b;
			return reliso;	
		}
		float getRelIsoR03( Lepton lepton )
		{
			float a = lepton.ChargedHadronIsoR03 + lepton.NeutralHadronIsoR03 + lepton.PhotonIsoR03;
			float b = fabs(lepton.Pt);
			float reliso = a/b;
			return reliso;	
		}
		void printCuts( string cutName, std::string parName )
		{
			std::string parMin=parName+"Min";
			std::string parMax=parName+"Max";
			if( getCut(parMin) <= -100000 && getCut(parMax) >= 100000 )      printf("%30s %10s %10s\n",   cutName.c_str(), "nan", "nan" );
			else if( getCut(parMin) <= -100000 && getCut(parMax) <  100000 ) printf("%30s %10s %10.2f\n", cutName.c_str(), "nan", getCut(parMax) );
			else if( getCut(parMin) >  -100000 && getCut(parMax) >= 100000 ) printf("%30s %10.2f %10s\n", cutName.c_str(),  getCut(parMin), "nan");
			else printf("%30s %10.2f %10.2f\n", cutName.c_str(),  getCut(parMin), getCut(parMax) );
		}

	private:	
		edm::ParameterSet iConfig;
		bool hasCuts;
		int max;
		int min;
		std::string leptype;
		bool   checkGlobalMuon;
		bool   checkTrackerMuon;
		map<std::string, double> mapPars;
};
#endif
