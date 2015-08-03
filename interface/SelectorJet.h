#ifndef SELECTORJET_H 
#define SELECTORJET_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Jet.h"

class SelectorJet{
	public:
		// constructor
		SelectorJet()
		{
			min=0;
			max=1;
			hasCuts=false;
		}
		SelectorJet( const edm::ParameterSet& iconfig, bool debug=false )
		{
			min=0;
			max=1;
			hasCuts=false;
			setCuts(iconfig, debug);
		}
		~SelectorJet(){};
 
		// functions
		void setCuts( const edm::ParameterSet& iconfig, bool debug=false )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;

			iConfig = iconfig;
	
			jettype	= iConfig.getParameter<std::string>("jetType");
			setPars( "jetPt"                 );
			setPars( "jetAbsEta"             );
			setPars( "jetNConstituents"      );
			setPars( "jetCHF"                );
			setPars( "jetNCH"                );
			setPars( "jetNEF"                );
			setPars( "jetNHF"                );
			setPars( "jetCEF"                );
			setPars( "jetCombinedSVBJetTags" );

			if( debug )
			{
				printf(">> [DEBUG] List current cuts in selector: %s\n", jettype.c_str());
				printf("%30s %10s %10s\n", "Selection", "Min", "Max");
				printCuts("pT",                 "jetPt"                 );
				printCuts("|Eta|",              "jetAbsEta"             );
				printCuts("NConstituents",      "jetNConstituents"      );
				printCuts("CHF",                "jetCHF"                );
				printCuts("NCH",                "jetNCH"                );
				printCuts("NEF",                "jetNEF"                );
				printCuts("NHF",                "jetNHF"                );
				printCuts("CEF",                "jetCEF"                );
				printCuts("CombinedSVBJetTags", "jetCombinedSVBJetTags" );
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

		bool isPass( Jet jet )
		{
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorJet::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return false;
			}

			// Jet selections
			if( !pass( jet.Pt,                 "jetPt"                 )) return false;
			if( !pass( fabs(jet.Eta),          "jetAbsEta"             )) return false;
			if( !pass( jet.NConstituents,      "jetNConstituents"      )) return false;
			if( !pass( jet.CHF,                "jetCHF"                )) return false;
			if( !pass( jet.NCH,                "jetNCH"                )) return false;
			if( !pass( jet.NEF,                "jetNEF"                )) return false;
			if( !pass( jet.NHF,                "jetNHF"                )) return false;
			if( !pass( jet.CEF,                "jetCEF"                )) return false;
			if( !pass( jet.CombinedSVBJetTags, "jetCombinedSVBJetTags" )) return false;
			return true;
		}

		// Helper
		template<typename parType>
		bool pass( parType value, double min, double max, bool isExclude=false )
		{
			bool ispass=false;
			if( isExclude ){
				if( value <= min || value >= max ) ispass=true; 
			}else{	
				if( value >  min && value <  max ) ispass=true; 
			}
			return ispass;
		}
		template<typename parType>
		bool pass( parType value, std::string parName, bool isExclude=false )
		{ 
			return pass( value, getCut(parName+"Min"), getCut(parName+"Max"), isExclude ); 
		}

		double getCut( std::string parName )
		{ 
			if( mapPars.find(parName) == mapPars.end() ){
				printf(">> [ERROR] %s is not found in SelectorJet::getCut(std::string)\n", parName.c_str());
			}
			return mapPars.find(parName)->second; 
		}

		void printCuts( string cutName, std::string parName )
		{
			std::string parMin=parName+"Min";
			std::string parMax=parName+"Max";
			if( getCut(parMin) <= -100000 && getCut(parMax) >= 100000 )      printf("%30s %10s %10s\n",   cutName.c_str(), "nan", "nan" );
			else if( getCut(parMin) <= -100000 && getCut(parMax) <  100000 ) printf("%30s %10s %10.3f\n", cutName.c_str(), "nan", getCut(parMax) );
			else if( getCut(parMin) >  -100000 && getCut(parMax) >= 100000 ) printf("%30s %10.3f %10s\n", cutName.c_str(),  getCut(parMin), "nan");
			else printf("%30s %10.3f %10.3f\n", cutName.c_str(),  getCut(parMin), getCut(parMax) );
		}

	private:	
		edm::ParameterSet iConfig;
		bool hasCuts;
		int max;
		int min;
		std::string jettype;
		map<std::string, double> mapPars;
};
#endif
