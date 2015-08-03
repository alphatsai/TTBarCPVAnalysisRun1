#ifndef SELECTORVERTEX_H 
#define SELECTORVERTEX_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Vertex.h"

class SelectorVertex{
	public:
		// constructor
		SelectorVertex()
		{
			min=0;
			max=1;
			hasCuts=false;
		}
		SelectorVertex( const edm::ParameterSet& iconfig, bool debug=false )
		{
			min=0;
			max=1;
			hasCuts=false;
			setCuts(iconfig, debug);
		}
		~SelectorVertex(){};
 
		// functions
		void setCuts( const edm::ParameterSet& iconfig, bool debug=false )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;

			iConfig = iconfig;
	
			vtxtype	         = iConfig.getParameter<std::string>("vtxType");
			CheckIsFake      = iConfig.getParameter<bool>("CheckIsFake");
			CheckIsOfflinePV = iConfig.getParameter<bool>("CheckIsOfflinePV");
			setPars( "vtxNdof" );
			setPars( "vtxAbsZ" );
			setPars( "vtxRho"  );

			if( debug )
			{
				printf(">> [DEBUG] List current cuts in selector: %s\n", vtxtype.c_str());
				printf("%30s %10s %10s\n", "Selection", "Min", "Max");
				printf("%30s %10d \n",     "CheckIsOfflinePV", CheckIsOfflinePV );
				printf("%30s %10d \n",     "CheckIsFake",      CheckIsFake      );
				printCuts("Ndof",          "vtxNdof"                            );
				printCuts("|z|",           "vtxAbsZ"                            );
				printCuts("Rho",           "vtxRho"                             );
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


		bool isPass( Vertex vtx )
		{
			// std::cout<<">> isPass() "<<std::endl;
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorVertex::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return false;
			}
			if( CheckIsFake )     { if( vtx.isFake    ) return false; }
			if( CheckIsOfflinePV ){ if( vtx.Type != 1 ) return false; }
			if( !pass( vtx.Ndof,    "vtxNdof" )) return false;
			if( !pass( fabs(vtx.z), "vtxAbsZ" )) return false;
			if( !pass( vtx.Rho,     "vtxRho"  )) return false;
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
			//std::cout<<">>   Pass(): "<<parName<<std::endl;
			return pass( value, getCut(parName+"Min"), getCut(parName+"Max")); 
		}

		double getCut( std::string parName )
		{ 
			if( mapPars.find(parName) == mapPars.end() ){
				printf(">> [ERROR] %s is not found in SelectorVertex::getCut(std::string)\n", parName.c_str());
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
		std::string vtxtype;
		bool   CheckIsOfflinePV;
		bool   CheckIsFake;
		map<std::string, double> mapPars;
};
#endif
