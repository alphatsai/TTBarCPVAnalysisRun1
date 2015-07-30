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
			setCuts(iconfig, debug);
		}
		~SelectorMuon(){ return; }
 
		// functions
		void setCuts( const edm::ParameterSet& iconfig, bool debug=false )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;

			iConfig = iconfig;
	
			leptype	                        = iConfig.getParameter<std::string>("lepType");
			pt[min]                         = iConfig.getParameter<double>("lepPtMin");
			pt[max]                         = iConfig.getParameter<double>("lepPtMax");	
			absEta[min]                     = iConfig.getParameter<double>("lepAbsEtaMin");	
			absEta[max]                     = iConfig.getParameter<double>("lepAbsEtaMax");	
			relIsoR04[min]                  = iConfig.getParameter<double>("lepRelIsoR40Min");	
			relIsoR04[max]                  = iConfig.getParameter<double>("lepRelIsoR40Max");	
			absMuInnerTrackDxy_PV[min]      = iConfig.getParameter<double>("MuAbsInnerTrackDxyPVMin");
			absMuInnerTrackDxy_PV[max]      = iConfig.getParameter<double>("MuAbsInnerTrackDxyPVMax");
			MuGlobalNormalizedChi2[min]     = iConfig.getParameter<double>("MuGlobalNormalizedChi2Min");
			MuGlobalNormalizedChi2[max]     = iConfig.getParameter<double>("MuGlobalNormalizedChi2Max");

			MuNMuonhits[min]                = iConfig.getParameter<int>("MuNMuonhitsMin");
			MuNMuonhits[max]                = iConfig.getParameter<int>("MuNMuonhitsMax");
		        MuNMatchedStations[min]         = iConfig.getParameter<int>("MuNMatchedStationsMin");
		        MuNMatchedStations[max]         = iConfig.getParameter<int>("MuNMatchedStationsMax");
			MuNTrackLayersWMeasurement[min] = iConfig.getParameter<int>("MuNTrackLayersWMeasurementMin");
			MuNTrackLayersWMeasurement[max] = iConfig.getParameter<int>("MuNTrackLayersWMeasurementMax");

			checkGlobalMuon              = iConfig.getParameter<bool>("CheckGlobalMuon");
			checkTrackerMuon             = iConfig.getParameter<bool>("CheckTrackerMuon");

			if( debug )
			{
				printf(">> [DEBUG] List current cuts in selector: %s", leptype.c_str());
				printf("%15s %10s %10s\n", "Selection", "Min", "Max");
				printCuts("pT", pt);
				printCuts("|Eta|", absEta);
				printCuts("relIsoR04", relIsoR04);
				printCuts("|MuInnerTrackDxy_PV|", absMuInnerTrackDxy_PV);
				printCuts("MuGlobalNormalizedChi2", MuGlobalNormalizedChi2);
				printCuts("MuNMuonhits", MuNMuonhits);
				printCuts("MuNMatchedStations", MuNMatchedStations);
				printCuts("MuNTrackLayersWMeasurement", MuNTrackLayersWMeasurement);
				printCuts("checkGlobalMuon", checkGlobalMuon);
				printCuts("checkTrackerMuon", checkTrackerMuon);
			}
		}
		template<typename ParType>
		void setPars( ParType* par, std::string parMin, std::string parMax )
		{
			par[min] = iConfig.getParameter<ParType*>(parMin.c_str());
			par[max] = iConfig.getParameter<ParType*>(parMax.c_str());
			mapPars[parMin]=par[min];
			mapPars[parMax]=par[max];
		}


		bool isPass( Lepton lepton )
		{
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorMuon::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return false;
			}

			// Muon selections
			if( !pass( lepton.Pt,				pt[min],     			pt[max] 			)) return false;
			if( !pass( getRelIsoR04(lepton),		relIsoR04[min],			relIsoR04[max] 			)) return false;
			if( !pass( fabs(lepton.Eta),			absEta[min], 			absEta[max] 			)) return false;
			if( !pass( fabs(lepton.MuInnerTrackDxy_PV), 	absMuInnerTrackDxy_PV[min], 	absMuInnerTrackDxy_PV[max] 	)) return false;
			if( !pass( lepton.MuNMuonhits,			MuNMuonhits[min],		MuNMuonhits[max] 		)) return false;
			if( !pass( lepton.MuNMatchedStations,		MuNMatchedStations[min],	MuNMatchedStations[max] 	)) return false;
			if( !pass( lepton.MuGlobalNormalizedChi2,	MuGlobalNormalizedChi2[min],	MuGlobalNormalizedChi2[max] 	)) return false;
			if( !pass( lepton.MuNTrackLayersWMeasurement,	MuNTrackLayersWMeasurement[min],MuNTrackLayersWMeasurement[max] )) return false;
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
		//bool pass( int value, int min, int max )
		//bool pass( int value, double min, double max )
		//{
		//	if( value > max || value < min ) return false;
		//	return true;
		//}
		//bool pass( float value, double min, double max )
		//{
		//	if( value > max || value < min ) return false;
		//	return true;
		//}
		template<typename parType>
		bool pass( parType value, double min, double max )
		{
			if( value > max || value < min ) return false;
			return true;
		}
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
		void printCuts( string cutName, double* cut )
		{
			printf("%25s %10.2f %10.2f\n", cutName.c_str(), cut[0], cut[1]);
		}
		void printCuts( string cutName, int* cut )
		{
			printf("%25s %10d %10d\n", cutName.c_str(), cut[0], cut[1]);
		}
		void printCuts( string cutName, bool cut )
		{
			printf("%25s %10d \n", cutName.c_str(), cut);
		}

	private:	
		edm::ParameterSet iConfig;
		bool hasCuts;
		int max;
		int min;
		std::string leptype;
		bool   checkGlobalMuon;
		bool   checkTrackerMuon;
		double pt[2];
		double relIsoR04[2];
		double absEta[2];
		double absMuInnerTrackDxy_PV[2];
		double MuGlobalNormalizedChi2[2];
		double MuNMuonhits[2];
		double MuNMatchedStations[2];
		double MuNTrackLayersWMeasurement[2];

		map<std::string, double> mapPars;
};
#endif
