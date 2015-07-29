#ifndef SELECTORMUON_H 
#define SELECTORMUON_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
		SelectorMuon( const edm::ParameterSet& iConfig )
		{
			min=0;
			max=1;
			setCuts(iConfig);
		}
 
		// functions
		void setCuts( const edm::ParameterSet& iConfig )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;
	
			pt[min] = iConfig.getParameter<double>("pt_Min");
			pt[max] = iConfig.getParameter<double>("pt_Max");	
			absEta[min] = iConfig.getParameter<double>("absEta_Min");	
			absEta[max] = iConfig.getParameter<double>("absEta_Max");	
			relIsoR04[min] = iConfig.getParameter<double>("relIsoR04_Min");	
			relIsoR04[max] = iConfig.getParameter<double>("relIsoR04_Max");	
			absMuInnerTrackDxy_PV[min] = iConfig.getParameter<double>("absMuInnerTrackDxy_PV_Min");
			absMuInnerTrackDxy_PV[max] = iConfig.getParameter<double>("absMuInnerTrackDxy_PV_Max");
			MuGlobalNormalizedChi2[min] = iConfig.getParameter<double>("MuGlobalNormalizedChi2_Min");
			MuGlobalNormalizedChi2[max] = iConfig.getParameter<double>("MuGlobalNormalizedChi2_Max");

			MuNMuonhits[min] = iConfig.getParameter<int>("MuNMuonhits_Min");
			MuNMuonhits[max] = iConfig.getParameter<int>("MuNMuonhits_Max");
		        MuNMatchedStations[min] = iConfig.getParameter<int>("MuNMatchedStations_Min");
		        MuNMatchedStations[max] = iConfig.getParameter<int>("MuNMatchedStations_Max");
			MuNTrackLayersWMeasurement[min] = iConfig.getParameter<int>("MuNTrackLayersWMeasurement_Min");
			MuNTrackLayersWMeasurement[max] = iConfig.getParameter<int>("MuNTrackLayersWMeasurement_Max");

			checkGlobalMuon = iConfig.getParameter<bool>("CheckGlobalMuon");
		}
		bool isPass( Lepton lepton )
		{
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorMuon::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return
			}

			if( !pass( lepton.Pt,		pt[min],     	pt[max] ))	  return false;
			if( !pass( relIsoR04(lepton),	relIsoR04[min],	relIsoR04[max] )) return false;
			if( !pass( absf(lepton.Eta),	absEta[min], 	absEta[max] ))	  return false;
			if( !pass( absf(lepton.MuInnerTrackDxy_PV), 	absMuInnerTrackDxy_PV[min], 	absMuInnerTrackDxy_PV[max] ))      return false;
			if( !pass( lepton.MuNMuonhits,			MuNMuonhits[min],		MuNMuonhits[max] ))		   return false;
			if( !pass( lepton.MuNMatchedStations,		MuNMatchedStations[min],	MuNMatchedStations[max] ))	   return false;
			if( !pass( lepton.MuGlobalNormalizedChi2,	MuGlobalNormalizedChi2[min],	MuGlobalNormalizedChi2[max] ))	   return false;
			if( !pass( lepton.MuNTrackLayersWMeasurement,	MuNTrackLayersWMeasurement[min],MuNTrackLayersWMeasurement[max] )) return false;
			if( checkGlobalMuon ){
				 if( (lepton.MuType&0x02) == 0 ) return false;
			}
			return true;
		}

		// Helper
		bool pass( int value, int min, int max )
		{
			if( value > max || value < min ) return false;
			return true;
		}
		bool pass( float value, double min, double max )
		{
			if( value > max || value < min ) return false;
			return true;
		}
		bool pass( bool value, double cut )
		{
			if( value != cut ) return false;
			return true;
		}
		float relIsoR04( Lepton lepton )
		{
			float a = lepton.ChargedHadronIsoR04 + lepton.NeutralHadronIsoR04 + lepton.PhotonIsoR04;
			float b = fabs(lepton.Pt);
			float reliso = a/b;
			return reliso;	
		}
		float relIsoR03( Lepton lepton )
		{
			float a = lepton.ChargedHadronIsoR03 + lepton.NeutralHadronIsoR03 + lepton.PhotonIsoR03;
			float b = fabs(lepton.Pt);
			float reliso = a/b;
			return reliso;	
		}

	private:		
		bool hasCuts;
		const int max;
		const int min;
		const double pt[2];
		const double relIsoR04[2]
		const double absEta[2];
		const double absMuInnerTrackDxy_PV[2];
		const int    MuNMuonhits[2];
		const int    MuNMatchedStations[2];
		const int    MuGlobalNormalizedChi2[2];
		const int    MuNTrackLayersWMeasurement[2];
		const bool   checkGlobalMuon;
};
#endif
