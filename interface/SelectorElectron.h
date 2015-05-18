#ifndef SELECTORLEP_H 
#define SELECTORLEP_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "Lepton.h"

class SelectorElectron{
	public:
		// constructor
		SelectorElectron()
		{
			hasCuts=false;
		}
		SelectorElectron( const edm::ParameterSet& iConfig )
		{
			setCuts(iConfig);
		}
 
		// functions
		void setCuts( const edm::ParameterSet& iConfig )
		{
			if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
			else hasCuts=true;
	
			pt_max = iConfig.getParameter<double>("ptMax");	
			pt_min = iConfig.getParameter<double>("ptMin");
			absEta_max = iConfig.getParameter<double>("absEtaMax");	
			absEta_min = iConfig.getParameter<double>("absEtaMin");	
			absGapEta_max = iConfig.getParameter<double>("absGapEtaMax");	
			absGapEta_min = iConfig.getParameter<double>("absGapEtaMin");	
			relIso_max = iConfig.getParameter<double>("relIsoMax");	
			relIso_min = iConfig.getParameter<double>("relIsoMin");	
		}
		bool isPass( Lepton lepton )
		{
			if( !hasCuts ){ 
				std::cout<<">> [ERROR] Not cut set yet"<<std::endl;
				std::cout<<">>         Please use SelectorElectron::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
				return
			}

			if( !pass( lepton.Pt,      pt_min,     pt_max ))     return false;
			if( !pass( lepton.Eta,     absEta_min, absEta_max )) return false;
			if( !pass( relIso(lepton), relIso_min, relIso_max )) return false;
			if( lepton.Eta > absGapEta_min && lepton.Eta < absGapEta_max ) return false;
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
		float relIso( Lepton lepton )
		{
			float a = lepton.ChargedHadronIso + lepton.NeutralHadronIso + lepton.PhotonIso;
			float b = fabs(lepton.Pt);
			float reliso = a/b;
			return reliso;	
		}

	private:		
		bool hasCuts;
		const double pt_max;
		const double pt_min;
		const double absEta_max;
		const double absEta_min;
		const double absGapEta_max;
		const double absGapEta_min;
		const double relIso_max;
		const double relIso_min;
};
#endif
