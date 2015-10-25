#ifndef SELECTORELECTRON_H 
#define SELECTORELECTRON_H 

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Lepton.h"

class SelectorElectron{
    public:
        // constructor
        SelectorElectron()
        {
            min=0;
            max=1;
            hasCuts=false;
        }
        SelectorElectron( const edm::ParameterSet& iconfig, bool debug=false )
        {
            min=0;
            max=1;
            hasCuts=false;
            setCuts(iconfig, debug);
        }
        ~SelectorElectron(){};

        // functions
        void setCuts( const edm::ParameterSet& iconfig, bool debug=false )
        {
            if( hasCuts ) std::cout<<">> [WARNING] Cuts has set before!"<<std::endl;
            else hasCuts=true;

            iConfig = iconfig;
            leptype = iConfig.getParameter<std::string>("lepType");
            ConversionVeto = iConfig.getParameter<bool>("ConversionVeto");
            setPars( "lepEt"                     );
            setPars( "lepAbsEta"                 );
            setPars( "lepAbsEtaExclude"          );
            setPars( "lepRelIsoR03"              );
            setPars( "ElAbsTrackDxyPV"           );
            setPars( "EgammaMVATrig"             );
            setPars( "NumberOfExpectedInnerHits" );

            if( debug )
            {
                printf(">> [DEBUG] List current cuts in selector: %s\n", leptype.c_str());
                printf("%30s %10s %10s\n", "Selection", "Min", "Max");
                printCuts("ET",                        "lepEt"                     );
                printCuts("|Eta|",                     "lepAbsEta"                 );
                printCuts("|Eta| Gap exclude",         "lepAbsEtaExclude"          );
                printCuts("relIsoR03",                 "lepRelIsoR03"              );
                printCuts("ElAbsTrackDxyPV",           "ElAbsTrackDxyPV"           );
                printCuts("EgammaMVATrig",             "EgammaMVATrig"             );
                printCuts("NumberOfExpectedInnerHits", "NumberOfExpectedInnerHits" );
                printf("%30s %10d \n", "ConversionVeto",  ConversionVeto );
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
                std::cout<<">>         Please use SelectorElectron::setCuts(const edm::ParameterSet& iConfig)"<<std::endl;
                return false;
            }

            // Electron selections
            bool exclude=true;
            if( !pass( lepton.Et,                        "lepEt"                             )) return false;
            //if( !pass( getRelIsoR03(lepton),             "lepRelIsoR03"                      )) return false;
            if( !pass( lepton.RelIsoR03,                 "lepRelIsoR03"                      )) return false;
            if( !pass( fabs(lepton.Eta),                 "lepAbsEta"                         )) return false;
            if( !pass( fabs(lepton.Eta),                 "lepAbsEtaExclude",         exclude )) return false;
            if( !pass( fabs(lepton.ElTrackDxy_PV),       "ElAbsTrackDxyPV"                   )) return false;
            if( !pass( lepton.EgammaMVATrig,             "EgammaMVATrig"                     )) return false;
            if( !pass( lepton.NumberOfExpectedInnerHits, "NumberOfExpectedInnerHits"         )) return false;
            if( ConversionVeto && lepton.ElhasConv ) return false; // = no passConversinoVeto;
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
                printf(">> [ERROR] %s is not found in SelectorElectron::getCut(std::string)\n", parName.c_str());
            }
            return mapPars.find(parName)->second; 
        }

        //float getRelIsoR04( Lepton lepton )
        //{
        //    float a = lepton.ChargedHadronIsoR04 + lepton.NeutralHadronIsoR04 + lepton.PhotonIsoR04;
        //    float b = fabs(lepton.Pt);
        //    float reliso = a/b;
        //    return reliso;
        //}
        //float getRelIsoR03( Lepton lepton )
        //{
        //    float a = lepton.ChargedHadronIsoR03 + lepton.NeutralHadronIsoR03 + lepton.PhotonIsoR03;
        //    float b = fabs(lepton.Pt);
        //    float reliso = a/b;
        //    return reliso;
        //}
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
        bool ConversionVeto;
        int max;
        int min;
        std::string leptype;
        map<std::string, double> mapPars;
};
#endif
