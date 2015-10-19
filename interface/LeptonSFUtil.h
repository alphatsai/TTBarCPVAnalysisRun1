#ifndef LEPTONSFUTIL_H 
#define LEPTONSFUTIL_H

#include <string>
#include "Lepton.h"

class LeptonSFUtil
{
    public:
        LeptonSFUtil(){};
        ~LeptonSFUtil(){};
        float getSF_TightElectron( Lepton electron, int shift=0 )
        {
            if( electron.LeptonType != 11 )
            { 
                std::cout<<">> [ERROR] LeptonSFUtil::getElectronSF only for Lepton::LeptonType==11(electron)"<<std::endl;
                return 1.;
            }

            float pt =electron.Pt;
            float eta=fabs(electron.Eta);
            float ptMin[]={10,15,20,30,40,50};
            float ptMax[]={15,20,30,40,50,200};
            float etaMin[]={0.0,0.8,1.442,1.556,2.0};
            float etaMax[]={0.8,1.442,1.556,2.0,2.5};
            const int sizePt  = int(sizeof(ptMin)/sizeof(ptMin[0]));
            const int sizeEta = int(sizeof(etaMin)/sizeof(etaMin[0]));

            float nominal[sizeEta][sizePt]={ 
            //pt     0,    1,    2,    3,    4,    5
                {0.827,0.924,0.960,0.978,0.981,0.982}, // eta 0
                {0.948,0.932,0.936,0.958,0.969,0.969}, // eta 1 
                {1.073,0.808,0.933,0.907,0.904,0.926}, // eta 2
                {0.854,0.853,0.879,0.909,0.942,0.957}, // eta 3
                {1.007,0.903,0.974,0.987,0.991,0.999}  // eta 4
            };
            float sigmaNeg[sizeEta][sizePt]={ 
            //pt     0,    1,    2,    3,    4,    5
                {0.021,0.010,0.003,0.001,0.001,0.002}, // eta 0
                {0.023,0.012,0.004,0.002,0.001,0.002}, // eta 1 
                {0.107,0.042,0.017,0.008,0.004,0.011}, // eta 2
                {0.047,0.022,0.007,0.003,0.002,0.004}, // eta 3
                {0.046,0.029,0.004,0.004,0.003,0.005}  // eta 4
            };
            float sigmaPos[sizeEta][sizePt]={ 
            //pt     0,    1,    2,    3,    4,    5
                {0.021,0.010,0.003,0.001,0.001,0.002}, // eta 0
                {0.024,0.012,0.004,0.002,0.001,0.002}, // eta 1 
                {0.117,0.045,0.015,0.008,0.004,0.011}, // eta 2
                {0.048,0.022,0.007,0.003,0.002,0.004}, // eta 3
                {0.047,0.029,0.004,0.004,0.003,0.005}  // eta 4
            };

            int iPt(-1);
            int iEta(-1);
            for( int i=0; i<sizeEta; i++ )
            {
                if( eta >= etaMin[i] && eta < etaMax[i] ){ iEta=i; break; }
                else iEta=i;
            } 
            for( int i=0; i<sizePt; i++ )
            {
                if(  pt >= ptMin[i]  &&  pt < ptMax[i] ){  iPt=i; break; }
                else iPt=i;
            } 

            float  weight = nominal[iEta][iPt];
            if( shift < 0 ) weight -= sigmaNeg[iEta][iPt];
            if( shift > 0 ) weight += sigmaPos[iEta][iPt];

            return weight;
        }
        float getSF_TightMuon( Lepton muon, int shift=0 )
        {
            if( muon.LeptonType != 13 )
            { 
                std::cout<<">> [ERROR] LeptonSFUtil::getElectronSF only for Lepton::LeptonType==13(muon)"<<std::endl;
                return 1.;
            }
            float weight=1;
            return weight;
        }
};
#endif
