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

            float pt = muon.Pt;
            float eta=fabs(muon.Eta);
            float ptMin[]={10,20,25,30,35,40,50,60,90,140};
            float ptMax[]={20,25,30,35,40,50,60,90,140,300};
            float etaMin[]={0.0,0.9,1.2,2.1};
            float etaMax[]={0.9,1.2,2.1,2.4};
            const int sizePt  = int(sizeof(ptMin)/sizeof(ptMin[0]));
            const int sizeEta = int(sizeof(etaMin)/sizeof(etaMin[0]));

            float nominal[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.970275,0.988865,0.992338,0.993283,0.993662,0.992392,0.991190,0.989417,1.003749,1.018503}, // eta 0
                {1.001731,0.993947,0.994722,0.993391,0.992285,0.991870,0.995010,0.990406,1.009028,1.010956}, // eta 1 
                {1.018018,1.000351,0.998486,0.996558,0.996026,0.996062,0.995183,0.992486,1.023129,0.974754}, // eta 2
                {1.005044,0.998089,0.996183,1.000551,0.992563,0.995144,0.993590,0.989484,1.060173,0.890547}, // eta 3
            };
            float sigmaNeg[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.004997,0.001638,0.000770,0.000535,0.000412,0.000266,0.000652,0.001032,0.003136,0.017382}, // eta 0
                {0.006844,0.002493,0.001378,0.001036,0.000764,0.000477,0.001225,0.001941,0.006195,0.032074}, // eta 1 
                {0.003795,0.001414,0.000856,0.000704,0.000585,0.000239,0.000922,0.001545,0.005351,0.030099}, // eta 2
                {0.007015,0.002920,0.001790,0.001504,0.001329,0.000982,0.002640,0.004972,0.010352,0.124062}, // eta 3
            };
            float sigmaPos[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.004982,0.001643,0.000774,0.000537,0.000415,0.000266,0.000655,0.001038,0.003164,0.017259}, // eta 0
                {0.006812,0.002504,0.001388,0.001044,0.000769,0.000479,0.001235,0.001960,0.006282,0.034515}, // eta 1 
                {0.003787,0.001418,0.000859,0.000707,0.000588,0.000239,0.000926,0.001553,0.005381,0.029173}, // eta 2
                {0.006984,0.002928,0.001799,0.001513,0.001336,0.000986,0.002655,0.004989,0.015680,0.160031}, // eta 3
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
};
#endif
