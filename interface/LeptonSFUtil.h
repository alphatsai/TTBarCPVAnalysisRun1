#ifndef LEPTONSFUTIL_H 
#define LEPTONSFUTIL_H

#include <string>
#include "Lepton.h"

class LeptonSFUtil
{
    public:
        LeptonSFUtil(){};
        ~LeptonSFUtil(){};
        float getSF_TightElectronID( Lepton electron, int shift=0 )
        {
            if( electron.LeptonType != 11 )
            { 
                std::cout<<">> [ERROR] LeptonSFUtil::getSF_TightElectronID only for Lepton::LeptonType==11(electron)"<<std::endl;
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
                else if( eta < etaMin[i] ) iEta=0;
                else iEta=i;
            } 
            for( int i=0; i<sizePt; i++ )
            {
                if(  pt >= ptMin[i]  &&  pt < ptMax[i] ){  iPt=i; break; }
                else if( pt < ptMin[i] ) iPt=0;
                else iPt=i;
            } 

            float  weight = nominal[iEta][iPt];
            if( shift < 0 ) weight -= sigmaNeg[iEta][iPt];
            if( shift > 0 ) weight += sigmaPos[iEta][iPt];

            return weight;
        }
        float getSF_TightMuonID( Lepton muon, int shift=0 )
        {
            if( muon.LeptonType != 13 )
            { 
                std::cout<<">> [ERROR] LeptonSFUtil::getSF_TightMuonID only for Lepton::LeptonType==13(muon)"<<std::endl;
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
            
            //// * DATA_over_MC_Tight_pt_*
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
                else if( eta < etaMin[i] ) iEta=0;
                else iEta=i;
            } 
            for( int i=0; i<sizePt; i++ )
            {
                if(  pt >= ptMin[i]  &&  pt < ptMax[i] ){  iPt=i; break; }
                else if( pt < ptMin[i] ) iPt=0;
                else iPt=i;
            } 

            float  weight = nominal[iEta][iPt];
            if( shift < 0 ) weight -= sigmaNeg[iEta][iPt];
            if( shift > 0 ) weight += sigmaPos[iEta][iPt];
            return weight;
        }
        float getSF_TightMuonIso( Lepton muon, int shift=0 )
        {
            if( muon.LeptonType != 13 )
            { 
                std::cout<<">> [ERROR] LeptonSFUtil::getSF_TightMuonIso only for Lepton::LeptonType==13(muon)"<<std::endl;
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

            //// * DATA_over_MC_combRelIsoPF04dBeta<012_*
            float nominal[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.940338,0.976738,0.996173,0.993228,0.993713,0.994084,0.996435,0.998658,1.000330,0.998813}, // eta 0
                {0.948434,0.986367,1.000356,1.000087,0.999092,0.998868,0.997890,0.999215,1.001421,1.001854}, // eta 1 
                {0.972437,0.990054,1.002820,1.004046,1.002147,1.000916,1.000404,1.000502,0.999158,0.996440}, // eta 2
                {1.116727,1.115540,1.096718,1.074812,1.060585,1.037707,1.025292,1.014950,1.008129,1.010620}, // eta 3
            };
            float sigmaNeg[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.004319,0.002177,0.001226,0.000798,0.000543,0.000273,0.000497,0.000578,0.001089,0.002471}, // eta 0
                {0.004799,0.003062,0.002012,0.001434,0.000924,0.000435,0.000794,0.000970,0.001708,0.004325}, // eta 1 
                {0.002315,0.001529,0.001033,0.000791,0.000562,0.000260,0.000480,0.000596,0.001189,0.002359}, // eta 2
                {0.005545,0.004065,0.002836,0.002155,0.001661,0.000918,0.001669,0.002094,0.004357,0.014975}, // eta 3
            };
            float sigmaPos[sizeEta][sizePt]={ 
            //pt        0,       1,       2,       3,       4,       5,       6,       7,       8,       9
                {0.004324,0.002182,0.001229,0.000800,0.000544,0.000274,0.000501,0.000585,0.001133,0.002807}, // eta 0
                {0.004821,0.003080,0.002021,0.001440,0.000930,0.000437,0.000805,0.000994,0.001854,0.005585}, // eta 1 
                {0.002323,0.001534,0.001037,0.000794,0.000564,0.000260,0.000485,0.000606,0.001261,0.003209}, // eta 2
                {0.005565,0.004078,0.002847,0.002164,0.001668,0.000923,0.001692,0.002146,0.004701,0.020880}, // eta 3
            };

            int iPt(-1);
            int iEta(-1);
            for( int i=0; i<sizeEta; i++ )
            {
                if( eta >= etaMin[i] && eta < etaMax[i] ){ iEta=i; break; }
                else if( eta < etaMin[i] ) iEta=0;
                else iEta=i;
            } 
            for( int i=0; i<sizePt; i++ )
            {
                if(  pt >= ptMin[i]  &&  pt < ptMax[i] ){ iPt=i; break; }
                else if( pt < ptMin[i] ) iPt=0;
                else iPt=i; 
            } 

            float  weight = nominal[iEta][iPt];
            if( shift < 0 ) weight -= sigmaNeg[iEta][iPt];
            if( shift > 0 ) weight += sigmaPos[iEta][iPt];
            return weight;
        }
};
#endif
