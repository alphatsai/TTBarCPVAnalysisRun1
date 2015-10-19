#ifndef JET_H 
#define JET_H 

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"
using namespace std;

class Jet
{
    public:
        // Constructor
        Jet()
        {
            appliedJER=false; shiftJER=0;
            appliedJES=false; shiftJES=0;
            filledInfo=false;
        }
        Jet( JetInfoBranches& JetInfo, int idx ){
            appliedJER=false; shiftJER=0;
            appliedJER=false; shiftJES=0;
            filledInfo=false;
            Fill( JetInfo, idx );
        }

        // Apply systematic unc and scale factors
        int applyJER( int shift=0 )
        {
            if( !filledInfo )
            {
                std::cout<<">> [ERROR] Please do Jet::fill(JetInfoBranches& JetCol, int index) first!"<<std::endl;
                return -1000;
            }

            appliedJER=true; shiftJER=shift;

            // Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
            int iEta=-1;
            float etaMin[]={0.0,0.5,1.1,1.7,2.3,2.8,3.2};
            float etaMax[]={0.5,1.1,1.7,2.3,2.8,3.2,5.0};
            float nominal[]={1.079,1.099,1.121,1.208,1.254,1.395,1.056}; 
            float sigmaStat[]={0.005,0.005,0.005,0.013,0.026,0.036,0.048}; 
            float sigmaSyst[]={0.026,0.028,0.029,0.045,0.056,0.051,0.185};

            // Decide eta region
            for( int i=0; i<int(sizeof(etaMin)/sizeof(etaMin[0])); i++)
            {
                if( fabs(Eta) >= etaMin[i] && fabs(Eta) < etaMax[i] ){ iEta=i; break; }
                else iEta=i;

            }

            // JER secale factor
            float sf = nominal[iEta];
            if( shiftJER > 0. ) sf += sqrt( sigmaStat[iEta]*sigmaStat[iEta] + sigmaSyst[iEta]*sigmaSyst[iEta] );
            if( shiftJER < 0. ) sf -= sqrt( sigmaStat[iEta]*sigmaStat[iEta] + sigmaSyst[iEta]*sigmaSyst[iEta] );

            //// Debug
            //std::cout<<">> fabs(eta): "<<fabs(Eta)<<", Eta: "<<Eta<<" , index: "<<iEta<<std::endl;
            //std::cout<<"       geneta: "<<fabs(GenJetEta)<<" , genPt: "<<GenJetPt<<", genphi "<<GenJetPhi<<std::endl;
            //std::cout<<"          Old pt: "<<Pt<<" ,Et: "<<Et<<std::endl;

            // Caclute pT relation between matched Gen-jet 
            float rescale=1.;
            if( GenJetEta != 0 && GenJetPt != 0 && GenJetPhi != 0 ) // Has matched GenJet
            {
                float caculatedPt = GenJetPt + ( Pt - GenJetPt )*sf;
                float maxPt = std::max( float(0.), caculatedPt );
                if( maxPt > 0. ) rescale = caculatedPt/Pt; 
            }

            Pt = Pt * rescale;
            Et = Et * rescale;
            Px = Px * rescale; 
            Py = Py * rescale; 
            Pz = Pz * rescale;
            Energy = Energy * rescale; 
            PtCorrRaw   = PtCorrRaw   * rescale;
            PtCorrL3    = PtCorrL3    * rescale;
            PtCorrL7g   = PtCorrL7g   * rescale;
            PtCorrL7uds = PtCorrL7uds * rescale;
            PtCorrL7c   = PtCorrL7c   * rescale;
            PtCorrL7b   = PtCorrL7b   * rescale;

            P3.SetXYZ( Px, Py, Pz );
            P4.SetPxPyPzE( Px, Py, Pz, Energy );

            Mass = P4.M();
            //std::cout<<"          new pt: "<<P4.Pt()<<" ,Et: "<<P4.Et()<<std::endl;

            return shiftJER;

        }

        int applyJES( int shift=0 )
        {
            appliedJES=true; shiftJES=shift;
            return shiftJES;
        }

        // Fill all infomations from bprimeKits 
        void Fill( JetInfoBranches& JetInfo, int idx )
        {
            filledInfo=true;
            Index = JetInfo.Index[idx];
            NTracks = JetInfo.NTracks[idx];
            Et = JetInfo.Et[idx];
            Pt = JetInfo.Pt[idx];
            Unc = JetInfo.Unc[idx];
            Eta = JetInfo.Eta[idx];
            Phi = JetInfo.Phi[idx];
            JetIDLOOSE = JetInfo.JetIDLOOSE[idx]; 
            JetCharge = JetInfo.JetCharge[idx];
            QGTagsMLP = JetInfo.QGTagsMLP[idx];
            QGTagsLikelihood = JetInfo.QGTagsLikelihood[idx];
            NConstituents = JetInfo.NConstituents[idx];
            NCH = JetInfo.NCH[idx];
            CEF = JetInfo.CEF[idx];
            NHF = JetInfo.NHF[idx];
            NEF = JetInfo.NEF[idx];
            CHF = JetInfo.CHF[idx];
            JVAlpha = JetInfo.JVAlpha[idx];
            JVBeta = JetInfo.JVBeta[idx];
            PtCorrRaw = JetInfo.PtCorrRaw[idx];  
            PtCorrL2 = JetInfo.PtCorrL2[idx];  
            PtCorrL3 = JetInfo.PtCorrL3[idx];  
            PtCorrL7g = JetInfo.PtCorrL7g[idx];
            PtCorrL7uds = JetInfo.PtCorrL7uds[idx];
            PtCorrL7c = JetInfo.PtCorrL7c[idx];  
            PtCorrL7b = JetInfo.PtCorrL7b[idx];  
            JetBProbBJetTags = JetInfo.JetBProbBJetTags[idx];
            JetProbBJetTags = JetInfo.JetProbBJetTags[idx];
            TrackCountHiPurBJetTags = JetInfo.TrackCountHiPurBJetTags[idx];  
            TrackCountHiEffBJetTags = JetInfo.TrackCountHiEffBJetTags[idx]; 

            SimpleSVBJetTags = JetInfo.SimpleSVBJetTags[idx];  
            SimpleSVHEBJetTags = JetInfo.SimpleSVHEBJetTags[idx];  
            SimpleSVHPBJetTags = JetInfo.SimpleSVHPBJetTags[idx];  
            CombinedSVBJetTags = JetInfo.CombinedSVBJetTags[idx];
            CombinedSVMVABJetTags = JetInfo.CombinedSVMVABJetTags[idx];
            SoftElecByIP3dBJetTags = JetInfo.SoftElecByIP3dBJetTags[idx];
            SoftElecByPtBJetTags = JetInfo.SoftElecByPtBJetTags[idx];  
            SoftMuonBJetTags = JetInfo.SoftMuonBJetTags[idx];      
            SoftMuonByIP3dBJetTags = JetInfo.SoftMuonByIP3dBJetTags[idx];
            SoftMuonByPtBJetTags = JetInfo.SoftMuonByPtBJetTags[idx];  
            DoubleSVHighEffBJetTags = JetInfo.DoubleSVHighEffBJetTags[idx]; 

            GenJetPt = JetInfo.GenJetPt[idx];
            GenJetEta = JetInfo.GenJetEta[idx];
            GenJetPhi = JetInfo.GenJetPhi[idx];
            GenPt = JetInfo.GenPt[idx];
            GenEta = JetInfo.GenEta[idx];
            GenPhi = JetInfo.GenPhi[idx];
            GenPdgID = JetInfo.GenPdgID[idx];
            GenFlavor = JetInfo.GenFlavor[idx];
            GenMCTag = JetInfo.GenMCTag[idx]; // 0: unknown, 1: decay from W, 2: decay from Z, (+10) from b', (+20) from t'

            Px = JetInfo.Px[idx]; 
            Py = JetInfo.Py[idx]; 
            Pz = JetInfo.Pz[idx]; 
            Energy = JetInfo.Energy[idx]; 

            Mass = JetInfo.Mass[idx];
            Area = JetInfo.Area[idx];

            NSubjets = JetInfo.NSubjets[idx];
            SubjetsIdxStart = JetInfo.SubjetsIdxStart[idx];

            P3.SetXYZ( Px, Py, Pz );
            P4.SetPxPyPzE( Px, Py, Pz, Energy );
        }

        TVector3 P3;
        TLorentzVector P4;

        int   Index;
        int   NTracks;
        float Et;
        float Pt;
        float Unc;
        float Eta;
        float Phi;
        int   JetIDLOOSE; 
        float JetCharge;
        float QGTagsMLP;
        float QGTagsLikelihood;
        int   NConstituents;
        int   NCH;
        float CEF;
        float NHF;
        float NEF;
        float CHF;
        float JVAlpha;
        float JVBeta;
        float PtCorrRaw;  
        float PtCorrL2;  
        float PtCorrL3;  
        float PtCorrL7g;
        float PtCorrL7uds;
        float PtCorrL7c;  
        float PtCorrL7b;  
        float JetBProbBJetTags;
        float JetProbBJetTags;
        float TrackCountHiPurBJetTags;  
        float TrackCountHiEffBJetTags; 

        float SimpleSVBJetTags;  
        float SimpleSVHEBJetTags;  
        float SimpleSVHPBJetTags;  
        float CombinedSVBJetTags;
        float CombinedSVMVABJetTags;
        float SoftElecByIP3dBJetTags;
        float SoftElecByPtBJetTags;  
        float SoftMuonBJetTags;      
        float SoftMuonByIP3dBJetTags;
        float SoftMuonByPtBJetTags;  
        float DoubleSVHighEffBJetTags; 

        float GenJetPt;
        float GenJetEta;
        float GenJetPhi;
        float GenPt;
        float GenEta;
        float GenPhi;
        int   GenPdgID;
        int   GenFlavor;
        int    GenMCTag; // 0: unknown, 1: decay from W, 2: decay from Z, (+10) from b', (+20) from t'

        float Px; 
        float Py; 
        float Pz; 
        float Energy; 

        float Mass;
        float Area;

        int NSubjets;
        int SubjetsIdxStart;

        std::vector<float> *SubjetMass;
        std::vector<float> *SubjetPt;
        std::vector<float> *SubjetEt;
        std::vector<float> *SubjetEta;
        std::vector<float> *SubjetPhi;
        std::vector<float> *SubjetCombinedSVBJetTags;
        std::vector<float> *SubjetPtUncorr;
        std::vector<float> *SubjetArea;

        bool appliedJER;
        bool appliedJES;
        int  shiftJER;
        int  shiftJES;

        private:
            bool filledInfo;
};
#endif
