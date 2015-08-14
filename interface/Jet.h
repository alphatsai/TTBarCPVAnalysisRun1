#ifndef JET_H 
#define JET_H 

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

class Jet{
    public:
        Jet(){};
        Jet( JetInfoBranches& JetInfo, int idx ){
            Fill( JetInfo, idx );
        }

        void Fill( JetInfoBranches& JetInfo, int idx ){
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
};

#endif
