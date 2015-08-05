#ifndef REREGISTEVT_H
#define REREGISTEVT_H

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

void reRegistEvt( EvtInfoBranches& OldEvt, EvtInfoBranches& NewEvt )
{
    NewEvt.RunNo = OldEvt.RunNo;
    NewEvt.EvtNo = OldEvt.EvtNo;
    NewEvt.BxNo = OldEvt.BxNo;
    NewEvt.LumiNo = OldEvt.LumiNo;
    NewEvt.Orbit = OldEvt.Orbit;
    NewEvt.McFlag = OldEvt.McFlag;
    NewEvt.McSigTag = OldEvt.McSigTag;
    //McbprimeMode[2];
    //MctprimeMode[2];
    //McWMode[4];
    //McZMode[2];
    //McbprimeMass[2];
    //MctprimeMass[2];
    //MctopMass[2];
    //McWMass[4];
    //McZMass[2];
    //McDauPt[14];
    //McDauEta[14];
    //McDauPhi[14];
    //RhoPU[2];
    //SigmaPU[2];
    //McDauPdgID[14];
    //PDFid1;
    //PDFid2;
    //PDFx1;
    //PDFx2;
    //PDFscale;
    //PDFv1;
    //PDFv2;

    //BeamSpotX;
    //BeamSpotY;
    //BeamSpotZ;

    //nBX;
    //nPU[MAX_BX];
    //BXPU[MAX_BX];
    //TrueIT[MAX_BX];

    NewEvt.PFMET = OldEvt.PFMET;
    NewEvt.PFMETType1CorrectedPFMetUnclusteredEnUp = OldEvt.PFMETType1CorrectedPFMetUnclusteredEnUp;
    NewEvt.PFMETType1CorrectedPFMetUnclusteredEnDown = OldEvt.PFMETType1CorrectedPFMetUnclusteredEnDown;
    NewEvt.PFMETPhi = OldEvt.PFMETPhi;
    NewEvt.PFRawMET = OldEvt.PFRawMET;
    NewEvt.PFRawMETPhi = OldEvt.PFRawMETPhi;
    NewEvt.PFSumEt = OldEvt.PFSumEt;
    NewEvt.PFMETSig = OldEvt.PFMETSig;
    NewEvt.PFMETlongitudinal = OldEvt.PFMETlongitudinal;
    NewEvt.PFMETRealSig = OldEvt.PFMETRealSig;
    NewEvt.PFGenMET = OldEvt.PFGenMET;
    NewEvt.PFGenMETPhi = OldEvt.PFGenMETPhi;

    NewEvt.PFMETx = OldEvt.PFMETx;
    NewEvt.PFMETy = OldEvt.PFMETy;

    //TrgCount;
    //nTrgBook;
    //TrgBook[N_TRIGGER_BOOKINGS];
    //nHLT;
    //HLTPrescaleFactor[512];
    //HLTName2enum[512];
    //HLTbits[N_TRIGGER_BOOKINGS];
    //L1[128];
    //TT[64];
    //HighPurityFraction;
    //NofTracks;
    //ptHat;
}

#endif
