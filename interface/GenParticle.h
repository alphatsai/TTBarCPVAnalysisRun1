#ifndef GENPARTICLE_H 
#define GENPARTICLE_H 

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

class GenParticle{
    public:
        GenParticle(){}
        GenParticle( GenInfoBranches& GenInfo, int idx ){
            Fill( GenInfo, idx );
        }
        void Fill( GenInfoBranches& GenInfo, int idx ){
            Pt   = GenInfo.Pt[idx];
            Eta  = GenInfo.Eta[idx];
            Phi  = GenInfo.Phi[idx];
            Mass = GenInfo.Mass[idx];
            nMo  = GenInfo.nMo[idx];
            nDa  = GenInfo.nDa[idx];
            Mo1  = GenInfo.Mo1[idx];
            Mo2  = GenInfo.Mo2[idx];
            Da1  = GenInfo.Da1[idx];
            Da2  = GenInfo.Da2[idx];

            PdgID      = GenInfo.PdgID[idx];
            Status     = GenInfo.Status[idx];
            PhotonFlag = GenInfo.PhotonFlag[idx];   // -1 : unknown or not photon, 0 : prompt photon, 1 : decay in flight, 2 : ISR, 3 : FSR

            Mo1PdgID  = GenInfo.Mo1PdgID[idx];
            Mo2PdgID  = GenInfo.Mo2PdgID[idx];
            Mo1Status = GenInfo.Mo1Status[idx];
            Mo2Status = GenInfo.Mo2Status[idx];
            Da1PdgID  = GenInfo.Da1PdgID[idx];
            Da2PdgID  = GenInfo.Da2PdgID[idx];
            GrandMo1PdgID  = GenInfo.GrandMo1PdgID[idx];
            GrandMo2PdgID  = GenInfo.GrandMo2PdgID[idx];
            GrandMo1Status = GenInfo.GrandMo1Status[idx];
            GrandMo2Status = GenInfo.GrandMo2Status[idx];

            P3.SetPtEtaPhi(  Pt, Eta, Phi );
            P4.SetPtEtaPhiM( Pt, Eta, Phi, Mass );
        }

        TVector3 P3;
        TLorentzVector P4;

        float Pt;
        float Eta;
        float Phi;
        float Mass;

        int PdgID;
        int PhotonFlag;   // -1 : unknown or not photon, 0 : prompt photon, 1 : decay in flight, 2 : ISR, 3 : FSR
        int Status;
        int nMo;
        int nDa;
        int Mo1;
        int Mo2;
        int Da1;
        int Da2;

        int Mo1PdgID;
        int Mo2PdgID;
        int Mo1Status;
        int Mo2Status;
        int Da1PdgID;
        int Da2PdgID;
        int GrandMo1PdgID;
        int GrandMo2PdgID;
        int GrandMo1Status;
        int GrandMo2Status;
};
#endif
