#ifndef TOPCANDIDATE_H 
#define TOPCANDIDATE_H 

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h"

class TopCandidate{
    public:
        TopCandidate(){};
        TopCandidate( Jet& bjet, Jet& jet1, Jet& jet2 ){
            Fill( bjet, jet1, jet2 );
        }
        TopCandidate( Jet& bjet, Lepton& lep, float& met, float& metPhi ){
            Fill( bjet, lep, met, metPhi );
        }

        // Hadronic top: t->bW( W->qq ) : jet1 and jet2 are from W boson decay
        void Fill( Jet& bjet, Jet& jet1, Jet& jet2, int isAntib=0 ) 
        { 
            isAntiB = isAntib;

            isHadronicTop = true;
            isLeptonicTop = false;
    
            W_P4   = jet1.P4 + jet2.P4;
            P4     = bjet.P4 + W_P4;
            P3     = P4.Vect(); 
            Pt     = P4.Pt(); 
            Eta    = P4.Eta(); 
            Phi    = P4.Phi(); 
            Mass   = P4.M(); 
            MassT  = P4.Mt(); 
            Energy = P4.Energy();

            bJet = &bjet; 
            Jet1 = &jet1; 
            Jet2 = &jet2;

            MET    = -100000;
            METPhi = -100000;
        }

        // Leptonic top: t->bW( W->ln ) : met is missing Et from nuetrino NOTE: ONLY result to transverse phase
        void Fill( Jet& bjet, Lepton& lep, float& met, float& metPhi, int isAntib=0 )
        {
            isAntiB = isAntib;

            isHadronicTop = false;
            isLeptonicTop = true;
            
            TLorentzVector nu_P4;
            nu_P4.SetPtEtaPhiM( met, lep.Eta, metPhi, 0. );

            W_P4   =  lep.P4 + nu_P4;
            P4     = bjet.P4 + W_P4;
            P3     = P4.Vect(); 
            Pt     = P4.Pt(); 
            Eta    = P4.Eta(); 
            Phi    = P4.Phi(); 
            Mass   = P4.M(); 
            //MassT  = P4.Mt(); 
            MassT  = TMath::Sqrt( lep.Et*lep.Et + met*met - 2*lep.Et*met*TMath::Cos(lep.Phi-metPhi));
            Energy = P4.Energy();

            Lept = &lep;
            MET    = met;
            METPhi = metPhi;
        }
 
        TVector3 P3;
        TLorentzVector P4;
        TLorentzVector W_P4;

        bool isAntiB;
        bool isHadronicTop;
        bool isLeptonicTop;

        float Pt;
        float Eta;
        float Phi;

        float Mass;
        float MassT;
        float Energy;

        float MET;
        float METPhi;

        Jet* bJet;
        Jet* Jet1; 
        Jet* Jet2;
        Lepton* Lept; 
};

#endif
