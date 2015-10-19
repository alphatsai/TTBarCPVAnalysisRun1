#ifndef SEMILEPTANICTREEBRANCHES_H 
#define SEMILEPTANICTREEBRANCHES_H 

#include <vector>
#include "Jet.h"
#include "Lepton.h"
#include "TopCandidate.h"

#define MAX_JETCOL 256
#define MAX_LEPCOL 256

// Jet
class newBranchesJet
{
    public:
        float Et;
        float Pt;
        float Px;
        float Py;
        float Pz;
        float Eta;
        float Phi;
        float Mass;
        float Energy;
        float bTag;
        void RegisterTree( TTree *root, std::string name="JetInfo")
        {
            root->Branch((name+".Et").c_str(),     &Et,     (name+".Et/F").c_str()     );
            root->Branch((name+".Pt").c_str(),     &Pt,     (name+".Pt/F").c_str()     );
            root->Branch((name+".Px").c_str(),     &Px,     (name+".Px/F").c_str()     );
            root->Branch((name+".Py").c_str(),     &Py,     (name+".Py/F").c_str()     );
            root->Branch((name+".Pz").c_str(),     &Pz,     (name+".Pz/F").c_str()     );
            root->Branch((name+".Eta").c_str(),    &Eta,    (name+".Eta/F").c_str()    );
            root->Branch((name+".Phi").c_str(),    &Phi,    (name+".Phi/F").c_str()    );
            root->Branch((name+".Mass").c_str(),   &Mass,   (name+".Mass/F").c_str()   );
            root->Branch((name+".Energy").c_str(), &Energy, (name+".Energy/F").c_str() );
            root->Branch((name+".bTag").c_str(),   &bTag,   (name+".bTag/F").c_str()   );
        }
        void Register( TTree *root, std::string name="JetInfo")
        {
            root->SetBranchAddress((name+".Et").c_str(),     &Et     );
            root->SetBranchAddress((name+".Pt").c_str(),     &Pt     );
            root->SetBranchAddress((name+".Px").c_str(),     &Px     );
            root->SetBranchAddress((name+".Py").c_str(),     &Py     );
            root->SetBranchAddress((name+".Pz").c_str(),     &Pz     );
            root->SetBranchAddress((name+".Eta").c_str(),    &Eta    );
            root->SetBranchAddress((name+".Phi").c_str(),    &Phi    );
            root->SetBranchAddress((name+".Mass").c_str(),   &Mass   );
            root->SetBranchAddress((name+".Energy").c_str(), &Energy );
            root->SetBranchAddress((name+".bTag").c_str(),   &bTag   );
        }

};
class newBranchesJetCol
{
    public:
        int Size;
        float Et[MAX_JETCOL];
        float Pt[MAX_JETCOL];
        float Px[MAX_JETCOL];
        float Py[MAX_JETCOL];
        float Pz[MAX_JETCOL];
        float Eta[MAX_JETCOL];
        float Phi[MAX_JETCOL];
        float Mass[MAX_JETCOL];
        float Energy[MAX_JETCOL];
        float bTag[MAX_JETCOL];
        void RegisterTree( TTree *root, std::string name="JetColInfo")
        {
            root->Branch((name+".Size").c_str(),   &Size,      (name+".Size/I").c_str()                  );
            root->Branch((name+".Et").c_str(),     &Et[0],     (name+".Et["+name+".Size]/F").c_str()     );
            root->Branch((name+".Pt").c_str(),     &Pt[0],     (name+".Pt["+name+".Size]/F").c_str()     );
            root->Branch((name+".Px").c_str(),     &Px[0],     (name+".Px["+name+".Size]/F").c_str()     );
            root->Branch((name+".Py").c_str(),     &Py[0],     (name+".Py["+name+".Size]/F").c_str()     );
            root->Branch((name+".Pz").c_str(),     &Pz[0],     (name+".Pz["+name+".Size]/F").c_str()     );
            root->Branch((name+".Eta").c_str(),    &Eta[0],    (name+".Eta["+name+".Size]/F").c_str()    );
            root->Branch((name+".Phi").c_str(),    &Phi[0],    (name+".Phi["+name+".Size]/F").c_str()    );
            root->Branch((name+".Mass").c_str(),   &Mass[0],   (name+".Mass["+name+".Size]/F").c_str()   );
            root->Branch((name+".Energy").c_str(), &Energy[0], (name+".Energy["+name+".Size]/F").c_str() );
            root->Branch((name+".bTag").c_str(),   &bTag[0],   (name+".bTag["+name+".Size]/F").c_str()   );
        }
        void Register( TTree *root, std::string name="JetColInfo")
        {
            root->SetBranchAddress((name+".Size").c_str(),   &Size      );
            root->SetBranchAddress((name+".Et").c_str(),     &Et[0]     );
            root->SetBranchAddress((name+".Pt").c_str(),     &Pt[0]     );
            root->SetBranchAddress((name+".Px").c_str(),     &Px[0]     );
            root->SetBranchAddress((name+".Py").c_str(),     &Py[0]     );
            root->SetBranchAddress((name+".Pz").c_str(),     &Pz[0]     );
            root->SetBranchAddress((name+".Eta").c_str(),    &Eta[0]    );
            root->SetBranchAddress((name+".Phi").c_str(),    &Phi[0]    );
            root->SetBranchAddress((name+".Mass").c_str(),   &Mass[0]   );
            root->SetBranchAddress((name+".Energy").c_str(), &Energy[0] );
            root->SetBranchAddress((name+".bTag").c_str(),   &bTag[0]   );
        }

};

// Lepton
class newBranchesLep
{
    public:
        float Et;
        float Pt;
        float Px;
        float Py;
        float Pz;
        float Eta;
        float Phi;
        int Charge;
        int LeptonType;
        void RegisterTree( TTree *root, std::string name="LepInfo")
        {
            root->Branch((name+".Et").c_str(),         &Et,         (name+".Et/F").c_str()         );
            root->Branch((name+".Pt").c_str(),         &Pt,         (name+".Pt/F").c_str()         );
            root->Branch((name+".Px").c_str(),         &Px,         (name+".Px/F").c_str()         );
            root->Branch((name+".Py").c_str(),         &Py,         (name+".Py/F").c_str()         );
            root->Branch((name+".Pz").c_str(),         &Pz,         (name+".Pz/F").c_str()         );
            root->Branch((name+".Eta").c_str(),        &Eta,        (name+".Eta/F").c_str()        );
            root->Branch((name+".Phi").c_str(),        &Phi,        (name+".Phi/F").c_str()        );
            root->Branch((name+".Charge").c_str(),     &Charge,     (name+".Charge/I").c_str()     );
            root->Branch((name+".LeptonType").c_str(), &LeptonType, (name+".LeptonType/I").c_str() );
        }
        void Register( TTree *root, std::string name="LepInfo")
        {
            root->SetBranchAddress((name+".Et").c_str(),         &Et         );
            root->SetBranchAddress((name+".Pt").c_str(),         &Pt         );
            root->SetBranchAddress((name+".Px").c_str(),         &Px         );
            root->SetBranchAddress((name+".Py").c_str(),         &Py         );
            root->SetBranchAddress((name+".Pz").c_str(),         &Pz         );
            root->SetBranchAddress((name+".Eta").c_str(),        &Eta        );
            root->SetBranchAddress((name+".Phi").c_str(),        &Phi        );
            root->SetBranchAddress((name+".Charge").c_str(),     &Charge     );
            root->SetBranchAddress((name+".LeptonType").c_str(), &LeptonType );
        }
};
class newBranchesLepCol
{
    public:
        int Size;
        float Et[MAX_LEPCOL];
        float Pt[MAX_LEPCOL];
        float Px[MAX_LEPCOL];
        float Py[MAX_LEPCOL];
        float Pz[MAX_LEPCOL];
        float Eta[MAX_LEPCOL];
        float Phi[MAX_LEPCOL];
        int Charge[MAX_LEPCOL];
        int LeptonType[MAX_LEPCOL];
        void RegisterTree( TTree *root, std::string name="LepColInfo")
        {
            root->Branch((name+".Size").c_str(),       &Size,          (name+".Size/I").c_str()                      );
            root->Branch((name+".Et").c_str(),         &Et[0],         (name+".Et["+name+".Size]/F").c_str()         );
            root->Branch((name+".Pt").c_str(),         &Pt[0],         (name+".Pt["+name+".Size]/F").c_str()         );
            root->Branch((name+".Px").c_str(),         &Px[0],         (name+".Px["+name+".Size]/F").c_str()         );
            root->Branch((name+".Py").c_str(),         &Py[0],         (name+".Py["+name+".Size]/F").c_str()         );
            root->Branch((name+".Pz").c_str(),         &Pz[0],         (name+".Pz["+name+".Size]/F").c_str()         );
            root->Branch((name+".Eta").c_str(),        &Eta[0],        (name+".Eta["+name+".Size]/F").c_str()        );
            root->Branch((name+".Phi").c_str(),        &Phi[0],        (name+".Phi["+name+".Size]/F").c_str()        );
            root->Branch((name+".Charge").c_str(),     &Charge[0],     (name+".Charge["+name+".Size]/I").c_str()     );
            root->Branch((name+".LeptonType").c_str(), &LeptonType[0], (name+".LeptonType["+name+".Size]/I").c_str() );
        }
        void Register( TTree *root, std::string name="LepColInfo")
        {
            root->SetBranchAddress((name+".Size").c_str(),       &Size          );
            root->SetBranchAddress((name+".Et").c_str(),         &Et[0]         );
            root->SetBranchAddress((name+".Pt").c_str(),         &Pt[0]         );
            root->SetBranchAddress((name+".Px").c_str(),         &Px[0]         );
            root->SetBranchAddress((name+".Py").c_str(),         &Py[0]         );
            root->SetBranchAddress((name+".Pz").c_str(),         &Pz[0]         );
            root->SetBranchAddress((name+".Eta").c_str(),        &Eta[0]        );
            root->SetBranchAddress((name+".Phi").c_str(),        &Phi[0]        );
            root->SetBranchAddress((name+".Charge").c_str(),     &Charge[0]     );
            root->SetBranchAddress((name+".LeptonType").c_str(), &LeptonType[0] );
        }
};

// Top
class newBranchesTopHadronic
{
    public:
        int isHadronicTop; 
        int isLeptonicTop; 
        int isAntiB;
        int idxJet1;
        int idxJet2;
        float Pt;
        float Px;
        float Py;
        float Pz;
        float Eta;
        float Phi;
        float Mass;
        float Energy;
        void RegisterTree( TTree *root, std::string name="TopInfo")
        {
            root->Branch((name+".isHadronicTop").c_str(), &isHadronicTop, (name+".isHadronicTop/I").c_str()   );
            root->Branch((name+".isLeptonicTop").c_str(), &isLeptonicTop, (name+".isLeptonicTop/I").c_str()   );
            root->Branch((name+".isAntiB").c_str(),       &isAntiB,       (name+".isAntiB/I").c_str()         );
            root->Branch((name+".idxJet1").c_str(),       &idxJet1,       (name+".idxJet1/I").c_str()         );
            root->Branch((name+".idxJet2").c_str(),       &idxJet2,       (name+".idxJet2/I").c_str()         );
            root->Branch((name+".Pt").c_str(),            &Pt,            (name+".Pt/F").c_str()              );
            root->Branch((name+".Px").c_str(),            &Px,            (name+".Px/F").c_str()              );
            root->Branch((name+".Py").c_str(),            &Py,            (name+".Py/F").c_str()              );
            root->Branch((name+".Pz").c_str(),            &Pz,            (name+".Pz/F").c_str()              );
            root->Branch((name+".Eta").c_str(),           &Eta,           (name+".Eta/F").c_str()             );
            root->Branch((name+".Phi").c_str(),           &Phi,           (name+".Phi/F").c_str()             );
            root->Branch((name+".Mass").c_str(),          &Mass,          (name+".Mass/F").c_str()            );
            root->Branch((name+".Energy").c_str(),        &Energy,        (name+".Energy/F").c_str()          );
        }
        void Register( TTree *root, std::string name="TopInfo")
        {
            root->SetBranchAddress((name+".isHadronicTop").c_str(), &isHadronicTop );
            root->SetBranchAddress((name+".isLeptonicTop").c_str(), &isLeptonicTop );
            root->SetBranchAddress((name+".isAntiB").c_str(),       &isAntiB       );
            root->SetBranchAddress((name+".idxJet1").c_str(),       &idxJet1       );
            root->SetBranchAddress((name+".idxJet2").c_str(),       &idxJet2       );
            root->SetBranchAddress((name+".Pt").c_str(),            &Pt            );
            root->SetBranchAddress((name+".Px").c_str(),            &Px            );
            root->SetBranchAddress((name+".Py").c_str(),            &Py            );
            root->SetBranchAddress((name+".Pz").c_str(),            &Pz            );
            root->SetBranchAddress((name+".Eta").c_str(),           &Eta           );
            root->SetBranchAddress((name+".Phi").c_str(),           &Phi           );
            root->SetBranchAddress((name+".Mass").c_str(),          &Mass          );
            root->SetBranchAddress((name+".Energy").c_str(),        &Energy        );
        }
};

class SemiLeptanicTreeBranches
{
    public:
        newBranchesJet    BJetNewBranches;  
        newBranchesJet    BbarJetNewBranches;  
        newBranchesJetCol nonBJetColNewBranches; 
        newBranchesLep    isoLepNewBranches;         
        newBranchesTopHadronic topHadronicNewBranches;

        // Tree builder
        void RegisterTree( TTree *root )
        {
            BJetNewBranches.RegisterTree(        root, "bJet"        );
            BbarJetNewBranches.RegisterTree(     root, "bbarJet"     );
            nonBJetColNewBranches.RegisterTree(  root, "nonBJetCol"  );
            isoLepNewBranches.RegisterTree(      root, "isoLepton"   );
            topHadronicNewBranches.RegisterTree( root, "topHadronic" );
        } 
        void Register( TTree *root )
        {
            BJetNewBranches.Register(        root, "bJet"        );
            BbarJetNewBranches.Register(     root, "bbarJet"     );
            nonBJetColNewBranches.Register(  root, "nonBJetCol"  );
            isoLepNewBranches.Register(      root, "isoLepton"   );
            topHadronicNewBranches.Register( root, "topHadronic" );
        }

        // Fill branches
        void fill_BJetNewBranches( Jet bjet )
        {
            BJetNewBranches.Et   = bjet.Et;
            BJetNewBranches.Pt   = bjet.Pt;
            BJetNewBranches.Px   = bjet.Px;
            BJetNewBranches.Py   = bjet.Py;
            BJetNewBranches.Pz   = bjet.Pz;
            BJetNewBranches.Eta  = bjet.Eta;
            BJetNewBranches.Phi  = bjet.Phi;
            BJetNewBranches.Mass = bjet.Mass;
            BJetNewBranches.Energy = bjet.Energy;
            BJetNewBranches.bTag   = bjet.CombinedSVBJetTags;
        }
        void fill_BbarJetNewBranches( Jet bbarjet )
        {
            BbarJetNewBranches.Et   = bbarjet.Et;
            BbarJetNewBranches.Pt   = bbarjet.Pt;
            BbarJetNewBranches.Px   = bbarjet.Px;
            BbarJetNewBranches.Py   = bbarjet.Py;
            BbarJetNewBranches.Pz   = bbarjet.Pz;
            BbarJetNewBranches.Eta  = bbarjet.Eta;
            BbarJetNewBranches.Phi  = bbarjet.Phi;
            BbarJetNewBranches.Mass = bbarjet.Mass;
            BbarJetNewBranches.Energy = bbarjet.Energy;
            BbarJetNewBranches.bTag   = bbarjet.CombinedSVBJetTags;
        }
        void fill_nonBJetColNewBranches( vector<Jet> nonBJetCol )
        {
            const int size=nonBJetCol.size();
            for( int i=0; i<size; i++ )
            {
                nonBJetColNewBranches.Et[i]   = nonBJetCol[i].Et;
                nonBJetColNewBranches.Pt[i]   = nonBJetCol[i].Pt;
                nonBJetColNewBranches.Px[i]   = nonBJetCol[i].Px;
                nonBJetColNewBranches.Py[i]   = nonBJetCol[i].Py;
                nonBJetColNewBranches.Pz[i]   = nonBJetCol[i].Pz;
                nonBJetColNewBranches.Eta[i]  = nonBJetCol[i].Eta;
                nonBJetColNewBranches.Phi[i]  = nonBJetCol[i].Phi;
                nonBJetColNewBranches.Mass[i] = nonBJetCol[i].Mass;
                nonBJetColNewBranches.Energy[i] = nonBJetCol[i].Energy;
                nonBJetColNewBranches.bTag[i]   = nonBJetCol[i].CombinedSVBJetTags;
            }
        } 
        void fill_isoLepNewBranches( Lepton lep )
        {
            isoLepNewBranches.Et = lep.Et;
            isoLepNewBranches.Pt = lep.Pt;
            isoLepNewBranches.Px = lep.Px;
            isoLepNewBranches.Py = lep.Py;
            isoLepNewBranches.Pz = lep.Pz;
            isoLepNewBranches.Eta = lep.Eta;
            isoLepNewBranches.Phi = lep.Phi;
            isoLepNewBranches.Charge = lep.Charge;
            isoLepNewBranches.LeptonType = lep.LeptonType;
        }
        void fill_topHadronicNewBranches( TopCandidate topHadronic, int idxJet1, int idxJet2 ) 
        {
            if( !topHadronic.isHadronicTop ) std::cout<<">> [WARNING] SemiLeptanicTreeBranches::fill_topHadronicNewBranches input with not hardronic top."<<std::endl;
            topHadronicNewBranches.isHadronicTop = topHadronic.isHadronicTop; 
            topHadronicNewBranches.isLeptonicTop = topHadronic.isLeptonicTop;
            topHadronicNewBranches.isAntiB = topHadronic.isAntiB;
            topHadronicNewBranches.idxJet1 = idxJet1;
            topHadronicNewBranches.idxJet2 = idxJet2;
            topHadronicNewBranches.Pt  = topHadronic.P4.Pt();
            topHadronicNewBranches.Px  = topHadronic.P4.Px();
            topHadronicNewBranches.Py  = topHadronic.P4.Py();
            topHadronicNewBranches.Pz  = topHadronic.P4.Pz();
            topHadronicNewBranches.Eta = topHadronic.P4.Eta();
            topHadronicNewBranches.Phi = topHadronic.P4.Phi();
            topHadronicNewBranches.Mass   = topHadronic.P4.M();
            topHadronicNewBranches.Energy = topHadronic.P4.Energy();
        } 
};

#endif 
