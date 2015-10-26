#ifndef FUNCTIONS_H 
#define FUNCTIONS_H

#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Jet.h" 
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/Lepton.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/GenParticle.h" 

const float mass_t  = 172.5;
const float mass_W  = 82.9;
const float width_t = 16.3;
const float width_W = 9.5;

namespace{

    using namespace std;

    ////* Functions lists -----
    template<class valuetype> 
        std::string num2str( valuetype i );

    template<class TH1>       
        void fillAsym( TH1* h, double value, double wrt=1 );

    template<class Object>    
        void getHighPtObject(  vector<Object> col, Object &obj );

    template<class Object>    
        void get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 );

    template <class Object, class matchingObject>
        bool matchObject( Object &obj, matchingObject &mobj, vector<matchingObject> col, double dR=0.5 );

    template <class Object, class matchingObject> 
        bool matchMultiObject( vector<Object> incol, vector<matchingObject> mcol, vector<matchingObject> &outcol );

    bool checkMo( GenParticle     GenInfo, int moPdgID );
    bool checkDa( GenParticle     GenInfo, int daPdgID );
    bool checkMo( GenInfoBranches GenInfo, int idx, int moPdgID );
    bool checkDa( GenInfoBranches GenInfo, int idx, int daPdgID );
    bool isIsoLeptonFromJets( Lepton lepton, vector<Jet> jetCol, double dR=0.5 );
    void rmElelectronOverlapeMuon( vector<Lepton> &elColl, vector<Lepton> muColl, double dR=0.3 );

    float getChi2( Jet jet1, Jet jet2, Jet bjet, float M_top=mass_t, float Wth_top=width_t, float M_W=mass_W, float Wth_W=width_W );

    double Obs2( TVector3 isoLep, TVector3 hardJet, TVector3 b,   TVector3 bbar );
    double Obs4( TVector3 isoLep, TVector3 hardJet, TVector3 b,   TVector3 bbar, int charge );
    double Obs7( TVector3 beam,   TVector3 b,       TVector3 bbar               );
    double Obs3( TLorentzVector isoLep, TLorentzVector hardJet, TLorentzVector b, TLorentzVector bbar, int charge ); // need boost to bb~ CM

    ////* Defination of Functions ----

    //* change number to string
    template<class valuetype>
        std::string num2str( valuetype i )
        {
            std::string s;
            stringstream ss(s);
            ss << i;
            return ss.str();
        }

    //* fill asymetric to histogram with 2 bins
    template<class TH1>
        void fillAsym( TH1* h, double value, double wrt )
        {
            if( value < 0. ) h->Fill(-1., wrt);
            else h->Fill(0., wrt);
        }

    //* get high pT object from collection
    template <class Object>
        void getHighPtObject( vector<Object> col, Object &obj )
        {
            int o1(-1);
            double pt(0);
            const int size=col.size();
            for( int i=0; i<size; i++)
            {
                if( pt < col[i].Pt )
                {
                    pt=col[i].Pt;
                    o1=i;
                }
            }
            obj=col[o1];
            obj.Index = o1;
        }

    //* get 2 high pT object from collection
    template <class Object>
        void get2HighPtObject( vector<Object> col, Object &obj1, Object &obj2 )
        {
            int o1(-1); 
            int o2(-1);
            double pt1(0); 
            double pt2(0);
            const int size=col.size();
            for( int i=0; i<size; i++)
            {
                if( pt1 < col[i].Pt ){
                    pt2 = pt1;
                    pt1 = col[i].Pt;
                    o2  = o1;
                    o1  = i;
                }else if( pt2 < col[i].Pt ){
                    pt2 = col[i].Pt;
                    o2  = i;
                }
            }
            obj1=col[o1]; obj1.Index = o1;
            obj2=col[o2]; obj2.Index = o2;
        }

    //* Matching object with dR 
    template <class Object, class matchingObject>
        bool matchObject( Object &obj, matchingObject &mobj, vector<matchingObject> col, double dR )
        {
            int     obsize = col.size();
            int       midx = -1;
            bool   matched = false;
            double      dr = 10E+6;
            for( int i=0; i<obsize; i++ )
            {
                double dr_ = col[i].P4.DeltaR(obj.P4); 
                if( dr_ > dR ) continue;
                if( dr_ < dr )
                {
                    dr = dr_;
                    midx = i; 
                    matched = true;
                }
            }
            if( matched ) mobj = col[midx];
            else std::cout<<"[WARING] Can't find matched particle in matchObject<"<<typeid(obj).name()<<"> within dR "<<dR<<std::endl;
            return matched;
        }

    //* Complex sorting method for matching 2 object collections 
    template <class Object, class matchingObject>
        bool matchMultiObject( vector<Object> incol, vector<matchingObject> mcol, vector<matchingObject> &outcol )
        {
            if( incol.size() > mcol.size() )
            {
                std::cout<<">> [ERROR] input size: "<<incol.size()<<" > matching size: "<<mcol.size()<<std::endl;
                return false;
            }
            const int inSize = incol.size();
            const int mSize  = mcol.size();
            double   dR[inSize][mSize];
            int    iObj[inSize][mSize];
            // Sort the dR for each objects
            for( int i=0; i<inSize; i++ ){
                for( int j=0; j<mSize; j++){
                    dR[i][j]=incol[i].P4.DeltaR( mcol[j].P4 );
                    iObj[i][j]=j;
                }
                for( int j=0; j<mSize; j++){
                    for( int k=0; k<mSize; k++){
                        if( dR[i][j] < dR[i][k] ) // NOTE: j is after k ( sort number small -> big )
                        {
                            double tmp1 =   dR[i][j];
                            int    tmp2 = iObj[i][j];
                            dR[i][j] = dR[i][k];
                            dR[i][k] = tmp1;
                            iObj[i][j] = iObj[i][k];
                            iObj[i][k] = tmp2;
                        }
                    }
                }
            }
            bool matched=false;
            int iMatchedObj[inSize];
            int iMatchedIdx[inSize];
            // Pick up the smallest dR objects for each input object
            for( int i=0; i<inSize; i++ )
            { 
                iMatchedObj[i] = iObj[i][0];
                iMatchedIdx[i] = 0; 
            }
            // Check same matched object between each input objects
            while( !matched ) // Matching
            {
                matched=true;
                for( int i=0; i<inSize; i++ ){
                    for( int j=0; j<i; j++ ){
                        if( iMatchedObj[i] == iMatchedObj[j] )
                        { 
                            double dRi = dR[i][iMatchedIdx[i]];
                            double dRj = dR[j][iMatchedIdx[j]];
                            if( dRi < dRj )
                            {
                                int originIdx = iMatchedIdx[j];
                                iMatchedIdx[j] = originIdx+1;             // move to next
                                iMatchedObj[j] = iObj[j][iMatchedIdx[j]]; // move to next
                            }else{
                                int originIdx = iMatchedIdx[i];
                                iMatchedIdx[i] = originIdx+1;             // move to next
                                iMatchedObj[i] = iObj[i][iMatchedIdx[i]]; // move to next
                            }   
                            matched=false;
                        }
                    }
                }
            }
            // output the matched objects with repceted to index of input collections 
            for( int i=0; i<inSize; i++ )
            {
                outcol.push_back( mcol[iMatchedObj[i]] );
            }    
            return true;
        }

    //* check if the mo1 or mo2 = moPdgID
    bool checkMo( GenParticle GenInfo, int moPdgID )
    {
        bool isMyMo=false;
        if( GenInfo.Mo1PdgID == moPdgID || GenInfo.Mo2PdgID == moPdgID ) isMyMo = true;
        return isMyMo;
    }
    bool checkMo( GenInfoBranches GenInfo, int idx, int moPdgID )
    {
        bool isMyMo=false;
        if( GenInfo.Mo1PdgID[idx] == moPdgID || GenInfo.Mo2PdgID[idx] == moPdgID ) isMyMo = true;
        return isMyMo;
    }

    //* check if the da1 or da2 = daPdgID
    bool checkDa( GenParticle GenInfo, int daPdgID )
    {
        bool isMyDa=false;
        if( GenInfo.Da1PdgID == daPdgID || GenInfo.Da2PdgID == daPdgID ) isMyDa = true;
        return isMyDa;
    }
    bool checkDa( GenInfoBranches GenInfo, int idx, int daPdgID )
    {
        bool isMyDa=false;
        if( GenInfo.Da1PdgID[idx] == daPdgID || GenInfo.Da2PdgID[idx] == daPdgID ) isMyDa = true;
        return isMyDa;
    }

    //* check if the lepton is isolated from jets
    bool isIsoLeptonFromJets( Lepton lepton, vector<Jet> jetCol, double dR )
    {
        bool isIsoLepFromJets = true;
        int const jetcolSize = jetCol.size();
        for( int i=0; i<jetcolSize; i++ )
        {
            double dr = lepton.P4.DeltaR( jetCol[i].P4 );
            if( dr < dR )
            { 
                isIsoLepFromJets = false;
                break;
            }
        }
        return isIsoLepFromJets;
    }

    //* Remove electron overpapping muon
    void rmElelectronOverlapeMuon( vector<Lepton> &elColl, vector<Lepton> muColl, double dR )
    {
        if( elColl.size() == 0 || muColl.size() == 0 ) return;
        vector<Lepton> newElColl;
        for( vector<Lepton>::const_iterator el = elColl.begin(); el != elColl.end(); el++ ){
            for( vector<Lepton>::const_iterator mu = muColl.begin(); mu != muColl.end(); mu++ ){
                if( el->P4.DeltaR( mu->P4 ) > dR ) newElColl.push_back(*el); 
            }
        }
        elColl.swap(newElColl);
    }


    //* Get chi2 from 2 non-b jets and 1 b-jet
    float getChi2( Jet jet1, Jet jet2, Jet bjet, float M_top, float Wth_top, float M_W, float Wth_W )
    {
        TLorentzVector qq_v  = jet1.P4 + jet2.P4;
        TLorentzVector bqq_v = bjet.P4 + qq_v;
        float iTop = ( bqq_v.M() - M_top )/Wth_top;
        float iW   = (  qq_v.M() - M_W   )/Wth_W;
        return iTop*iTop + iW*iW;
    }

    //* Observable 2
    double Obs2( TVector3 isoLep, TVector3 hardJet, TVector3 b, TVector3 bbar )
    {
        TVector3 O2_1v = b + bbar;
        TVector3 O2_2v = isoLep.Cross( hardJet );
        double O2 = O2_1v.Dot( O2_2v );
        return O2;
    }

    //* Observable 3 (boost to bbbar C.M.)
    double Obs3( TLorentzVector isoLep, TLorentzVector hardJet, TLorentzVector b, TLorentzVector bbar, int charge )
    {
        if( abs(charge) != 1 ){ std::cout<<">> [ERROR] SemiLeptanicAnalysis::Obs3 input strange charge: "<<charge<<", it shall be 1 or -1"<<std::endl; return 0.; }

        // In CM(b,b~), b and b~ will back-to-back, i.e vector(b)=-vector(b~)
        TVector3 bbCM = -( b + bbar ).BoostVector();
        b.Boost(bbCM); isoLep.Boost(bbCM); hardJet.Boost(bbCM);

        TVector3 bV       = b.Vect(); 
        TVector3 isoLepV  = isoLep.Vect();
        TVector3 hardJetV = hardJet.Vect(); 

        double O3 = double(charge)*(bV.Dot(isoLepV.Cross(hardJetV)));
        return O3;
    }

    //* Observable 4
    double Obs4( TVector3 isoLep, TVector3 hardJet, TVector3 b, TVector3 bbar, int charge )
    {
        if( abs(charge) != 1 ){ std::cout<<">> [ERROR] SemiLeptanicAnalysis::Obs4 input strange charge: "<<charge<<", it shall be 1 or -1"<<std::endl; return 0.; }
        TVector3 O4_1v = b - bbar;
        TVector3 O4_2v = isoLep.Cross( hardJet );
        double O4 = double(charge)*(O4_1v.Dot( O4_2v ));
        return O4;
    }

    //* Observable 7
    double Obs7( TVector3 beam, TVector3 b, TVector3 bbar )
    {
        double O7_1z = beam.Dot( b - bbar );
        double O7_2z = beam.Dot( b.Cross( bbar ));
        double O7 = O7_1z * O7_2z;
        return O7;
    }
}
#endif
