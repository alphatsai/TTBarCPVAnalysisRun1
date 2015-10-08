#ifndef BTAGSFUTIL_H 
#define BTAGSFUTIL_H 

#include <string>
#include "Jet.h"
#include "SFlightFuncs_EPS2013.h"

class BTagSFUtil
{
    public:
        BTagSFUtil(){};
        ~BTagSFUtil(){};
        float getSF( std::string tagType, Jet jet, int shift=0 )
        {
            return getSF( tagType, jet.Pt, jet.Eta, jet.GenFlavor, shift );
        }
        float getSF( std::string tagType, float pt, float eta, int flavor, int shift=0 ) 
        {   
            float weight=1; 
            if( abs(flavor) == 5 )
                weight=getSFb( tagType, pt, shift );
            else if( abs(flavor) == 4 ) // c use the same SF from b, but twice of unc.
                weight=getSFb( tagType, pt, shift, true );
            else 
                weight=getSFl( tagType, pt, eta, shift );
            return weight;
        }
        float getSFb( std::string tagType, float pt, int shift=0, bool isC=false );
        float getSFl( std::string tagType, float pt, float eta, int shift=0 );

    private:

};

float BTagSFUtil::getSFb( std::string tagType, float pt, int shift, bool isC )
{
    if( abs(shift) != 1 && shift != 0 )
    {
        std::cout<<">> [ERROR] BTagSFUtil::GetSFb shift should not be "<<shift<<std::endl; 
        std::cout<<">>         Accept 0(nominal), 1(sigma), -1(-sigma)"<<std::endl; 
        return 1.; 
    }

    float weight=1;
    float uncScale=1;
    float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
    float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    float error=0;
    float SFb;

    int bin=-1;
    if( isC ) uncScale=2.;
    if( pt < 20  ){ pt = 20.;   uncScale=2; }
    if( pt > 800 ){ pt = 799.9; uncScale=2; }
    for( int i=0; i<int(sizeof(ptmin)/sizeof(ptmin[0])); ++i)
    {
        if( pt >= ptmin[i] && pt < ptmax[i] )
        { 
            bin=i;
            break;
        }
    }
    if( bin == -1 ) std::cout<<">> [ERROR] BTagSFUtil::GetSFb Wrong bin value!"<<std::endl;

    // Tagger: CSVL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    if( tagType.compare("CSVL") == 0 )  
    {
        SFb = 0.997942*((1.+(0.00923753*pt))/(1.+(0.0096119*pt)));
        float SFb_error[] = {
            0.033299,
            0.0146768,
            0.013803,
            0.0170145,
            0.0166976,
            0.0137879,
            0.0149072,
            0.0153068,
            0.0133077,
            0.0123737,
            0.0157152,
            0.0175161,
            0.0209241,
            0.0278605,
            0.0346928,
            0.0350099 };
        error = SFb_error[bin];
    }
    // Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    else if( tagType.compare("CSVM") == 0 ) 
    { 
        SFb = (0.938887+(0.00017124*pt))+(-2.76366e-07*(pt*pt));
        float SFb_error[] = {
            0.0415707,
            0.0204209,
            0.0223227,
            0.0206655,
            0.0199325,
            0.0174121,
            0.0202332,
            0.0182446,
            0.0159777,
            0.0218531,
            0.0204688,
            0.0265191,
            0.0313175,
            0.0415417,
            0.0740446,
            0.0596716 };
        error = SFb_error[bin];
    }
    // Tagger: CSVT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
    else if( tagType.compare("CSVT") == 0 ) 
    { 
        SFb = (0.927563+(1.55479e-05*pt))+(-1.90666e-07*(pt*pt));
        float SFb_error[] = {
            0.0515703,
            0.0264008,
            0.0272757,
            0.0275565,
            0.0248745,
            0.0218456,
            0.0253845,
            0.0239588,
            0.0271791,
            0.0273912,
            0.0379822,
            0.0411624,
            0.0786307,
            0.0866832,
            0.0942053,
            0.102403 };
        error = SFb_error[bin];
    }
    else
    { 
        std::cout<<">> [ERROR] BTagSFUtil::GetSFb tagType not found "<<tagType<<std::endl; 
        std::cout<<">>         Accept CSVL CSVM or CSMT"<<std::endl; 
        return weight; 
    }

    weight = SFb*( 1 + shift*uncScale*error );  
    return weight;

}

float BTagSFUtil::getSFl( std::string tagType, float pt, float eta, int shift )
{
    if( abs(shift) != 1 && shift != 0 )
    {
        std::cout<<">> [ERROR] BTagSFUtil::GetSFl shift should not be "<<shift<<std::endl; 
        std::cout<<">>         Accept 0(nominal), 1(sigma), -1(-sigma)"<<std::endl; 
        return 1.; 
    }
    if( abs(eta) > 2.4 )
    {
        std::cout<<">> [ERROR] BTagSFUtil::GetSFl eta should not be "<<eta<<std::endl; 
        std::cout<<">>         Accept abs(eta) < 2.4"<<std::endl; 
        return 1.; 
    }
    eta = abs(eta);

    float weight=1;
    if( tagType.compare("CSVL") == 0 )
    {  
        for ( int i=0; i<4; ++i) 
        {
            if( eta >= SFlight_CSVL_etamin[i] && eta < SFlight_CSVL_etamax[i] ) 
            {
                     if( shift > 0 ) weight = GetSFlmax("CSV","L",SFlight_CSVL_etamin[i], SFlight_CSVL_etamax[i], "ABCD")->Eval(pt); 
                else if( shift < 0 ) weight = GetSFlmin("CSV","L",SFlight_CSVL_etamin[i], SFlight_CSVL_etamax[i], "ABCD")->Eval(pt); 
                else weight = (GetSFlmean("CSV","L",SFlight_CSVL_etamin[i], SFlight_CSVL_etamax[i], "ABCD"))->Eval(pt); 
                break ;
            }
        }
    }
    else if( tagType.compare("CSVM") == 0 ) 
    { 
        for (int i= 0; i< 4; ++i) 
        {
            if( eta >= SFlight_CSVM_etamin[i] && eta < SFlight_CSVM_etamax[i] ) 
            {
                     if( shift > 0) weight = GetSFlmax("CSV","M",SFlight_CSVM_etamin[i], SFlight_CSVM_etamax[i], "ABCD")->Eval(pt); 
                else if( shift < 0) weight = GetSFlmin("CSV","M",SFlight_CSVM_etamin[i], SFlight_CSVM_etamax[i], "ABCD")->Eval(pt); 
                else weight = (GetSFlmean("CSV","M",SFlight_CSVM_etamin[i], SFlight_CSVM_etamax[i], "ABCD"))->Eval(pt); 
                break ;
            }
        }
    }
    else if( tagType.compare("CSVT") == 0 )
    { 
             if( shift > 0) weight = GetSFlmax("CSV","T",0.0, 2.4, "ABCD")->Eval(pt); 
        else if( shift < 0) weight = GetSFlmin("CSV","T",0.0, 2.4, "ABCD")->Eval(pt);  
        else weight = GetSFlmean("CSV","T",0.0, 2.4, "ABCD")->Eval(pt); 
    }
    else 
    {
        std::cout<<">> [ERROR] BTagSFUtil::GetSFl tagType not found "<<tagType<<std::endl; 
        std::cout<<">>         Accept CSVL CSVM or CSMT"<<std::endl; 
        return weight; 
    }
    return weight ; 
}


#endif
