#include <sstream>
#include <string>
std::string double2str( double i )
{
    std::string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}
////* Caculate ACP
double caculateACPerror( TH1* Oh, bool murmur=false)
{
    int binOg0=2; 
    int binOl0=1; 

    double Og0 = Oh->GetBinContent(binOg0);
    double Ol0 = Oh->GetBinContent(binOl0);

    double n3  = (Og0+Ol0)*(Og0+Ol0)*(Og0+Ol0);
    double ACPe = 2*sqrt(Og0*Ol0/n3);

    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
    printf("Call caculateACPerror: %f\n", ACPe);
    return ACPe;
}
    template<class th1>
double caculateACPerrorWrt( th1* Oh, bool murmur=false)
{
    int binOg0=2; 
    int binOl0=1; 

    double Og0 = Oh->GetBinContent(binOg0);
    double Ol0 = Oh->GetBinContent(binOl0);
    double Og0e = Oh->GetBinError(binOg0);
    double Ol0e = Oh->GetBinError(binOl0);

    double t1 = 2*Ol0/((Ol0+Og0)*(Ol0+Og0));  
    double t2 = 2*Og0/((Ol0+Og0)*(Ol0+Og0)); 
    double ACPe = sqrt( t1*t1*Og0e*Og0e + t2*t2*Ol0e*Ol0e ); 

    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
    printf("Call caculateACPerror: %f\n", ACPe);
    return ACPe;
}
    template<class th1>
double caculateACP( th1* Oh, bool murmur=false)
{
    int binOg0=2; 
    int binOl0=1; 

    double Og0 = Oh->GetBinContent(binOg0);
    double Ol0 = Oh->GetBinContent(binOl0);

    double ACP = (Og0-Ol0)/(Og0+Ol0);
    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", Ol0, Og0);
    printf("Call caculateACP: %f\n", ACP);
    return ACP;
}

template<class th1>
std::string caculateACPDetail( th1* h, bool wrtError=true, std::string Oname="O" )
{
    int binOg0=2; 
    int binOl0=1; 

    double Og0 = h->GetBinContent(binOg0);
    double Ol0 = h->GetBinContent(binOl0);
    double Og0e = h->GetBinError(binOg0);
    double Ol0e = h->GetBinError(binOl0);

    double ACP = caculateACP(h);
    double ACPe; 
    if( wrtError ) ACPe = caculateACPerrorWrt(h);
    else ACPe = caculateACPerror(h);

    printf("%s>0=%5.0f, %s<0=%5.0f, ACP=%6.3f +/-%6.3f \n", Oname.c_str(), Og0, Oname.c_str(), Ol0, ACP, ACPe);
    std::string out = Oname+"+ = "+double2str(Og0)+"("+double2str(Og0e)+"), "+Oname+"- ="+double2str(Ol0)+"("+double2str(Ol0e)+"), ACP="+double2str(ACP)+" +/- "+double2str(ACPe);
    return out; 
}
void caculateACPDetail( TFile* f, std::string histName, std::string Oname="O")
{
    TH1D* Oh = (TH1D*)f->Get(histName.c_str());
    printf("%-15s: ", histName.c_str());
    std::cout<<caculateACPDetail(Oh, Oname)<<std::endl;;
}

////* Caculate Cutflow
void getCutFlowNum( TH1D* h_cutflow, std::string name="" )
{
    int minBin = 1;
    int maxBin = h_cutflow->GetNbinsX();
    printf("%s\n", name.c_str());
    printf("%15s %10s %10s\n", "Selection", "Events", "Eff(%)");
    printf("-----------------------------------------\n");
    double lastEvents=0;
    for( int bin=minBin; bin<=maxBin; bin++)
    {
        if( bin == minBin ) lastEvents = h_cutflow->GetBinContent(bin);
        double eff = h_cutflow->GetBinContent(bin)/lastEvents*100;
        char label[100];
        sprintf( label, "%s", h_cutflow->GetXaxis()->GetBinLabel(bin));
        printf("%15s %10.0f %10.1f\n", label, h_cutflow->GetBinContent(bin), eff);
        lastEvents = h_cutflow->GetBinContent(bin);
    }
}
void getCutFlowNum( TFile* f, std::string h_name="", bool printOut=false, std::string output="." )
{
    std::string name = h_name;
    TH1D* h_cutflow = (TH1D*)f->Get(name.c_str());

    FILE * out;
    if( printOut )
    {
        printf("Writing to %s...", (output+"/Table_"+name+".txt").c_str());
        out = fopen ((output+"/Table_"+name+".txt").c_str(),"w");
    }

    int minBin = 1;
    int maxBin = h_cutflow->GetNbinsX();
    printf("\nHistogram: %s\n", name.c_str());
    printf("%15s %10s %10s\n", "Selection", "Events", "Eff(%)");
    printf("-----------------------------------------\n");
    if( printOut )
    {
        out = fopen ((output+"/Table_"+name+".txt").c_str(),"w");
        fprintf(out, "\nHistogram: %s\n", name.c_str());
        fprintf(out, "%15s %10s %10s\n", "Selection", "Events", "Eff(%)");
        fprintf(out, "-----------------------------------------\n");
    }
    double lastEvents=0;
    for( int bin=minBin; bin<=maxBin; bin++)
    {
        if( bin == minBin ) lastEvents = h_cutflow->GetBinContent(bin);
        double eff = h_cutflow->GetBinContent(bin)/lastEvents*100;
        char label[100];
        sprintf( label, "%s", h_cutflow->GetXaxis()->GetBinLabel(bin));
        printf("%15s %10.0f %10.1f\n", label, h_cutflow->GetBinContent(bin), eff);
        if( printOut ) fprintf(out, "%15s %10.0f %10.1f\n", label, h_cutflow->GetBinContent(bin), eff);
        lastEvents = h_cutflow->GetBinContent(bin);
    }
    printf("-----------------------------------------\n");
    if( printOut )
    { 
        fprintf(out, "-----------------------------------------\n");
        fclose(out);
    }

}
