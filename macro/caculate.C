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
    int binObsPos=2; 
    int binObsNeg=1; 

    double ObsPos = Oh->GetBinContent(binObsPos);
    double ObsNeg = Oh->GetBinContent(binObsNeg);

    double n3  = (ObsPos+ObsNeg)*(ObsPos+ObsNeg)*(ObsPos+ObsNeg);
    double ACPe = 2*sqrt(ObsPos*ObsNeg/n3);

    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", ObsNeg, ObsPos);
    printf("Call caculateACPerror: %f\n", ACPe);
    return ACPe;
}
    template<class th1>
double caculateACPerrorWrt( th1* Oh, bool murmur=false)
{
    int binObsPos=2; 
    int binObsNeg=1; 

    double ObsPos = Oh->GetBinContent(binObsPos);
    double ObsNeg = Oh->GetBinContent(binObsNeg);
    double ObsPosUnc = Oh->GetBinError(binObsPos);
    double ObsNegUnc = Oh->GetBinError(binObsNeg);

    double t1 = 2*ObsNeg/((ObsNeg+ObsPos)*(ObsNeg+ObsPos));  
    double t2 = 2*ObsPos/((ObsNeg+ObsPos)*(ObsNeg+ObsPos)); 
    double ACPe = sqrt( t1*t1*ObsPosUnc*ObsPosUnc + t2*t2*ObsNegUnc*ObsNegUnc ); 

    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", ObsNeg, ObsPos);
    printf("Call caculateACPerror: %f\n", ACPe);
    return ACPe;
}
    template<class th1>
double caculateACP( th1* Oh, bool murmur=false)
{
    int binObsPos=2; 
    int binObsNeg=1; 

    double ObsPos = Oh->GetBinContent(binObsPos);
    double ObsNeg = Oh->GetBinContent(binObsNeg);

    double ACP = (ObsPos-ObsNeg)/(ObsPos+ObsNeg);
    if ( murmur ) printf("O<0 : %.1f, O>0 : %.1f\n", ObsNeg, ObsPos);
    printf("Call caculateACP: %f\n", ACP);
    return ACP;
}

template<class th1>
std::string caculateACPDetail( th1* h, bool wrtError=true, std::string Oname="O" )
{
    int binObsPos=2; 
    int binObsNeg=1; 

    double ObsPos = h->GetBinContent(binObsPos);
    double ObsNeg = h->GetBinContent(binObsNeg);
    double ObsPosUnc = h->GetBinError(binObsPos);
    double ObsNegUnc = h->GetBinError(binObsNeg);

    double ACP = caculateACP(h);
    double ACPe; 
    if( wrtError ) ACPe = caculateACPerrorWrt(h);
    else ACPe = caculateACPerror(h);

    printf("%s>0=%5.0f, %s<0=%5.0f, ACP=%6.3f +/-%6.3f \n", Oname.c_str(), ObsPos, Oname.c_str(), ObsNeg, ACP, ACPe);
    std::string out = Oname+"+ = "+double2str(ObsPos)+"("+double2str(ObsPosUnc)+"), "+Oname+"- ="+double2str(ObsNeg)+"("+double2str(ObsNegUnc)+"), ACP="+double2str(ACP)+" +/- "+double2str(ACPe);
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
void getCutFlowNum( TFile* f, std::string h_name="", bool printOut=false, std::string output=".", std::string h_name_add="" )
{
    std::string name = h_name;
    TH1D* h_cutflow = (TH1D*)((TH1D*)f->Get(name.c_str()))->Clone();
    if( h_name_add.compare("") != 0 )
    {
        TH1D* h_cutflow_add = (TH1D*)f->Get(h_name_add.c_str());
        h_cutflow->Add(h_cutflow_add);
        name.append(h_name_add);
    }

    FILE * out;
    if( printOut )
    {
        printf("Writing to %s...", (output+"/Table_"+name+".txt").c_str());
        out = fopen ((output+"/Table_"+name+".txt").c_str(),"w");
    }

    int minBin = 1;
    int maxBin = h_cutflow->GetNbinsX();
    printf("\nHistogram: %s\n", name.c_str());
    printf("%15s %10s %10s %10s\n", "Selection", "Events", "Unc.", "Eff(%)");
    printf("------------------------------------------------------\n");
    if( printOut )
    {
        out = fopen ((output+"/Table_"+name+".txt").c_str(),"w");
        fprintf(out, "\nHistogram: %s\n", name.c_str());
        fprintf(out, "%15s %10s %10s %10s\n", "Selection", "Events", "Unc.", "Eff(%)");
        fprintf(out, "------------------------------------------------------\n");
    }
    double lastEvents=0;
    for( int bin=minBin; bin<=maxBin; bin++)
    {
        if( bin == minBin ) lastEvents = h_cutflow->GetBinContent(bin);
        double eff = h_cutflow->GetBinContent(bin)/lastEvents*100;
        char label[100];
        sprintf( label, "%s", h_cutflow->GetXaxis()->GetBinLabel(bin));
        printf("%15s %10.0f %10.0f %10.1f\n", label, h_cutflow->GetBinContent(bin), h_cutflow->GetBinError(bin), eff);
        if( printOut ) fprintf(out, "%15s %10.0f %10.0f %10.1f\n", label, h_cutflow->GetBinContent(bin), h_cutflow->GetBinError(bin), eff);
        lastEvents = h_cutflow->GetBinContent(bin);
    }
    printf("------------------------------------------------------\n");
    if( printOut )
    { 
        fprintf(out, "------------------------------------------------------\n");
        fclose(out);
    }

}
