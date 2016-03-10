#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TColor.h"
#include "../caculate.C"
#include "../CMS_lumi.C"
#include "../help.C"
const int cmslumi=0;
const float topMassUncRescale=6;
void drawDifference( TH1D* h_data, TH1D* h_fit, float yRange=0.5, std::string xtitle="M_{eb} [GeV]", std::string ytitle="#frac{Data}{Fitted} / 10 GeV")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TH1D* h_diff = (TH1D*)h_data->Clone();
    h_diff->Divide(h_fit);

    TCanvas *c1 = new TCanvas("c1", "Histos1",64,226,790,429);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c1->Range(-119.469,0.2445972,538.9381,1.570727);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1814516);
    c1->SetRightMargin(0.05913978);
    c1->SetTopMargin(0.05333334);
    c1->SetBottomMargin(0.1925926);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    h_diff->GetYaxis()->SetRangeUser( 1-yRange, 1+yRange );
    h_diff->SetLineWidth(2);
    h_diff->GetXaxis()->SetTitle(xtitle.c_str());
    h_diff->GetXaxis()->SetNdivisions(509);
    h_diff->GetXaxis()->SetLabelFont(42);
    h_diff->GetXaxis()->SetLabelSize(0.06);
    h_diff->GetXaxis()->SetTitleSize(0.06);
    h_diff->GetXaxis()->SetTitleFont(42);
    h_diff->GetYaxis()->SetTitle(ytitle.c_str());
    h_diff->GetYaxis()->SetNdivisions(505);
    //h_diff->GetYaxis()->SetNdivisions(504);
    h_diff->GetYaxis()->SetLabelFont(42);
    h_diff->GetYaxis()->SetLabelSize(0.06);
    h_diff->GetYaxis()->SetTitleSize(0.07);
    h_diff->GetYaxis()->SetTitleFont(42);
    h_diff->Draw();

    TLine* l = new TLine( h_diff->GetXaxis()->GetBinLowEdge(1), 1, h_diff->GetXaxis()->GetBinUpEdge(h_diff->GetXaxis()->GetLast()), 1);
    l->SetLineColor(2);
    l->SetLineWidth(2);
    l->Draw("SAME");
    h_diff->Draw("SAME");
}

float* getEvtPDF( TFile* fin, std::string hName, int nPDF )
{
    float *values = new float[2];
    TH1D* h_mean = (TH1D*)((TH1D*)fin->Get(hName.c_str()))->Clone(); 
    TH1D* h_pdf  = (TH1D*)((TH1D*)fin->Get((hName+"_PDF"+int2str(nPDF)).c_str()))->Clone();
    values[0]=h_mean->Integral();
    values[1]=h_pdf->Integral();
    return values;
}
float* drawSyst( TFile* fin, std::string hName, std::string systName, std::string output=".", std::string xTitle="", std::string yTitle="Events", bool logy=false )
{
    float *values = new float[3];
    TH1D* h_mean = (TH1D*)((TH1D*)fin->Get(hName.c_str()))->Clone(); 
    TH1D* h_up   = (TH1D*)((TH1D*)fin->Get((hName+"_"+systName+"up").c_str()))->Clone();
    TH1D* h_down = (TH1D*)((TH1D*)fin->Get((hName+"_"+systName+"down").c_str()))->Clone();

    values[0]=h_mean->Integral();
    values[1]=h_up->Integral();
    values[2]=h_down->Integral();
    cout<<values[0]<<" +"<<values[1]<<" -"<<values[2]<<endl;

    if( logy ){
        if( h_mean->GetMaximum() < h_up->GetMaximum())   h_mean->SetMaximum(h_up->GetMaximum()*10);
        if( h_mean->GetMaximum() < h_down->GetMaximum()) h_mean->SetMaximum(h_down->GetMaximum()*10);
        if( h_mean->GetMinimum() > 0 )
            h_mean->SetMinimum(h_mean->GetMinimum()/10);
        else
            h_mean->SetMinimum(10);
    }else{
        if( h_mean->GetMaximum() < h_up->GetMaximum() )  h_mean->SetMaximum(h_up->GetMaximum()+h_up->GetMaximum()/10);
        if( h_mean->GetMaximum() < h_down->GetMaximum()) h_mean->SetMaximum(h_down->GetMaximum()+h_down->GetMaximum()/10);
    }
    h_mean->GetXaxis()->SetTitle(xTitle.c_str());
    h_mean->GetXaxis()->SetRange(1,500);
    h_mean->GetXaxis()->SetNdivisions(508);
    h_mean->GetXaxis()->SetLabelFont(42);
    h_mean->GetXaxis()->SetLabelSize(0.06);
    h_mean->GetXaxis()->SetTitleSize(0.07);
    h_mean->GetXaxis()->SetTitleOffset(0.96);
    h_mean->GetXaxis()->SetTitleFont(42);
    h_mean->GetYaxis()->SetTitle(yTitle.c_str());
    h_mean->GetYaxis()->SetNdivisions(503);
    h_mean->GetYaxis()->SetLabelFont(42);
    h_mean->GetYaxis()->SetLabelSize(0.05);
    h_mean->GetYaxis()->SetTitleSize(0.07);
    h_mean->GetYaxis()->SetTitleOffset(1.01);
    h_mean->SetLineColor(1);
    h_mean->SetLineWidth(3);
    h_up->SetLineColor(2);
    h_up->SetLineStyle(2);
    h_up->SetLineWidth(3);
    h_down->SetLineColor(4);
    h_down->SetLineStyle(2);
    h_down->SetLineWidth(3);

    TCanvas* c1;
    if( logy ) 
        c1 = new TCanvas( ("C_Log_"+hName).c_str(), "", 15,94,1016,824);
    else
        c1 = new TCanvas( ("C_Linear_"+hName).c_str(), "", 15,94,1016,824);
    c1->Range(-123.385,-1839.168,569.1214,8340.115);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1781716);
    c1->SetTopMargin(0.06022585);
    c1->SetBottomMargin(0.1806775);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if( logy ) c1->SetLogy(1);
    else c1->SetLogy(0);

    char text[50];
    TLegend *leg;
    leg = new TLegend(0.5382463,0.6900878,0.8544776,0.9209536,NULL,"brNDC");
    leg->SetHeader(systName.c_str());
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03764115);
    leg->SetLineColor(1);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    sprintf( text, "Nominal, %5.0f events", values[0] ); cout<<text<<endl;
    leg->AddEntry(h_mean, text, "l");
    sprintf( text, "+1 #sigma, %5.0f events", values[1] ); cout<<text<<endl;
    leg->AddEntry(h_up,   text, "l");
    sprintf( text, "-1 #sigma, %5.0f events", values[2] ); cout<<text<<endl;
    leg->AddEntry(h_down, text, "l");

    TPaveText* t_title;
    t_title = new TPaveText(0.1427239,0.944793,0.7723881,0.9899624,"brNDC");
    t_title->AddText("CMS #sqrt{s} = 8TeV, L = 19.7/fb");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);

    h_mean->Draw("HIST");
    h_up->Draw("HISTSAME");
    h_down->Draw("HISTSAME");
    leg->Draw();
    t_title->Draw();

    if( logy )
        c1->SaveAs((output+"/Syst_Log_"+hName+"_"+systName+".pdf").c_str());
    else
        c1->SaveAs((output+"/Syst_Linear_"+hName+"_"+systName+".pdf").c_str());
    return values;
}

void drawFittedStack( TFile* f, std::string hName, std::string* systName, int systNi, int nPDFup=0, int nPDFdown=0, int chN=0, float minMlb=0, float maxMlb=0, std::string output=".", std::string xTitle="", std::string yTitle="Events", bool logy=false, bool doTopMassUncRescale=true, bool cuthist=false )
{
    int nCh=2;
    int lineWidth=3;

    std::string ch;
    std::string channel;
    if( chN == 1 ){
        ch="_El";
        channel="Electron channel";
    }else if( chN == 2){
        ch="_Mu";
        channel="Muon channel";
    }else{
        ch="";
        channel="Combined channel";
    }

    const int systN = systNi+1;
    TH1D* hs;
    TH1D* h_data;
    TH1D* h_bkg;
    TH1D* h_sig;
    TH1D* h_all;
    TH1D* h_allUnc[systN][2];
    double uncSig[systN][2]; 
    double uncBkg[systN][2]; 
    double uncAll[systN][2];
    double meanData; 
    double meanSig; 
    double meanBkg; 
    double meanAll; 

    h_data = (TH1D*)((TH1D*)f->Get((("DATA"+ch).c_str())))->Clone("DATA");
    h_sig  = (TH1D*)((TH1D*)f->Get((("SigFitted"+ch).c_str())))->Clone("SIG"); 
    h_bkg  = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch).c_str())))->Clone("BKG"); 
    h_all  = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch).c_str())))->Clone("ALL"); 
    h_all->Add(h_sig);

    int bins = h_data->GetNbinsX();
    if( minMlb == maxMlb )
    {
        minMlb = 1;
        maxMlb = bins;
    }
    else
    {
        for( int b=1; b<=h_data->GetNbinsX(); b++ )
        {
            if( h_data->GetXaxis()->GetBinLowEdge(b) == minMlb ) minMlb=b;
            if( h_data->GetXaxis()->GetBinLowEdge(b) == maxMlb ) maxMlb=b-1;
        }
    }
    meanData = h_data->Integral(minMlb,maxMlb);
    meanSig  = h_sig ->Integral(minMlb,maxMlb);
    meanBkg  = h_bkg ->Integral(minMlb,maxMlb);
    meanAll  = h_all ->Integral(minMlb,maxMlb);

    int iTopMass=-1;
    std::string tune[2] = {"up", "down"};
    for( int i=0; i<systN; i++ )
    {   
        //cout<<systName[i]<<endl;
        TH1D* h_sigUnc[2];
        if( i != systN-1 )
        {
            if( systName[i].find("TopMass")!=std::string::npos && doTopMassUncRescale ) iTopMass=i;
            for( int t=0; t<2; t++)
            {
                std::string name = systName[i]+tune[t];
                cout<<name<<endl;
                h_sigUnc[t]    = (TH1D*)f->Get((("SigFitted"+ch+"_"+name).c_str()));
                h_allUnc[i][t] = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch+"_"+name).c_str())))->Clone(("UNC_"+name).c_str());
                uncSig[i][t]=h_sigUnc[t]->Integral(minMlb,maxMlb);
                uncBkg[i][t]=h_allUnc[i][t]->Integral(minMlb,maxMlb); h_allUnc[i][t]->Add(h_sigUnc[t]);
                uncAll[i][t]=h_allUnc[i][t]->Integral(minMlb,maxMlb);
                if( iTopMass == i )
                {
                    uncSig[i][t]=(uncSig[i][t]+(topMassUncRescale-1)*meanSig)/topMassUncRescale;
                    uncBkg[i][t]=(uncBkg[i][t]+(topMassUncRescale-1)*meanBkg)/topMassUncRescale;
                    uncAll[i][t]=(uncAll[i][t]+(topMassUncRescale-1)*meanAll)/topMassUncRescale;
                }
            }
        }
        else
        { 
            // For PDF
            cout<<"PDF"<<endl;
            h_sigUnc[0]    = (TH1D*)f->Get((("SigFitted"+ch+"_PDF"+int2str(nPDFup)).c_str()));
            h_sigUnc[1]    = (TH1D*)f->Get((("SigFitted"+ch+"_PDF"+int2str(nPDFdown)).c_str()));
            h_allUnc[i][0] = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch+"_PDF"+int2str(nPDFup)).c_str())))->Clone("UNC_PDFup");
            h_allUnc[i][1] = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch+"_PDF"+int2str(nPDFdown)).c_str())))->Clone("UNC_PDFdown");
            uncSig[i][0]=h_sigUnc[0]->Integral(minMlb,maxMlb);
            uncSig[i][1]=h_sigUnc[1]->Integral(minMlb,maxMlb);
            uncBkg[i][0]=h_allUnc[i][0]->Integral(minMlb,maxMlb); h_allUnc[i][0]->Add(h_sigUnc[0]);
            uncBkg[i][1]=h_allUnc[i][1]->Integral(minMlb,maxMlb); h_allUnc[i][1]->Add(h_sigUnc[1]);
            uncAll[i][0]=h_allUnc[i][0]->Integral(minMlb,maxMlb);
            uncAll[i][1]=h_allUnc[i][1]->Integral(minMlb,maxMlb);
        }
        delete *h_sigUnc;
    }

    //// Draw plots
    TGraphAsymmErrors* h_allAsymErr = new TGraphAsymmErrors(h_all);
    for( int b=1; b<=bins; b++ )
    {
        float nomV = h_all->GetBinContent(b);
        float sigUp   = h_all->GetBinContent(b); // Has been cross check with h_all->IntegralAndError(0,50,error); 
        float sigDown = h_all->GetBinContent(b); // cout<<error; The error is after weighting
        //float sigUp   = 0;
        //float sigDown = 0;
        for( int i=0; i<systN; i++ ){
            for( int t=0; t<2; t++ )
            {
                float err = h_allUnc[i][t]->GetBinContent(b) - nomV;
                if( iTopMass == i ) err = err/topMassUncRescale;
                if( err > 0. ) sigUp +=err*err;
                else sigDown += err*err;
            }
        }
        //cout<<b<<": "<<sqrt(sigUp)<<" "<<sqrt(sigDown)<<endl;
        h_allAsymErr->SetPointEYhigh( b-1, sqrt(sigUp));
        h_allAsymErr->SetPointEYlow(  b-1, sqrt(sigDown));
    }

    h_data->SetLineWidth(lineWidth);
    h_data->SetLineColor(1);
    h_data->SetMarkerColor(1);
    h_data->SetMarkerStyle(8);

    h_sig->SetLineWidth(lineWidth-1);
    h_sig->SetLineColor(1);
    h_sig->SetFillColor(50);

    h_bkg->SetLineWidth(lineWidth-1);
    h_bkg->SetLineColor(1);
    h_bkg->SetFillColor(8);

    h_allAsymErr->SetLineColor(13);
    h_allAsymErr->SetFillColor(13);
    h_allAsymErr->SetFillStyle(3001);

    float xMin = h_all->GetXaxis()->GetBinLowEdge(1);
    float xMax = h_all->GetXaxis()->GetBinUpEdge(bins);
    if( logy )
        hs = new TH1D(("TH1DinStackLog"+hName+channel).c_str(), "", bins, xMin, xMax);
    else
        hs = new TH1D(("TH1DinStackLinear"+hName+channel).c_str(), "", bins, xMin, xMax);

    if( cuthist ) hs->GetXaxis()->SetRange(minMlb, maxMlb);
    hs->GetXaxis()->SetTitle(xTitle.c_str());
    hs->GetYaxis()->SetTitle(yTitle.c_str());
    hs->GetXaxis()->SetNdivisions(509);
    hs->GetXaxis()->SetLabelFont(42);
    hs->GetXaxis()->SetLabelSize(0.06);
    hs->GetXaxis()->SetTitleSize(0.07);
    hs->GetXaxis()->SetTitleOffset(1.03);
    hs->GetXaxis()->SetTitleFont(42);
    hs->GetYaxis()->SetLabelFont(42);
    hs->GetYaxis()->SetLabelSize(0.06);
    hs->GetYaxis()->SetTitleSize(0.079);
    hs->GetYaxis()->SetTitleOffset(1.11);
    hs->GetYaxis()->SetTitleFont(42);

    THStack* h_stack = new THStack("THStcak", "");
    h_stack->SetHistogram(hs);

    h_stack->Add(h_bkg);
    h_stack->Add(h_sig);

    TCanvas* c1;
    TPad *p1, *p2; 
    if( logy ){
        c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",59,72,W,H);
    }else{
        c1 = new TCanvas( ("C_Linear_"+hName).c_str(), "",59,72,W,H);
    }
    c1->Range(-123.0964,-841.8412,557.1066,3883.14);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1997487);
    c1->SetRightMargin(0.06909548);
    c1->SetTopMargin(0.08900524);
    c1->SetBottomMargin(0.1797557);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if( logy ) c1->SetLogy(1);
    else c1->SetLogy(0);

    TLegend *leg;
    leg = new TLegend(0.6256281,0.617801,0.9258794,0.904014,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "Data", "lpe");
    leg->AddEntry(h_sig,  "Estimated sig.", "f");
    leg->AddEntry(h_bkg,  "Estimated bkg.", "f");
    leg->AddEntry(h_allAsymErr, "1#sigma, Stat.+Syst.", "f");
    //leg->AddEntry(h_allAsymErr, "1#sigma, Syst.", "f");

    //TPaveText* t_title;
    //t_title = new TPaveText(0.1455224,0.9134253,0.7751866,0.9974906,"brNDC");
    ////t_title = new TPaveText(0.07088487,0.9153846,0.7007612,1,"brNDC");
    //t_title->AddText("CMS #sqrt{s} = 8TeV, L = 19.7/fb");
    //t_title->SetTextColor(kBlack);
    //t_title->SetFillColor(kWhite);
    //t_title->SetFillStyle(0);
    //t_title->SetBorderSize(0);
    //t_title->SetTextAlign(11);
    //t_title->SetTextSize(0.04805273);

    h_stack->Draw("HIST");
    h_allAsymErr->Draw("E2SAME");
    h_data->Draw("ESAME");
    leg->Draw();
    //t_title->Draw();

    writeExtraText=true;
    CMS_lumi( c1, 2, cmslumi, true );

    if( logy )
        c1->SaveAs((output+"/StackFitted_Log_"+hName+ch+".pdf").c_str());
    else
        c1->SaveAs((output+"/StackFitted_Linear_"+hName+ch+".pdf").c_str());

    //// Summery table 
    float sumw2[2]={0,0}, sumw2sig[2]={0,0}, sumw2bkg[2]={0,0};
    FILE* outTxt;
    outTxt = fopen((output+"/SystUncertainties_"+hName+ch+".txt").c_str(),"w");
    fprintf( outTxt, "%s data: %.0f ( sig: %.f, bkg: %.f )\n", channel.c_str(), meanData, meanSig, meanBkg );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<systNi; i++ ){ fprintf( outTxt, "%9s", systName[i].c_str()); } fprintf( outTxt, "%9s", "PDF" );  
    fprintf( outTxt, "\nTotal:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncAll[i][0]-meanAll)/meanAll*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2[0]+=v*v; } }   
    fprintf( outTxt, "%+9.2f %+9.f", sqrt(sumw2[0]), sqrt(sumw2[0])*meanAll/100. );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncAll[i][1]-meanAll)/meanAll*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2[1]+=v*v; } }  
    fprintf( outTxt, "%+9.2f %+9.f", -sqrt(sumw2[1]), -sqrt(sumw2[1])*meanAll/100. );
    fprintf( outTxt, "\n");
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncSig[i][0]-meanSig)/meanSig*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sig[0]+=v*v; } }   
    fprintf( outTxt, "%+9.2f %+9.f %5.2f(%+5.2f)", sqrt(sumw2sig[0]), sqrt(sumw2sig[0])*meanSig/100., meanSig/meanAll*100., sqrt(sumw2sig[0])*meanSig/meanAll);
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncSig[i][1]-meanSig)/meanSig*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sig[1]+=v*v; } }   
    fprintf( outTxt, "%+9.2f %+9.f %5.2f(%+5.2f)", -sqrt(sumw2sig[1]), -sqrt(sumw2sig[1])*meanSig/100., meanSig/meanAll*100., -sqrt(sumw2sig[1])*meanSig/meanAll );
    fprintf( outTxt, "\n");
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncBkg[i][0]-meanBkg)/meanBkg*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkg[0]+=v*v; }}   
    fprintf( outTxt, "%+9.2f %+9.f", sqrt(sumw2bkg[0]), sqrt(sumw2bkg[0])*meanBkg/100. );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systN; i++ ){ float v=(uncBkg[i][1]-meanBkg)/meanBkg*100.; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkg[1]+=v*v; }}   
    fprintf( outTxt, "%+9.2f %+9.f", -sqrt(sumw2bkg[1]), -sqrt(sumw2bkg[1])*meanBkg/100. );
    fprintf( outTxt, "\n\n");
    fclose( outTxt );
}

void fillFinalAcpText( FILE* outTxt, std::string oName, std::string chName, double Acp, double eAcp, int systN, std::string* systName, double* eAcpSystUp, double* eAcpSystDown )
{
    char text[128];
    sprintf( text, "Nominal Acp(%s) = %+.2lf\n", oName.c_str(), Acp*100. ); 
    fprintf( outTxt, text ); 
     printf( text );
    fprintf( outTxt, "%9s %9s ", " ", "Stat." ); 
     printf( "%9s %9s ", " ", "Stat." ); 
    for( int s=0; s<systN; s++)
    { 
         printf( "%9s ", systName[s].c_str()); 
        fprintf( outTxt, "%9s ", systName[s].c_str() );
    }
    fprintf( outTxt, "\n" );
     printf( "\n" );
    fprintf( outTxt, "%9s", "+1Sigma");
     printf( "%9s ", "+1Sigma");
    sprintf( text, "%+9.2lf ", eAcp*100. ); 
    fprintf( outTxt, text ); 
     printf( text );
     for( int s=0; s<systN; s++ )
     { 
         sprintf( text, "%+9.2lf ", eAcpSystUp[s]*100. ); 
         fprintf( outTxt, text ); 
          printf( text );
     }
    fprintf( outTxt, "\n" );
     printf( "\n" );
    fprintf( outTxt, "%9s", "-1Sigma");
     printf( "%9s ", "-1Sigma");
    sprintf( text, "%+9.2lf ", -1*eAcp*100. ); 
    fprintf( outTxt, text ); 
     printf( text );
     for( int s=0; s<systN; s++ )
     { 
         sprintf( text, "%+9.2lf ", eAcpSystDown[s]*100. ); 
         fprintf( outTxt, text ); 
          printf( text );
     }
    fprintf( outTxt, "\n" );
     printf( "\n" );
}

void getSubtractBkgResults( TFile* f, FILE* outTxt, std::string oName, int systNi, std::string* systName, int nPDFup=0, int nPDFdown=0, int channel=0, float minMlb=0, float maxMlb=0, bool unBlind=false, bool doTopMassUncRescale=true, TFile* fCR=NULL )
{
    std::string chName="";
    std::string chFullName="Combined channel";
    if( channel == 0 ){ chName="_El"; chFullName="Electron channel";} 
    if( channel == 1 ){ chName="_Mu"; chFullName="Muon channel";} 
    std::string tuneName[2]={"up","down"};
    std::string dataName="DATA";
    if( !unBlind ){ dataName="MC"; std::cout<<"[INFO] Blinding! Using MC"<<std::endl;}

    bool usingCR = (fCR==NULL)? false:true;

    const int systN = systNi+1;
    float evtBkgMean;
    float evtBkgFitMean;
    float bkgWrt;
    double  Acp;
    double eAcp;
    double  AcpSyst[2][systN];
    double eAcpSyst[2][systN];

    TH1D* hO_data;
    TH1D* hO_bkg; 
    TH1D* hO_bkgEst; 
    TH1D* hT_bkg;
    TH1D* hO_dataSyst[2][systN];
    TH1D* hO_bkgSyst[2][systN];
    TH1D* hT_bkgSyst[2][systN];   
    float bkgWrtSyst[2][systN];
    
    hO_data = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName).c_str()))->Clone();
    hT_bkg  = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName).c_str()))->Clone();
    if( usingCR )
    {
        // Subtract O from CR
        std::cout<<"[INFO] Bkg O from data CR..."<<std::endl;
        if( channel == 0 ){      hO_bkg  = (TH1D*)((TH1D*)fCR->Get(("DATA_Electron__EvtChi2_"+oName+"Asym_El").c_str()))->Clone(); }
        else if( channel == 1 ){ hO_bkg  = (TH1D*)((TH1D*)fCR->Get(("DATA_Muon__EvtChi2_"+oName+"Asym_Mu").c_str()))->Clone();     }
        else{
            TH1D* hO_tmp = (TH1D*)fCR->Get(("DATA_Electron__EvtChi2_"+oName+"Asym_El").c_str()); 
            hO_bkg  = (TH1D*)((TH1D*)fCR->Get(("DATA_Electron__EvtChi2_"+oName+"Asym_Mu").c_str()))->Clone();
            hO_bkg->Add(hO_tmp);
            delete hO_tmp; 
        }
    }
    else{ hO_bkg  = (TH1D*)((TH1D*)f->Get(("BkgMC_"+oName+chName).c_str()))->Clone(); } // Subtract O from bkg MC

    int bins = hT_bkg->GetNbinsX();
    if( minMlb == maxMlb )
    {
        minMlb = 1;
        maxMlb = bins;
    }
    else
    {
        for( int b=1; b<=bins; b++ )
        {
            if( hT_bkg->GetXaxis()->GetBinLowEdge(b) == minMlb ) minMlb=b;
            if( hT_bkg->GetXaxis()->GetBinLowEdge(b) == maxMlb ) maxMlb=b-1;
        }
    }

    if( !unBlind )
    {
        // Make fake data base on MC 
        TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName).c_str());
        hO_data->Scale(tmpData->Integral()/hO_data->Integral());
        hO_data->SetBinError( 1, sqrt(hO_data->GetBinContent(1)));
        hO_data->SetBinError( 2, sqrt(hO_data->GetBinContent(2)));
        delete tmpData;
    }

    // Weight bkg O to estimated bkg events 
    evtBkgMean  = hO_bkg->Integral();
    evtBkgFitMean = hT_bkg->Integral(minMlb,maxMlb);
    bkgWrt = evtBkgFitMean/evtBkgMean;
    hO_bkgEst = (TH1D*)hO_bkg->Clone();
    hO_bkgEst->Scale(bkgWrt);

    fprintf( outTxt, "* %s (%%)\n", chFullName.c_str() );
     printf( "* %s (%%)\n", chFullName.c_str() );
    fprintf( outTxt, "* Events - Data:%.1f\n", hO_data->Integral() );
     printf( "* Events - Data:%.1f\n", hO_data->Integral() );
    fprintf( outTxt, "* Events - BkgMC:%.1f, EstBkg:%.1f, wrt:%.3f\n", evtBkgMean, evtBkgFitMean, bkgWrt);
     printf( "* Events - BkgMC:%.1f, EstBkg:%.1f, wrt:%.3f\n", evtBkgMean, evtBkgFitMean, bkgWrt);

    fprintf( outTxt, "* ACP(%s) - Data:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_data)*100., caculateACPerrorWrt(hO_data)*100., hO_data->GetBinContent(1), hO_data->GetBinContent(2));
     printf( "* ACP(%s) - Data:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_data)*100., caculateACPerrorWrt(hO_data)*100., hO_data->GetBinContent(1), hO_data->GetBinContent(2));
    fprintf( outTxt, "* ACP(%s) - BkgMC:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_bkg)*100., caculateACPerrorWrt(hO_bkg)*100., hO_bkg->GetBinContent(1), hO_bkg->GetBinContent(2));
     printf( "* ACP(%s) - BkgMC:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_bkg)*100., caculateACPerrorWrt(hO_bkg)*100., hO_bkg->GetBinContent(1), hO_bkg->GetBinContent(2));

    fprintf( outTxt, "* ACP(%s) - EstBkg:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_bkgEst)*100., caculateACPerrorWrt(hO_bkgEst)*100., hO_bkgEst->GetBinContent(1), hO_bkgEst->GetBinContent(2));
     printf( "* ACP(%s) - EstBkg:%.2f(%.2f) n:%.1f, p:%.1f\n", oName.c_str(), caculateACP(hO_bkgEst)*100., caculateACPerrorWrt(hO_bkgEst)*100., hO_bkgEst->GetBinContent(1), hO_bkgEst->GetBinContent(2));

    // Subtract bkg
    hO_data->Add(hO_bkgEst,-1);

    fprintf( outTxt, "* Events - DataSig:%.1f (n:%.1f, p:%.1f)\n", hO_data->Integral(), hO_data->GetBinContent(1), hO_data->GetBinContent(2));
     printf( "* Events - DataSig:%.1f (n:%.1f, p:%.1f)\n", hO_data->Integral(), hO_data->GetBinContent(1), hO_data->GetBinContent(2));

    Acp = caculateACP( hO_data, 1 );
    eAcp = caculateACPerrorWrt( hO_data );
    // * Systmatic uncs.
    int iTopMass=-1;
    string *systName1 = new string[systN];
    for( int s=0; s<systN; s++)
    {
        if( s != systN-1 )
        {
            if( systName[s].find("TopMass")!=std::string::npos && doTopMassUncRescale ) iTopMass=s;
            systName1[s] = systName[s];
            for( int t=0; t<2; t++ ){
                hT_bkgSyst[t][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName+"_"+systName[s]+tuneName[t]).c_str()))->Clone();  
                bkgWrtSyst[t][s] = hT_bkgSyst[t][s]->Integral(minMlb,maxMlb)/evtBkgMean; 
                if( iTopMass == s ){ bkgWrtSyst[t][s] = (bkgWrtSyst[t][s]+(topMassUncRescale-1)*bkgWrt)/topMassUncRescale; }
                hO_bkgSyst[t][s] = (TH1D*)hO_bkg->Clone();
                hO_bkgSyst[t][s]->Scale(bkgWrtSyst[t][s]); 
                hO_dataSyst[t][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName).c_str()))->Clone();
                hO_dataSyst[t][s]->Add(hO_bkgSyst[t][s],-1);
                AcpSyst[t][s] = caculateACP(hO_dataSyst[t][s]);
                eAcpSyst[t][s] = (AcpSyst[t][s]-Acp);
            }
        }
        else
        { 
            // PDF
            systName1[s]="PDF";
            hT_bkgSyst[0][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName+"_PDF"+int2str(nPDFup)).c_str()))->Clone();   
            hT_bkgSyst[1][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName+"_PDF"+int2str(nPDFdown)).c_str()))->Clone();   
            bkgWrtSyst[0][s] = hT_bkgSyst[0][s]->Integral(minMlb,maxMlb)/evtBkgMean; 
            bkgWrtSyst[1][s] = hT_bkgSyst[1][s]->Integral(minMlb,maxMlb)/evtBkgMean; 
            hO_bkgSyst[0][s] = (TH1D*)hO_bkg->Clone();
            hO_bkgSyst[1][s] = (TH1D*)hO_bkg->Clone();
            hO_bkgSyst[0][s]->Scale(bkgWrtSyst[0][s]); 
            hO_bkgSyst[1][s]->Scale(bkgWrtSyst[1][s]); 
            hO_dataSyst[0][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName).c_str()))->Clone();
            hO_dataSyst[1][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName).c_str()))->Clone();
            hO_dataSyst[0][s]->Add(hO_bkgSyst[0][s],-1);
            hO_dataSyst[1][s]->Add(hO_bkgSyst[1][s],-1);
            AcpSyst[0][s] = caculateACP(hO_dataSyst[0][s]);
            AcpSyst[1][s] = caculateACP(hO_dataSyst[1][s]);
            eAcpSyst[0][s] = (AcpSyst[0][s]-Acp);
            eAcpSyst[1][s] = (AcpSyst[1][s]-Acp);
        }
    }
    // * Get data ACP
    fillFinalAcpText( outTxt, oName, chFullName, Acp, eAcp, systN, systName1, eAcpSyst[0], eAcpSyst[1] );
    fprintf( outTxt, "\n" );
    printf( "\n" );
    delete [] systName1;
}

void getSubtractBkgResultsCombined( TFile* f, FILE* outTxt, std::string* systName, int systNi, int nPDFup=0, int nPDFdown=0, float minMlb_=0, float maxMlb_=0, std::string output=".", std::string oName="O2", std::string xTitle="O_{2}", bool unBlind=false, float assumpAcp=0, bool doTopMassUncRescale=true, TFile* fCR=NULL )
{
    std::string tuneName[2]={"up","down"};
    std::string chName[3]={"_El","_Mu",""};
    std::string chFullName[3]={"Electron channel","Muon channel","Combined channel"};
    std::string dataName="DATA";
    if( !unBlind ) dataName="MC";

    bool usingCR = (fCR==NULL)? false:true;
    const int systN = systNi + 1;
    string *systName1 = new string[systN];
    for( int s=0; s<systN; s++){
        if( s != systN-1 ) systName1[s] = systName[s];
        else systName1[s] = "PDF";
    }

    float evtBkgMean[3];
    float evtBkgFitMean[3];
    float bkgWrt[3];
    double  Acp0[3];
    double eAcp0[3];
    double  Acp[3];
    double eAcp[3];
    double  AcpSyst[3][2][systN];
    double eAcpSyst[3][2][systN];
    double eAcpSumw2Syst[3][2];
    double eAcpTotalUnc[3][2];

    TH1D* hO_data[3];
    TH1D* hO_bkg[3]; 
    TH1D* hO_bkgEst[3]; 
    TH1D* hT_bkg[3];
    TH1D* hO_dataSyst[3][2][systN];
    TH1D* hO_bkgSyst[3][2][systN];
    TH1D* hT_bkgSyst[3][2][systN];   
    float bkgWrtSyst[3][2][systN];
    for( int ch=0; ch<3; ch++ )
    {
        hO_data[ch] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();

        if( !unBlind )
        {
            // Make fake data base on MC 
            TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
            float total = tmpData->Integral();
            hO_data[ch]->Scale(total/hO_data[ch]->Integral());
            if( assumpAcp != 0 )
            {
                hO_data[ch]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                hO_data[ch]->SetBinContent( 2, (1+assumpAcp)*total/2 );
            }
            hO_data[ch]->SetBinError( 1, sqrt(hO_data[ch]->GetBinContent(1)));
            hO_data[ch]->SetBinError( 2, sqrt(hO_data[ch]->GetBinContent(2)));
            delete tmpData;
        }
        Acp0[ch] = caculateACP( hO_data[ch] );
        eAcp0[ch] = caculateACPerrorWrt( hO_data[ch] );

        if( ch != 2 ) // Electron or muon channel
        {
            if( usingCR )
            {
                // Subtract O from CR
                std::cout<<"[INFO] Bkg O from data CR..."<<std::endl;
                if( ch == 0 ){      hO_bkg[ch] = (TH1D*)((TH1D*)fCR->Get(("DATA_Electron__EvtChi2_"+oName+"Asym_El").c_str()))->Clone(); }
                else if( ch == 1 ){ hO_bkg[ch] = (TH1D*)((TH1D*)fCR->Get(("DATA_Muon__EvtChi2_"+oName+"Asym_Mu").c_str()))->Clone();     }
            }
            else{ hO_bkg[ch]  = (TH1D*)((TH1D*)f->Get(("BkgMC_"+oName+chName[ch]).c_str()))->Clone(); } // Subtract O from bkg MC

            hT_bkg[ch] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]).c_str()))->Clone();

            int bins = hT_bkg[ch]->GetNbinsX();
            int binW = hT_bkg[ch]->GetBinWidth(1);
            int minMlb, maxMlb;
            if( minMlb_ == maxMlb_ )
            {
                minMlb = 1;
                maxMlb = bins;
            }
            else
            {
                minMlb = minMlb_/binW+1;
                maxMlb = maxMlb_/binW;
            }

            evtBkgMean[ch]  = hO_bkg[ch]->Integral();
            evtBkgFitMean[ch] = hT_bkg[ch]->Integral(minMlb,maxMlb); //cout<<"Here ch "<<ch<<" evtBkgMean[ch] "<<evtBkgMean[ch]<<" evtBkgFitMean[ch] "<<evtBkgFitMean[ch]<<" min "<<minMlb<<" max "<<maxMlb<<endl;
            bkgWrt[ch] = evtBkgFitMean[ch]/evtBkgMean[ch];
            hO_bkgEst[ch] = (TH1D*)hO_bkg[ch]->Clone();  
            hO_bkgEst[ch]->Scale(bkgWrt[ch]);   // weight bkg with fitted bkg
            hO_data[ch]->Add(hO_bkgEst[ch],-1); // subtract bkg

            Acp[ch] = caculateACP( hO_data[ch] );
            eAcp[ch] = caculateACPerrorWrt( hO_data[ch] );
            // * Systmatic uncs.
            int iTopMass=-1;
            for( int s=0; s<systN; s++)
            {
                int pdft[2]={nPDFup, nPDFdown};
                for( int t=0; t<2; t++ ){
                    if( s != systN-1 ){
                        if( systName[s].find("TopMass")!=std::string::npos && doTopMassUncRescale ) iTopMass=s;
                        hT_bkgSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]+"_"+systName[s]+tuneName[t]).c_str()))->Clone();   
                    }else{
                        hT_bkgSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]+"_PDF"+int2str(pdft[t])).c_str()))->Clone();   
                    }
                    bkgWrtSyst[ch][t][s] = hT_bkgSyst[ch][t][s]->Integral(minMlb,maxMlb)/evtBkgMean[ch]; 
                    if( iTopMass == s ){ bkgWrtSyst[ch][t][s] = (bkgWrtSyst[ch][t][s]+(topMassUncRescale-1)*bkgWrt[ch])/topMassUncRescale; }
                    hO_bkgSyst[ch][t][s] = (TH1D*)hO_bkg[ch]->Clone();
                    hO_bkgSyst[ch][t][s]->Scale(bkgWrtSyst[ch][t][s]); 
                    hO_dataSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();
                    if( !unBlind )
                    {
                        // Make fake data base on MC 
                        TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
                        float total = tmpData->Integral();
                        hO_dataSyst[ch][t][s]->Scale(total/hO_dataSyst[ch][t][s]->Integral());
                        if( assumpAcp != 0 )
                        {
                            hO_dataSyst[ch][t][s]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                            hO_dataSyst[ch][t][s]->SetBinContent( 2, (1+assumpAcp)*total/2 );
                        }
                        delete tmpData;
                    }
                    hO_dataSyst[ch][t][s]->Add(hO_bkgSyst[ch][t][s],-1);
                    AcpSyst[ch][t][s] = caculateACP(hO_dataSyst[ch][t][s]);
                    eAcpSyst[ch][t][s] = (AcpSyst[ch][t][s]-Acp[ch]);
                }
            }
        }
        // Combined channel
        else
        {
            hO_bkgEst[ch] = (TH1D*)hO_bkgEst[0]->Clone();
            hO_bkgEst[ch] ->Add(hO_bkgEst[1]);     // Combined electron and muon
            hO_data[ch]->Add(hO_bkgEst[ch],-1); // subtract bkg
             Acp[ch] = caculateACP( hO_data[ch] );
            eAcp[ch] = caculateACPerrorWrt( hO_data[ch] );
            // * Systmatic uncs.
            for( int s=0; s<systN; s++){
                for( int t=0; t<2; t++ ){
                    hO_bkgSyst[ch][t][s] = (TH1D*)hO_bkgSyst[0][t][s]->Clone();
                    hO_bkgSyst[ch][t][s] ->Add(hO_bkgSyst[1][t][s]); 
                    hO_dataSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();
                    if( !unBlind )
                    {
                        // Make fake data base on MC 
                        TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
                        float total = tmpData->Integral();
                        hO_dataSyst[ch][t][s]->Scale(total/hO_dataSyst[ch][t][s]->Integral());
                        if( assumpAcp != 0 )
                        {
                            hO_dataSyst[ch][t][s]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                            hO_dataSyst[ch][t][s]->SetBinContent( 2, (1+assumpAcp)*total/2 );
                        }
                        delete tmpData;
                    }
                    hO_dataSyst[ch][t][s]->Add(hO_bkgSyst[ch][t][s],-1);
                     AcpSyst[ch][t][s] = caculateACP(hO_dataSyst[ch][t][s]);
                    eAcpSyst[ch][t][s] = (AcpSyst[ch][t][s]-Acp[ch]);
                }
            }
        }

        fprintf( outTxt, "* %s (%%)\n", chFullName[ch].c_str() );
         printf( "* %s (%%)\n", chFullName[ch].c_str() );

        // * Sum Syst.
        eAcpSumw2Syst[ch][0]=0.;
        eAcpSumw2Syst[ch][1]=0.;
        for( int s=0; s<systN; s++ )
        {
            if( eAcpSyst[ch][0][s]*eAcpSyst[ch][1][s] < 0 )
            {
                if( eAcpSyst[ch][0][s] > 0 ){
                    eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                    eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                }else{
                    eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                    eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                } 
            }
            else
            {
                if( fabs(eAcpSyst[ch][0][s]) > fabs(eAcpSyst[ch][1][s]) )
                {
                    if( eAcpSyst[ch][0][s] > 0 ) eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                    if( eAcpSyst[ch][0][s] < 0 ) eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                }
                else
                {
                    if( eAcpSyst[ch][1][s] > 0 ) eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                    if( eAcpSyst[ch][1][s] < 0 ) eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                }
            }
        }
        eAcpTotalUnc[ch][0] = sqrt(eAcpSumw2Syst[ch][0]+eAcp[ch]*eAcp[ch]);
        eAcpTotalUnc[ch][1] = sqrt(eAcpSumw2Syst[ch][1]+eAcp[ch]*eAcp[ch]);

        // * Get data ACP
        fprintf( outTxt, "Before subtracting: %+.2f (%.2f)\n", Acp0[ch]*100, eAcp0[ch]*100 );
        fillFinalAcpText( outTxt, oName, chFullName[ch], Acp[ch], eAcp[ch], systN, systName1, eAcpSyst[ch][0], eAcpSyst[ch][1] );
        fprintf( outTxt, "%9s %+9.2f %+9.2f\n", "Sum Syst.", sqrt(eAcpSumw2Syst[ch][0])*100., -1*sqrt(eAcpSumw2Syst[ch][1])*100. );
         printf( "%9s %+9.2f %+9.2f\n", "Sum Syst.", sqrt(eAcpSumw2Syst[ch][0])*100., -1*sqrt(eAcpSumw2Syst[ch][1])*100. );
        fprintf( outTxt, "%9s %+9.2f %+9.2f\n", "Sum Unc.", eAcpTotalUnc[ch][0]*100., -1*eAcpTotalUnc[ch][1]*100. );
         printf( "%9s %+9.2f %+9.2f\n", "Sum Unc.", eAcpTotalUnc[ch][0]*100., -1*eAcpTotalUnc[ch][1]*100. );
        fprintf( outTxt, "\n" );
         printf( "\n" );
    }
    fprintf( outTxt, "\n\n\n" );
     printf( "\n\n\n" );


    // * Plots!
    TH1D* h0 = new TH1D(("ACP_"+oName).c_str(), "", 3, 0, 3);

    h0->GetXaxis()->SetBinLabel(1, (xTitle+"^{e+#mu}").c_str());
    h0->GetXaxis()->SetBinLabel(2, (xTitle+"^{e}").c_str()    );
    h0->GetXaxis()->SetBinLabel(3, (xTitle+"^{#mu}").c_str()  );

    double percent=100; 
    h0->Fill( 0., Acp[2]*percent );
    h0->Fill( 1., Acp[0]*percent );
    h0->Fill( 2., Acp[1]*percent );

    h0->SetBinError( 1, eAcp[2]*percent );
    h0->SetBinError( 2, eAcp[0]*percent );
    h0->SetBinError( 3, eAcp[1]*percent );

    TGraphAsymmErrors* h1s = new TGraphAsymmErrors(h0);
    h1s->SetPointEYhigh( 0, eAcpTotalUnc[2][0]*percent );
    h1s->SetPointEYlow(  0, eAcpTotalUnc[2][1]*percent );
    h1s->SetPointEYhigh( 1, eAcpTotalUnc[0][0]*percent );
    h1s->SetPointEYlow(  1, eAcpTotalUnc[0][1]*percent );
    h1s->SetPointEYhigh( 2, eAcpTotalUnc[1][0]*percent );
    h1s->SetPointEYlow(  2, eAcpTotalUnc[1][1]*percent );
    h1s->SetFillColor(kGreen);

    TGraphAsymmErrors* h2s = new TGraphAsymmErrors(h0);
    h2s->SetPointEYhigh( 0, 2*eAcpTotalUnc[2][0]*percent );
    h2s->SetPointEYlow(  0, 2*eAcpTotalUnc[2][1]*percent );
    h2s->SetPointEYhigh( 1, 2*eAcpTotalUnc[0][0]*percent );
    h2s->SetPointEYlow(  1, 2*eAcpTotalUnc[0][1]*percent );
    h2s->SetPointEYhigh( 2, 2*eAcpTotalUnc[1][0]*percent );
    h2s->SetPointEYlow(  2, 2*eAcpTotalUnc[1][1]*percent );
    h2s->SetFillColor(kYellow);

    TCanvas *c1 = new TCanvas(("CAcp_"+oName).c_str(), "c1", 41,133,859,639); c1->Clear();
    gStyle->SetOptStat(0);
    c1->Range(-0.5518248,-13.80753,3.192701,11.71548);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1473684);
    c1->SetRightMargin(0.05146199);
    c1->SetTopMargin(0.06721312);
    c1->SetBottomMargin(0.1491803);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
    
    double wrt=1;
    h0->SetMaximum( 0.1*percent*wrt);
    h0->SetMinimum(-0.1*percent*wrt);
    h0->GetXaxis()->SetLabelOffset(0.01);
    h0->GetXaxis()->SetLabelSize(0.12);
    h0->GetXaxis()->SetLabelFont(62);
    h0->GetXaxis()->SetTitleSize(0.035);
    h0->GetYaxis()->SetTitle("A'_{CP} [%]");
    h0->GetYaxis()->SetLabelOffset(0.01);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetTitleSize(0.07);
    h0->GetYaxis()->SetTitleOffset(0.84);
    h0->GetYaxis()->SetTitleFont(42);
    h0->GetZaxis()->SetLabelSize(0.035);
    h0->GetZaxis()->SetTitleSize(0.035);
    h0->Draw();

    h2s->Draw("E2SAME");
    h1s->Draw("E2SAME");

    //h0->SetMarkerStyle(22);
    h0->SetMarkerStyle(8);
    h0->SetMarkerSize(2);
    h0->Draw("psame");

    TLine* line = new TLine(0,0,3,0);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();

    TLegend *leg;
    //if( legX == 0 ) //Left 
    //    leg = new TLegend(0.173516,0.6726768,0.4310502,0.8363384,NULL,"brNDC");
    //else //right
    if( assumpAcp <= 0 )
        leg = new TLegend(0.580117,0.6944444,0.8643275,0.9068627,NULL,"brNDC");
    else
        leg = new TLegend(0.580117,0.1944444,0.8643275,0.4068627,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.06535948);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h0, "1#sigma stat.","pel");
    leg->AddEntry(h1s,"1#sigma stat.+syst.","f");
    leg->AddEntry(h2s,"2#sigma stat.+syst.","f");
    leg->Draw();

    TPaveText* t_title;
    t_title = new TPaveText(0.1134503,0.9393443,0.7426901,0.9852459,"brNDC");
    if( unBlind )
        t_title->AddText("CMS #sqrt{s} = 8TeV, L = 19.7fb^{-1}");
    else
        t_title->AddText("CMS Simulation #sqrt{s} = 8TeV, L = 19.7fb^{-1}");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);
    t_title->Draw();

    if( assumpAcp == 0 ){
        if( unBlind )
            c1->SaveAs((output+"/FinalACP_"+oName+".pdf").c_str());
        else
            c1->SaveAs((output+"/FinalBlindACP_"+oName+".pdf").c_str());
    }else{
        char text[100];
        sprintf( text, "%s/FinalBlindACP_%s_Assumed%.0f.pdf", output.c_str(), oName.c_str(), assumpAcp*100 );
        c1->SaveAs(text);
    }
    delete h0;
    delete [] systName1;
}

void drawSubtractBkgResultsCombined( TFile* f, std::string* systName, int systNi, int nPDFup=0, int nPDFdown=0, float minMlb_=0, float maxMlb_=0, std::string output=".", std::string oName="O2", std::string xTitle="O_{2}", bool unBlind=false, float assumpAcp=0, bool doTopMassUncRescale=true, TFile* fCR=NULL )
{
    std::string tuneName[2]={"up","down"};
    std::string chName[3]={"_El","_Mu",""};
    std::string chFullName[3]={"Electron channel","Muon channel","Combined channel"};
    std::string dataName="DATA";
    if( !unBlind ) dataName="MC";

    bool usingCR = (fCR==NULL)? false:true;
    const int systN = systNi + 1;
    string *systName1 = new string[systN];
    for( int s=0; s<systN; s++){
        if( s != systN-1 ) systName1[s] = systName[s];
        else systName1[s] = "PDF";
    }

    float evtBkgMean[3];
    float evtBkgFitMean[3];
    float bkgWrt[3];
    double  Acp[3];
    double eAcp[3];
    double  AcpSyst[3][2][systN];
    double eAcpSyst[3][2][systN];
    double eAcpSumw2Syst[3][2];
    double eAcpTotalUnc[3][2];
    double  AcpBkg[3];
    double eAcpBkg[3];
    double  Acp0[3];
    double eAcp0[3];

    TH1D* hO_data[3];
    TH1D* hO_bkg[3]; 
    TH1D* hO_bkgEst[3]; 
    TH1D* hT_bkg[3];
    TH1D* hO_dataSyst[3][2][systN];
    TH1D* hO_bkgSyst[3][2][systN];
    TH1D* hT_bkgSyst[3][2][systN];   
    float bkgWrtSyst[3][2][systN];
    for( int ch=0; ch<3; ch++ )
    {
        hO_data[ch] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();
        if( !unBlind )
        {
            // Make fake data base on MC 
            TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
            float total = tmpData->Integral();
            hO_data[ch]->Scale(total/hO_data[ch]->Integral());
            if( assumpAcp != 0 )
            {
                hO_data[ch]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                hO_data[ch]->SetBinContent( 2, (1+assumpAcp)*total/2 );
            }
            hO_data[ch]->SetBinError( 1, sqrt(hO_data[ch]->GetBinContent(1)));
            hO_data[ch]->SetBinError( 2, sqrt(hO_data[ch]->GetBinContent(2)));
            delete tmpData;
        }
        Acp0[ch]  = caculateACP(hO_data[ch]);
        eAcp0[ch] = caculateACPerrorWrt(hO_data[ch]);

        if( ch != 2 ) // Electron or muon channel
        {
            if( usingCR )
            {
                // Subtract O from CR
                std::cout<<"[INFO] Bkg O from data CR..."<<std::endl;
                if( ch == 0 ){      hO_bkg[ch] = (TH1D*)((TH1D*)fCR->Get(("DATA_Electron__EvtChi2_"+oName+"Asym_El").c_str()))->Clone(); }
                else if( ch == 1 ){ hO_bkg[ch] = (TH1D*)((TH1D*)fCR->Get(("DATA_Muon__EvtChi2_"+oName+"Asym_Mu").c_str()))->Clone();     }
            }
            else{ hO_bkg[ch]  = (TH1D*)((TH1D*)f->Get(("BkgMC_"+oName+chName[ch]).c_str()))->Clone(); } // Subtract O from bkg MC

            hT_bkg[ch] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]).c_str()))->Clone();

            int bins = hT_bkg[ch]->GetNbinsX();
            int binW = hT_bkg[ch]->GetBinWidth(1);
            int minMlb, maxMlb;
            if( minMlb_ == maxMlb_ )
            {
                minMlb = 1;
                maxMlb = bins;
            }
            else
            {
                minMlb = minMlb_/binW+1;
                maxMlb = maxMlb_/binW;
            }

            evtBkgMean[ch]  = hO_bkg[ch]->Integral();
            evtBkgFitMean[ch] = hT_bkg[ch]->Integral(minMlb,maxMlb);
            bkgWrt[ch] = evtBkgFitMean[ch]/evtBkgMean[ch];
            hO_bkgEst[ch] = (TH1D*)hO_bkg[ch]->Clone();  
            hO_bkgEst[ch]->Scale(bkgWrt[ch]);   // weight bkg with fitted bkg
            hO_data[ch]->Add(hO_bkgEst[ch],-1); // subtract bkg

             Acp[ch] = caculateACP( hO_data[ch] );
            eAcp[ch] = caculateACPerrorWrt( hO_data[ch] );
             AcpBkg[ch] = caculateACP( hO_bkg[ch] );
            eAcpBkg[ch] = caculateACPerrorWrt( hO_bkg[ch] );
            // * Systmatic uncs.
            int iTopMass=-1;
            for( int s=0; s<systN; s++)
            {
                int pdft[2]={nPDFup, nPDFdown};
                for( int t=0; t<2; t++ ){
                    if( s != systN-1 ){
                        if( systName[s].find("TopMass")!=std::string::npos && doTopMassUncRescale ) iTopMass=s;
                        hT_bkgSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]+"_"+systName[s]+tuneName[t]).c_str()))->Clone();   
                    }else{
                        hT_bkgSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get(("BkgFitted"+chName[ch]+"_PDF"+int2str(pdft[t])).c_str()))->Clone();   
                    }
                    bkgWrtSyst[ch][t][s] = hT_bkgSyst[ch][t][s]->Integral(minMlb,maxMlb)/evtBkgMean[ch]; 
                    if( iTopMass == s ){ bkgWrtSyst[ch][t][s] = (bkgWrtSyst[ch][t][s]+(topMassUncRescale-1)*bkgWrt[ch])/topMassUncRescale; }
                    hO_bkgSyst[ch][t][s] = (TH1D*)hO_bkg[ch]->Clone();
                    hO_bkgSyst[ch][t][s]->Scale(bkgWrtSyst[ch][t][s]); 
                    hO_dataSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();
                    if( !unBlind )
                    {
                        // Make fake data base on MC 
                        TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
                        float total = tmpData->Integral();
                        hO_dataSyst[ch][t][s]->Scale(total/hO_dataSyst[ch][t][s]->Integral());
                        if( assumpAcp != 0 )
                        {
                            hO_dataSyst[ch][t][s]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                            hO_dataSyst[ch][t][s]->SetBinContent( 2, (1+assumpAcp)*total/2 );
                        }
                        delete tmpData;
                    }
                    hO_dataSyst[ch][t][s]->Add(hO_bkgSyst[ch][t][s],-1);
                    AcpSyst[ch][t][s] = caculateACP(hO_dataSyst[ch][t][s]);
                    eAcpSyst[ch][t][s] = (AcpSyst[ch][t][s]-Acp[ch]);
                }
            }
        }
        // Combined channel
        else
        {
            std::cout<<">> Combined channel"<<std::endl;
            hO_bkgEst[ch] = (TH1D*)hO_bkgEst[0]->Clone();
            hO_bkgEst[ch] ->Add(hO_bkgEst[1]);     // Combined electron and muon
            hO_data[ch]->Add(hO_bkgEst[ch],-1); // subtract bkg
             Acp[ch] = caculateACP( hO_data[ch] );
            eAcp[ch] = caculateACPerrorWrt( hO_data[ch] );
             AcpBkg[ch] = caculateACP( hO_bkgEst[ch] );
            eAcpBkg[ch] = caculateACPerrorWrt( hO_bkgEst[ch] );
            // * Systmatic uncs.
            for( int s=0; s<systN; s++)
            {
                for( int t=0; t<2; t++ ){
                    hO_bkgSyst[ch][t][s] = (TH1D*)hO_bkgSyst[0][t][s]->Clone();
                    hO_bkgSyst[ch][t][s] ->Add(hO_bkgSyst[1][t][s]); 
                    hO_dataSyst[ch][t][s] = (TH1D*)((TH1D*)f->Get((dataName+"_"+oName+chName[ch]).c_str()))->Clone();
                    if( !unBlind )
                    {
                        // Make fake data base on MC 
                        TH1D* tmpData = (TH1D*)f->Get(("DATA_"+oName+chName[ch]).c_str());
                        float total = tmpData->Integral();
                        hO_dataSyst[ch][t][s]->Scale(total/hO_dataSyst[ch][t][s]->Integral());
                        if( assumpAcp != 0 )
                        {
                            hO_dataSyst[ch][t][s]->SetBinContent( 1, (1-assumpAcp)*total/2 );
                            hO_dataSyst[ch][t][s]->SetBinContent( 2, (1+assumpAcp)*total/2 );
                        }
                        delete tmpData;
                    }
                    hO_dataSyst[ch][t][s]->Add(hO_bkgSyst[ch][t][s],-1);
                     AcpSyst[ch][t][s] = caculateACP(hO_dataSyst[ch][t][s]);
                    eAcpSyst[ch][t][s] = (AcpSyst[ch][t][s]-Acp[ch]);
                }
            }
        }

        // * Sum Syst.
        eAcpSumw2Syst[ch][0]=0.;
        eAcpSumw2Syst[ch][1]=0.;
        for( int s=0; s<systN; s++ )
        {
            if( eAcpSyst[ch][0][s]*eAcpSyst[ch][1][s] < 0 )
            {
                if( eAcpSyst[ch][0][s] > 0 ){
                    eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                    eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                }else{
                    eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                    eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                } 
            }
            else
            {
                if( fabs(eAcpSyst[ch][0][s]) > fabs(eAcpSyst[ch][1][s]) )
                {
                    if( eAcpSyst[ch][0][s] > 0 ) eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                    if( eAcpSyst[ch][0][s] < 0 ) eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][0][s]*eAcpSyst[ch][0][s];
                }
                else
                {
                    if( eAcpSyst[ch][1][s] > 0 ) eAcpSumw2Syst[ch][0]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                    if( eAcpSyst[ch][1][s] < 0 ) eAcpSumw2Syst[ch][1]+=eAcpSyst[ch][1][s]*eAcpSyst[ch][1][s];
                }
            }
        }
        eAcpTotalUnc[ch][0] = sqrt(eAcpSumw2Syst[ch][0]+eAcp[ch]*eAcp[ch]);
        eAcpTotalUnc[ch][1] = sqrt(eAcpSumw2Syst[ch][1]+eAcp[ch]*eAcp[ch]);

    }

    // * Plots!
    TH1D* h0   = new TH1D(("ACP0_"+oName).c_str(),   "", 3, 0, 3); h0->Sumw2();
    TH1D* hbkg = new TH1D(("ACPbkg_"+oName).c_str(), "", 3, 0, 3); hbkg->Sumw2();
    TH1D* hsig = new TH1D(("ACP_"+oName).c_str(),    "", 3, 0, 3); hsig->Sumw2();

    hbkg->GetXaxis()->SetBinLabel(1, (xTitle+"^{e+#mu}").c_str());
    hbkg->GetXaxis()->SetBinLabel(2, (xTitle+"^{e}").c_str()    );
    hbkg->GetXaxis()->SetBinLabel(3, (xTitle+"^{#mu}").c_str()  );

    double percent=100; 
    h0->Fill( 0., Acp0[2]*percent ); 
    h0->Fill( 1., Acp0[0]*percent );
    h0->Fill( 2., Acp0[1]*percent );
    h0->SetBinError( 1, eAcp0[2]*percent );
    h0->SetBinError( 2, eAcp0[0]*percent );
    h0->SetBinError( 3, eAcp0[1]*percent );
    hbkg->Fill( 0., AcpBkg[2]*percent );
    hbkg->Fill( 1., AcpBkg[0]*percent );
    hbkg->Fill( 2., AcpBkg[1]*percent );
    hbkg->SetBinError( 1, eAcpBkg[2]*percent );
    hbkg->SetBinError( 2, eAcpBkg[0]*percent );
    hbkg->SetBinError( 3, eAcpBkg[1]*percent );
    
    int maxBin=1;
    if( fabs(AcpBkg[2]) < fabs(AcpBkg[0]) ) maxBin=2;
    if( fabs(AcpBkg[0]) < fabs(AcpBkg[1]) ) maxBin=3;
    if( fabs(AcpBkg[2]) > fabs(AcpBkg[1]) && fabs(AcpBkg[2]) > fabs(AcpBkg[0]) ) maxBin=1;

    hsig->Fill( 0., Acp[2]*percent );
    hsig->Fill( 1., Acp[0]*percent );
    hsig->Fill( 2., Acp[1]*percent );
    hsig->SetBinError( 1, eAcp[2]*percent );
    hsig->SetBinError( 2, eAcp[0]*percent );
    hsig->SetBinError( 3, eAcp[1]*percent );

    TGraphAsymmErrors* h1s = new TGraphAsymmErrors(hsig);
    h1s->SetPointEYhigh( 0, eAcpTotalUnc[2][0]*percent );
    h1s->SetPointEYlow(  0, eAcpTotalUnc[2][1]*percent );
    h1s->SetPointEYhigh( 1, eAcpTotalUnc[0][0]*percent );
    h1s->SetPointEYlow(  1, eAcpTotalUnc[0][1]*percent );
    h1s->SetPointEYhigh( 2, eAcpTotalUnc[1][0]*percent );
    h1s->SetPointEYlow(  2, eAcpTotalUnc[1][1]*percent );
    //h1s->SetFillColor(kGreen);

    TCanvas *c1 = new TCanvas(("CAcp_"+oName).c_str(), "c1", 41,133,W,H); c1->Clear();
    gStyle->SetOptStat(0);
    c1->Range(-0.5518248,-13.80753,3.192701,11.71548);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);

    c1->SetLeftMargin(0.1520101);
    c1->SetRightMargin(0.03015076);
    c1->SetTopMargin(0.08726004);
    c1->SetBottomMargin(0.1291448);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    double wrt=0.3;
    TH1D* hbkgLine = (TH1D*)hbkg->Clone();
    hbkg->SetMaximum( 0.1*percent*wrt);
    hbkg->SetMinimum(-0.1*percent*wrt);
    hbkg->GetXaxis()->SetLabelOffset(0.01);
    hbkg->GetXaxis()->SetLabelSize(0.12);
    hbkg->GetXaxis()->SetLabelFont(62);
    hbkg->GetXaxis()->SetTitleSize(0.035);
    hbkg->GetYaxis()->SetTitle("A'_{CP} [%]");
    hbkg->GetYaxis()->SetLabelOffset(0.01);
    hbkg->GetYaxis()->SetLabelSize(0.06);
    hbkg->GetYaxis()->SetTitleSize(0.07);
    hbkg->GetYaxis()->SetTitleOffset(0.84);
    hbkg->GetYaxis()->SetTitleFont(42);
    hbkg->GetZaxis()->SetLabelSize(0.035);
    hbkg->GetZaxis()->SetTitleSize(0.035);
    hbkg->SetLineWidth(2);
    hbkg->SetLineColor(kGreen+2);
    hbkg->SetFillColor(kGreen-8);
    hbkg->SetFillStyle(3008);
    hbkg->Draw("e2");

    hbkgLine->SetLineWidth(2);
    hbkgLine->SetLineColor(kGreen+2);
    hbkgLine->Draw("histsame");

    h0->SetMarkerStyle(4);
    h0->SetMarkerSize(1.5);
    h0->SetMarkerColor(2);
    h0->SetLineColor(2);
    h0->SetLineWidth(2);
    h0->Draw("pe1same");

    h1s->SetLineColor(1);
    h1s->SetLineWidth(2);
    h1s->SetMarkerColor(1);
    h1s->SetMarkerStyle(8);
    h1s->SetMarkerSize(1.5);
    h1s->Draw("PESAME");

    TLine* line = new TLine(0,0,3,0);
   line->SetLineStyle(6);
   line->SetLineWidth(2);
    //line->SetLineColor(2);
    //line->SetLineWidth(3);
    line->Draw();

    TLegend *leg;
    if( assumpAcp == 0 ){
        leg = new TLegend(0.4145729,0.6901571,0.6984925,0.9030716,NULL,"brNDC");
        //if( hbkg->GetBinContent(maxBin) > 0 ) // Only for big error band case
        //    leg = new TLegend(0.1695906,0.1748366,0.4538012,0.3872549,NULL,"brNDC");
    }else if( assumpAcp < 0 ){
        leg = new TLegend(0.4145729,0.6701571,0.6984925,0.8830716,NULL,"brNDC");
    }else{
        leg = new TLegend(0.1695906,0.1748366,0.4538012,0.3872549,NULL,"brNDC");
    }
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04901961);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1s,  "t#bar{t} events #pm 1#sigma (stat.+syst.)","pel");
    leg->AddEntry(h0,   "Before background subtraction",       "pel");
    leg->AddEntry(hbkg, "Estimated background",                "fl");
    leg->Draw();

    //TPaveText* t_title;
    //t_title = new TPaveText(0.1134503,0.9393443,0.7426901,0.9852459,"brNDC");
    //if( unBlind )
    //    t_title->AddText("CMS #sqrt{s} = 8TeV, L = 19.7fb^{-1}");
    //else
    //    t_title->AddText("CMS Simulation #sqrt{s} = 8TeV, L = 19.7fb^{-1}");
    //t_title->SetTextColor(kBlack);
    //t_title->SetFillColor(kWhite);
    //t_title->SetFillStyle(0);
    //t_title->SetBorderSize(0);
    //t_title->SetTextAlign(12);
    //t_title->SetTextSize(0.04);
    //t_title->Draw();

    writeExtraText=true;
    if( !unBlind ) extraText="Simulation";

    writeExtraText=true;
    CMS_lumi( c1, 2, cmslumi, true );

    if( assumpAcp == 0 ){
        if( unBlind )
            c1->SaveAs((output+"/FinalACPSum_"+oName+".pdf").c_str());
        else
            c1->SaveAs((output+"/FinalBlindACPSum_"+oName+".pdf").c_str());
    }else{
        char text[100];
        sprintf( text, "%s/FinalBlindACPSum_%s_Assumed%.0f.pdf", output.c_str(), oName.c_str(), assumpAcp*100 );
        c1->SaveAs(text);
    }
    delete h0;
    delete hsig;
    delete hbkg;
}
