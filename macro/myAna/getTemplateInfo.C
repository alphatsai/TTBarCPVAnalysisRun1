#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TColor.h"
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

void drawFittedStack( TFile* f, std::string hName, std::string* systName, int systNi, int chN=0, std::string output=".", std::string xTitle="", std::string yTitle="Events", bool logy=false )
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

    const int systN = systNi;
    TH1D* hs;
    TH1D* h_data;
    TH1D* h_bkg;
    TH1D* h_sig;
    TH1D* h_all;
    TH1D* h_allUnc[systN][2];
    double uncSig[systN][2]; 
    double uncBkg[systN][2]; 
    double uncAll[systN][2];
    double meanSig; 
    double meanBkg; 
    double meanAll; 

    h_data = (TH1D*)((TH1D*)f->Get((("DATA"+ch).c_str())))->Clone("DATA");
    h_sig  = (TH1D*)((TH1D*)f->Get((("SigFitted"+ch).c_str())))->Clone("SIG"); meanSig = h_sig->Integral();
    h_bkg  = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch).c_str())))->Clone("BKG"); meanBkg = h_bkg->Integral();
    h_all  = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch).c_str())))->Clone("ALL"); 
    h_all->Add(h_sig);
    meanAll = h_all->Integral();

    std::string tune[2] = {"up", "down"};
    for( int i=0; i<systN; i++ )
    {   
        //cout<<systName[i]<<endl;
        TH1D* h_sigUnc[2];
        for( int t=0; t<2; t++)
        {
            std::string name = systName[i]+tune[t];
            cout<<name<<endl;
            h_sigUnc[t]    = (TH1D*)f->Get((("SigFitted"+ch+"_"+name).c_str()));
            h_allUnc[i][t] = (TH1D*)((TH1D*)f->Get((("BkgFitted"+ch+"_"+name).c_str())))->Clone(("UNC_"+name).c_str());
            uncSig[i][t]=h_sigUnc[t]->Integral();
            uncBkg[i][t]=h_allUnc[i][t]->Integral(); h_allUnc[i][t]->Add(h_sigUnc[t]);
            uncAll[i][t]=h_allUnc[i][t]->Integral();
        }
        delete *h_sigUnc;
    }

    TGraphAsymmErrors* h_allAsymErr = new TGraphAsymmErrors(h_all);
    int bins = h_all->GetNbinsX();
    for( int b=1; b<=bins; b++ )
    {
        float nomV = h_all->GetBinContent(b);
        float sigUp   = 0;
        float sigDown = 0;
        for( int i=0; i<systN; i++ ){
            for( int t=0; t<2; t++ )
            {
                float err = h_allUnc[i][t]->GetBinContent(b) - nomV;
                if( err > 0. ) sigUp +=err*err;
                else sigDown += err*err;
            }
        }
        //cout<<b<<": "<<sqrt(sigUp)<<" "<<sqrt(sigDown)<<endl;
        h_allAsymErr->SetPointEYhigh( b, sqrt(sigUp));
        h_allAsymErr->SetPointEYlow(  b, sqrt(sigDown));
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
        c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",59,72,1076,826);
    }else{
        c1 = new TCanvas( ("C_Linear_"+hName).c_str(), "",59,72,1076,826);
    }
    c1->Range(-123.0964,-841.8412,557.1066,3883.14);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1809701);
    c1->SetRightMargin(0.08395522);
    c1->SetTopMargin(0.06148055);
    c1->SetBottomMargin(0.1781681);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if( logy ) c1->SetLogy(1);
    else c1->SetLogy(0);

    TLegend *leg;
    leg = new TLegend(0.6091418,0.6398996,0.9085821,0.9259724,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "Data", "lpe");
    leg->AddEntry(h_sig,  "Estimated sig.", "f");
    leg->AddEntry(h_bkg,  "Estimated bkg.", "f");
    leg->AddEntry(h_allAsymErr, "1#sigma, Stat.+Syst.", "f");

    TPaveText* t_title;
    t_title = new TPaveText(0.1455224,0.9134253,0.7751866,0.9974906,"brNDC");
    //t_title = new TPaveText(0.07088487,0.9153846,0.7007612,1,"brNDC");
    t_title->AddText("CMS #sqrt{s} = 8TeV, L = 19.7/fb");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(11);
    t_title->SetTextSize(0.04805273);

    h_stack->Draw("HIST");
    h_allAsymErr->Draw("E2SAME");
    h_data->Draw("ESAME");
    leg->Draw();
    t_title->Draw();

    if( logy )
        c1->SaveAs((output+"/StackFitted_Log_"+hName+ch+".pdf").c_str());
    else
        c1->SaveAs((output+"/StackFitted_Linear_"+hName+ch+".pdf").c_str());

    float sumw2[2]={0,0}, sumw2sig[2]={0,0}, sumw2bkg[2]={0,0};
    FILE* outTxt;
    outTxt = fopen((output+"/SystUncertainties_"+hName+ch+".txt").c_str(),"w");
    fprintf( outTxt, "%s\n", channel.c_str() );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<systNi; i++ ){ fprintf( outTxt, "%9s", systName[i].c_str()); }   
    fprintf( outTxt, "\nTotal:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncAll[i][0]-meanAll)/meanAll*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2[0]+=v*v; } }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncAll[i][1]-meanAll)/meanAll*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2[1]+=v*v; } }  
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2[1]) );
    fprintf( outTxt, "\n");
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncSig[i][0]-meanSig)/meanSig*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sig[0]+=v*v; } }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2sig[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncSig[i][1]-meanSig)/meanSig*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sig[1]+=v*v; } }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2sig[1]) );
    fprintf( outTxt, "\n");
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncBkg[i][0]-meanBkg)/meanBkg*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkg[0]+=v*v; }}   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2bkg[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<systNi; i++ ){ float v=(uncBkg[i][1]-meanBkg)/meanBkg*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkg[1]+=v*v; }}   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2bkg[1]) );
    fprintf( outTxt, "\n\n");
    fclose( outTxt );
}
