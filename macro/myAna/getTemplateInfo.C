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
    sprintf( text, "Nominal, %5d events", values[0] ); cout<<text<<endl;
    leg->AddEntry(h_mean, text, "l");
    sprintf( text, "+1 #sigma, %5d events", values[1] ); cout<<text<<endl;
    leg->AddEntry(h_up,   text, "l");
    sprintf( text, "-1 #sigma, %5d events", values[2] ); cout<<text<<endl;
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
