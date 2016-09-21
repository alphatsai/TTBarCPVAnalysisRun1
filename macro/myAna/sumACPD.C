#include "/Users/Alpha/myAna/TTBarCPV/TTBarCPVAnalysisRun1/macro/CMS_lumi.C"
void sumACPD()
{
    float ACP2[3][3] = {{-0.19 , 0.61 , 0.59}, { 0.46 , 0.57 , 0.65}, { 0.16 , 0.42 , 0.44}}; 
    float ACP3[3][3] = {{ 0.02 , 0.61 , 0.59}, {-0.59 , 0.57 , 0.65}, {-0.31 , 0.42 , 0.44}};  
    float ACP4[3][3] = {{-0.17 , 0.61 , 0.59}, {-0.10 , 0.57 , 0.65}, {-0.13 , 0.42 , 0.44}};
    float ACP7[3][3] = {{-0.38 , 0.61 , 0.59}, { 0.43 , 0.57 , 0.65}, { 0.06 , 0.42 , 0.44}};

    float D2 = 0.575;
    float D3 = 0.385;
    float D4 = 0.369;
    float D7 = 0.730;

    float cACP2[3][2]; 
    float cACP3[3][2];  
    float cACP4[3][2];
    float cACP7[3][2];

    for( int i=0; i<3; i++)
    {
            cACP2[i][0]=ACP2[i][0]/D2;
            cACP3[i][0]=ACP3[i][0]/D3;
            cACP4[i][0]=ACP4[i][0]/D4;
            cACP7[i][0]=ACP7[i][0]/D7;
        
            cACP2[i][1]=sqrt(ACP2[i][1]*ACP2[i][1]+ACP2[i][2]*ACP2[i][2])/D2;
            cACP3[i][1]=sqrt(ACP3[i][1]*ACP3[i][1]+ACP3[i][2]*ACP3[i][2])/D3;
            cACP4[i][1]=sqrt(ACP4[i][1]*ACP4[i][1]+ACP4[i][2]*ACP4[i][2])/D4;
            cACP7[i][1]=sqrt(ACP7[i][1]*ACP7[i][1]+ACP7[i][2]*ACP7[i][2])/D7;
    }

    printf("%12s &", "$\\Otwo$");
    for( int i=0; i<3; i++){ printf("$%+.2f \\pm %.2f $ &", cACP2[i][0], cACP2[i][1] );}
    printf("\n%12s &", "$\\Othree$");
    for( int i=0; i<3; i++){ printf("$%+.2f \\pm %.2f $ &", cACP3[i][0], cACP3[i][1] );}
    printf("\n%12s &", "$\\Ofour$");
    for( int i=0; i<3; i++){ printf("$%+.2f \\pm %.2f $ &", cACP4[i][0], cACP4[i][1] );}
    printf("\n%12s &", "$\\Oseven$");
    for( int i=0; i<3; i++){ printf("$%+.2f \\pm %.2f $ &", cACP7[i][0], cACP7[i][1] );}
    printf("\n");

    TH1D* hall = new TH1D("All", "", 12, 0, 12 ); hall->Sumw2();

    for( int i=0; i<3; i++)
    {
        hall->GetXaxis()->SetBinLabel(1+4*i, "O_{2}");
        hall->GetXaxis()->SetBinLabel(2+4*i, "O_{3}");
        hall->GetXaxis()->SetBinLabel(3+4*i, "O_{4}");
        hall->GetXaxis()->SetBinLabel(4+4*i, "O_{7}");
    
        hall->SetBinContent(1+4*i, cACP2[i][0]);
        hall->SetBinContent(2+4*i, cACP3[i][0]);
        hall->SetBinContent(3+4*i, cACP4[i][0]);
        hall->SetBinContent(4+4*i, cACP7[i][0]);
        hall->SetBinError(1+4*i, cACP2[i][1]);
        hall->SetBinError(2+4*i, cACP3[i][1]);
        hall->SetBinError(3+4*i, cACP4[i][1]);
        hall->SetBinError(4+4*i, cACP7[i][1]);
    }

    TCanvas *c1 = new TCanvas("CSumAcp", "c1", 41,133,800,600); c1->Clear();
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

    float width=4;
    hall->SetMaximum(    width );
    hall->SetMinimum( -1*width );
    hall->GetXaxis()->SetLabelOffset(0.01);
    hall->GetXaxis()->SetLabelSize(0.07);
    hall->GetXaxis()->SetLabelFont(62);
    hall->GetXaxis()->SetTitleSize(0.035);
    hall->GetYaxis()->SetTitle("A'_{CP} [%]");
    hall->GetYaxis()->SetLabelOffset(0.01);
    hall->GetYaxis()->SetLabelSize(0.06);
    hall->GetYaxis()->SetTitleSize(0.07);
    hall->GetYaxis()->SetTitleOffset(0.84);
    hall->GetYaxis()->SetTitleFont(42);
    hall->GetZaxis()->SetLabelSize(0.035);
    hall->GetZaxis()->SetTitleSize(0.035);
    hall->SetMarkerStyle(4);
    hall->SetMarkerSize(1.5);
    hall->SetMarkerColor(2);
    hall->SetLineColor(2);
    hall->SetLineWidth(2);
    hall->Draw("pe");

    const int Osize=4;
    const int Bins=Osize*3;
    TLine* line0 = new TLine(0,0,Bins,0);
    line0->SetLineStyle(6);
    line0->SetLineWidth(2);
    line0->Draw();

    TLine* line1 = new TLine(Osize,-0.625*width,Osize,0.625*width);
    line1->SetLineStyle(6);
    line1->SetLineColor(38);
    line1->SetLineWidth(2);
    line1->Draw();

    TLine* line2 = new TLine(Osize*2,-0.625*width,Osize*2,0.625*width);
    line2->SetLineStyle(6);
    line2->SetLineColor(38);
    line2->SetLineWidth(2);
    line2->Draw();

    TPaveText* chEl, *chMu, *chCo;
    chEl = new TPaveText(1,         -1*width+0.5, 3,         -1*width+1, "br");
    chMu = new TPaveText(1+Osize,   -1*width+0.5, 3+Osize,   -1*width+1, "br");
    chCo = new TPaveText(1+2*Osize, -1*width+0.5, 3+2*Osize, -1*width+1, "br");
    chEl->SetTextColor(kBlack);
    chMu->SetTextColor(kBlack);
    chCo->SetTextColor(kBlack);
    chEl->SetFillColor(kWhite);
    chMu->SetFillColor(kWhite);
    chCo->SetFillColor(kWhite);
    chEl->SetFillStyle(0);
    chMu->SetFillStyle(0);
    chCo->SetFillStyle(0);
    chEl->SetBorderSize(0);
    chMu->SetBorderSize(0);
    chCo->SetBorderSize(0);
    chEl->SetTextAlign(21);
    chMu->SetTextAlign(21);
    chCo->SetTextAlign(21);
    chEl->SetTextFont(42);
    chMu->SetTextFont(42);
    chCo->SetTextFont(42);
    chEl->SetTextSize(0.05235602);
    chMu->SetTextSize(0.05235602);
    chCo->SetTextSize(0.05235602);
    chEl->AddText("e+jets");
    chMu->AddText("#mu+jets");
    chCo->AddText("l+jets");
    chEl->Draw();
    chMu->Draw();
    chCo->Draw();

    TLegend *leg;
    leg = new TLegend(0.4924623,0.6858639,0.7763819,0.8987784,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04363002);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hall, "Applied correction",            "pel");
    leg->Draw();

    writeExtraText=false;
    CMS_lumi( c1, 2, 10, writeExtraText );
    c1->SaveAs(("fig/FinalACPD.pdf").c_str());
    c1->SaveAs(("fig/root/FinalACPD.root").c_str());



}


