void integralFromLowerBins( TH1D* h_in, TH1D* h_out ) 
{
   if( h_in->GetNbinsX() != h_out->GetNbinsX() ) printf(">> [ERROR] Deffirent bin size between h_in(%d) and h_out(%d)\n", h_in->GetNbinsX(), h_out->GetNbinsX()); 
    int minBin = 1;
    int maxBin = h_in->GetNbinsX();
    double sum   = 0;
    double sumw2 = 0;
    for( int b=minBin; b<=maxBin; b++ )
    {
        sum   += h_in->GetBinContent(b);
        sumw2 += h_in->GetBinError(b)*h_in->GetBinError(b);
        h_out->SetBinContent(b, sum);
        h_out->SetBinError(b, sqrt(sumw2));
    }
}
void drawBBbarTagEff( TFile* f, std::string savePath, int rebin=1, std::string ana="bbarTagEff" )
{
    TH1D* h_topHadronChi2 = (TH1D*)f->Get((ana+"/Evt_Top_Hadronic_Chi2").c_str()); 
    TH1D* h_bbMatchedChi2 = (TH1D*)f->Get((ana+"/Evt_bMatched_Chi2").c_str()); 
    h_topHadronChi2->Rebin(rebin);
    h_bbMatchedChi2->Rebin(rebin);

    TH1D* h_topHadronChi2Int = (TH1D*)h_topHadronChi2->Clone("h_topHadronChi2Int");
    TH1D* h_bbMatchedChi2Int = (TH1D*)h_bbMatchedChi2->Clone("h_bbMatchedChi2Int");
    integralFromLowerBins( h_topHadronChi2, h_topHadronChi2Int ); 
    integralFromLowerBins( h_bbMatchedChi2, h_bbMatchedChi2Int ); 
   
    double totalEvts = h_topHadronChi2->Integral();
    double nomalized = 1/totalEvts; 
    TH1D* h_topHadronChi2Eff = (TH1D*)h_topHadronChi2Int->Clone("h_topHadronChi2Eff"); 
    TH1D* h_bbMatchedChi2Eff = (TH1D*)h_bbMatchedChi2Int->Clone("h_bbMatchedChi2Eff"); 
    TH1D* h_bbTotalChi2Eff   = (TH1D*)h_bbMatchedChi2Int->Clone("h_bbMatchedChi2Eff");

    h_topHadronChi2Eff->Scale(nomalized); 
    h_bbMatchedChi2Eff->Divide(h_topHadronChi2Int); 
    h_bbTotalChi2Eff->Scale(nomalized); // The same h_topHadronChi2Eff->Multiply(h_bbMatchedChi2Eff);

    h_topHadronChi2Eff->GetXaxis()->SetTitle("#chi^{2} cut"); 
    h_topHadronChi2Eff->GetXaxis()->SetLabelFont(42);
    h_topHadronChi2Eff->GetXaxis()->SetLabelSize(0.05);
    h_topHadronChi2Eff->GetXaxis()->SetTitleSize(0.07);
    h_topHadronChi2Eff->GetXaxis()->SetTitleOffset(0.83);
    h_topHadronChi2Eff->GetXaxis()->SetTitleFont(42);
    h_topHadronChi2Eff->GetYaxis()->SetTitle("Efficiency");
    h_topHadronChi2Eff->GetYaxis()->SetLabelFont(42);
    h_topHadronChi2Eff->GetYaxis()->SetLabelSize(0.05);
    h_topHadronChi2Eff->GetYaxis()->SetNdivisions(410);
    h_topHadronChi2Eff->GetYaxis()->SetTitleSize(0.06);
    h_topHadronChi2Eff->GetYaxis()->SetTickLength(0.02);
    h_topHadronChi2Eff->GetYaxis()->SetTitleOffset(0.83);
    h_topHadronChi2Eff->GetYaxis()->SetTitleFont(42);
    h_topHadronChi2Eff->SetMinimum(0);
    h_topHadronChi2Eff->SetLineColor(2); 
    h_topHadronChi2Eff->SetLineWidth(3); 
    h_bbMatchedChi2Eff->SetLineColor(4); 
    h_bbMatchedChi2Eff->SetLineWidth(3); 
    h_bbTotalChi2Eff->SetLineColor(1); 
    h_bbTotalChi2Eff->SetLineWidth(3); 
    h_bbTotalChi2Eff->SetLineStyle(7); 

    TCanvas *c1 = new TCanvas("c1_drawBBbarTagEff", "",1443,136,956,719);
    c1->Clear();
    gStyle->SetOptStat(0);
    c1->Range(-30.12658,0.01462481,210.8861,1.099256);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.125);
    c1->SetRightMargin(0.04516807);
    c1->SetTopMargin(0.05362319);
    c1->SetBottomMargin(0.1507246);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
   
    h_topHadronChi2Eff->Draw("HIST"); 
    h_bbMatchedChi2Eff->Draw("HISTSAME"); 
    h_bbTotalChi2Eff->Draw("HISTSAME");

    TLegend *leg = new TLegend(0.6186975,0.2217391,0.9978992,0.3985507,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry( h_topHadronChi2Eff, "#chi^{2} cut eff.",   "l");
    leg->AddEntry( h_bbMatchedChi2Eff, "b#bar{b}-tagging eff.","l");
    leg->AddEntry( h_bbTotalChi2Eff,   "Total eff.","l");
    leg->Draw();

    c1->SaveAs((savePath+"/bbbartagEff.pdf").c_str());

    delete h_topHadronChi2;
    delete h_bbMatchedChi2;
}
void drawBBbarTagPDF( TFile* f, std::string savePath, bool log=true, int rebin=1, int xMax=200, std::string ana="bbarTagEff" )
{
    TH1D* h_bbMatchedChi2     = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2"      ).c_str()))->Clone("h_bbMatchedChi2"    ); h_bbMatchedChi2->Rebin(rebin);
    TH1D* h_missTagChi2_lbbar = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_lbbar").c_str()))->Clone("h_missTagChi2_lbbar"); h_missTagChi2_lbbar->Rebin(rebin);
    TH1D* h_missTagChi2_bl    = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_bl"   ).c_str()))->Clone("h_missTagChi2_bl"   ); h_missTagChi2_bl->Rebin(rebin);
    TH1D* h_missTagChi2_bbarb = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_bbarb").c_str()))->Clone("h_missTagChi2_bbarb"); h_missTagChi2_bbarb->Rebin(rebin);
    TH1D* h_missTagChi2_bbarl = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_bbarl").c_str()))->Clone("h_missTagChi2_bbarl"); h_missTagChi2_bbarl->Rebin(rebin);
    TH1D* h_missTagChi2_lb    = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_lb"   ).c_str()))->Clone("h_missTagChi2_lb"   ); h_missTagChi2_lb->Rebin(rebin);
    TH1D* h_missTagChi2_ll    = (TH1D*)((TH1D*)f->Get((ana+"/Evt_bMatched_Chi2_ll"   ).c_str()))->Clone("h_missTagChi2_ll"   ); h_missTagChi2_ll->Rebin(rebin);

    TH1D* h_oneTagChi2 = (TH1D*)h_missTagChi2_lbbar->Clone("h_oneTagChi2"); 
    h_oneTagChi2->Add(h_missTagChi2_bl);

    TH1D* h_invTagChi2 = (TH1D*)h_missTagChi2_bbarb->Clone("h_invTagChi2");
    h_invTagChi2->Add(h_missTagChi2_bbarl);
    h_invTagChi2->Add(h_missTagChi2_lb);

    double n_bb    = h_bbMatchedChi2->GetEntries();
    double n_lbbar = h_missTagChi2_lbbar->GetEntries();
    double n_bl    = h_missTagChi2_bl->GetEntries();
    double n_bbarb = h_missTagChi2_bbarb->GetEntries();
    double n_bbarl = h_missTagChi2_bbarl->GetEntries();
    double n_lb    = h_missTagChi2_lb->GetEntries();
    double n_ll    = h_missTagChi2_ll->GetEntries();

    double evts = n_bb + n_lbbar + n_bl + n_bbarb + n_bbarl + n_lb + n_ll;
 
    h_bbMatchedChi2->GetXaxis()->SetTitle("#chi^{2}"); 
    h_bbMatchedChi2->GetXaxis()->SetLabelFont(42);
    h_bbMatchedChi2->GetXaxis()->SetLabelSize(0.05);
    h_bbMatchedChi2->GetXaxis()->SetTitleSize(0.07);
    h_bbMatchedChi2->GetXaxis()->SetTitleOffset(0.83);
    h_bbMatchedChi2->GetXaxis()->SetTitleFont(42);
    h_bbMatchedChi2->GetXaxis()->SetRangeUser(0,xMax);
    h_bbMatchedChi2->GetYaxis()->SetTitle("PDF");
    h_bbMatchedChi2->GetYaxis()->SetLabelFont(42);
    h_bbMatchedChi2->GetYaxis()->SetLabelSize(0.05);
    h_bbMatchedChi2->GetYaxis()->SetNdivisions(410);
    h_bbMatchedChi2->GetYaxis()->SetTitleSize(0.06);
    h_bbMatchedChi2->GetYaxis()->SetTickLength(0.02);
    h_bbMatchedChi2->GetYaxis()->SetTitleOffset(0.83);
    h_bbMatchedChi2->GetYaxis()->SetTitleFont(42);
    h_bbMatchedChi2->SetLineColor(2);   
    h_bbMatchedChi2->SetLineWidth(3); 
    h_oneTagChi2->SetLineColor(4); 
    h_oneTagChi2->SetLineWidth(3); 
    h_invTagChi2->SetLineColor(8); 
    h_invTagChi2->SetLineWidth(3); 
    h_invTagChi2->SetLineStyle(7); 
    h_missTagChi2_ll->SetLineColor(1); 
    h_missTagChi2_ll->SetLineWidth(3); 
    h_missTagChi2_ll->SetLineStyle(7);
 
    TCanvas *c1 = new TCanvas("c1_drawBBbarTagPDF", "",1443,136,956,719);
    c1->Clear();
    gStyle->SetOptStat(0);
    c1->Range(-30.12658,0.01462481,210.8861,1.099256);
    //c1->SetGridx();
    //c1->SetGridy();
    if( log ) c1->SetLogy();
    else c1->SetLogy(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.125);
    c1->SetRightMargin(0.04516807);
    c1->SetTopMargin(0.05362319);
    c1->SetBottomMargin(0.1507246);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    h_bbMatchedChi2->DrawNormalized("HIST"); 
    h_oneTagChi2->DrawNormalized("HISTSAME"); 
    h_invTagChi2->DrawNormalized("HISTSAME"); 
    h_missTagChi2_ll->DrawNormalized("HISTSAME"); 

    TLegend *leg = new TLegend(0.4516807,0.6676301,0.8308824,0.9089595,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04913295);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry( h_bbMatchedChi2,  "Matched b#bar{b}",              "l");
    leg->AddEntry( h_oneTagChi2,     "Matched b or #bar{b}",          "l");
    leg->AddEntry( h_invTagChi2,     "Inversed matching b or #bar{b}","l");
    leg->AddEntry( h_missTagChi2_ll, "Both matched light quarks",     "l");
    leg->Draw();

    c1->SaveAs((savePath+"/bbbartagPDF.pdf").c_str());   

    printf("- Total events : %.0d\n",     evts             );
    printf("- Gen-Matching :  Eff[%s]\n", "%"              );
    printf(" Matched b  b~ : %5.2f\n",    100*n_bb/evts    );
    printf(""                                              );
    printf(" Matched lq b~ : %5.2f\n",    100*n_lbbar/evts );
    printf(" Matched b  lq : %5.2f\n",    100*n_bl/evts    );
    printf(""                                              );
    printf(" Matched b~ b  : %5.2f\n",    100*n_bbarb/evts );
    printf(" Matched b~ lq : %5.2f\n",    100*n_bbarl/evts );
    printf(" Matched lq b  : %5.2f\n",    100*n_lb/evts    );
    printf(" Matched lq lq : %5.2f\n",    100*n_ll/evts    );
    printf(">> Save to %s\n", (savePath+"/bbbartagGenMatchingTable.txt").c_str());

    FILE * out = fopen((savePath+"/bbbartagGenMatchingTable.txt").c_str(),"w");
    fprintf( out, "- Total events : %.0d\n",     evts             );
    fprintf( out, "- Gen-Matching :  Eff[%s]\n", "%"              );
    fprintf( out, " Matched b  b~ : %5.2f\n",    100*n_bb/evts    );
    fprintf( out, ""                                              );
    fprintf( out, " Matched lq b~ : %5.2f\n",    100*n_lbbar/evts );
    fprintf( out, " Matched b  lq : %5.2f\n",    100*n_bl/evts    );
    fprintf( out, ""                                              );
    fprintf( out, " Matched b~ b  : %5.2f\n",    100*n_bbarb/evts );
    fprintf( out, " Matched b~ lq : %5.2f\n",    100*n_bbarl/evts );
    fprintf( out, " Matched lq b  : %5.2f\n",    100*n_lb/evts    );
    fprintf( out, " Matched lq lq : %5.2f\n",    100*n_ll/evts    );
    fclose(out);

    double n_bb    = h_bbMatchedChi2->GetEntries();
    double n_lbbar = h_missTagChi2_lbbar->GetEntries();
    double n_bl    = h_missTagChi2_bl->GetEntries();
    double n_bbarb = h_missTagChi2_bbarb->GetEntries();
    double n_bbarl = h_missTagChi2_bbarl->GetEntries();
    double n_lb    = h_missTagChi2_lb->GetEntries();
    double n_ll    = h_missTagChi2_ll->GetEntries();
}
