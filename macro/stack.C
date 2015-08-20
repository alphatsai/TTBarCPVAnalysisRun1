#include "help.C"
void drawStack( TFile* f, std::string hName, std::string xtitle="", std::string ytitle="Events", std::string output=".", int rebin=1, bool logy=false ){

    int lineWidth=3;
    Color_t c_tt    = kOrange+1;
    Color_t c_ttbkg = kOrange+2;
    Color_t c_t     = kSpring+9;
    Color_t c_b     = kAzure+2;
    Color_t c_unc   = kRed-7;
    Color_t c_allunc = kRed-1;

    TH1D *h_tt, *h_ttbkg, *h_t, *h_b, *h_bkg, *h_all, *hs;
    h_tt     = (TH1D*)((TH1D*)f->Get(("TTJets_SemiLeptMGDecays__"+hName).c_str()))->Clone("TTBarSemiLept");
    h_ttbkg  = (TH1D*)((TH1D*)f->Get(("TTJets_NonSemiLeptMGDecays__"+hName).c_str()))->Clone("TTBarNonSemiLept");
    h_t      = (TH1D*)((TH1D*)f->Get(("SingleT__"+hName).c_str()))->Clone("SingleTop");
    h_b      = (TH1D*)((TH1D*)f->Get(("Boson__"+hName).c_str()))->Clone("Boson");
    h_bkg    = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__"+hName).c_str()))->Clone("BkgUnc");
    //h_bkg = (TH1D*)((TH1D*)f->Get(("BkgMC__"+hName).c_str()))->Clone("BkgUnc");
    h_all    = (TH1D*)((TH1D*)f->Get(("MC__"+hName).c_str()))->Clone("AllUnc");

    fix(h_tt);    h_tt->Rebin(rebin);
    fix(h_ttbkg); h_ttbkg->Rebin(rebin);
    fix(h_t);     h_t->Rebin(rebin);
    fix(h_b);     h_b->Rebin(rebin);
    fix(h_bkg);   h_bkg->Rebin(rebin);
    fix(h_all);   h_all->Rebin(rebin);

    int bins = h_all->GetXaxis()->GetLast();
    int xMin = h_all->GetXaxis()->GetBinLowEdge(1);
    int xMax = h_all->GetXaxis()->GetBinUpEdge(bins);

    h_tt->SetLineWidth(lineWidth);
    h_tt->SetLineColor(c_tt);
    h_tt->SetFillColor(c_tt);

    h_ttbkg->SetLineWidth(lineWidth);
    h_ttbkg->SetLineColor(c_ttbkg);
    h_ttbkg->SetFillColor(c_ttbkg);

    h_t->SetLineWidth(lineWidth);
    h_t->SetLineColor(c_t);
    h_t->SetFillColor(c_t);

    h_b->SetLineWidth(lineWidth);
    h_b->SetLineColor(c_b);
    h_b->SetFillColor(c_b);

    h_bkg->SetFillStyle(3244);
    h_bkg->SetFillColor(c_unc);

    h_all->SetFillStyle(3244);
    h_all->SetFillColor(c_allunc);

    if( logy )
        hs = new TH1D(("TH1DinStackLog"+hName).c_str(), "", bins, xMin, xMax);
    else
        hs = new TH1D(("TH1DinStackLinear"+hName).c_str(), "", bins, xMin, xMax);

    hs->SetMinimum(0.8995655);
    hs->GetXaxis()->SetTitle(xtitle.c_str());
    hs->GetYaxis()->SetTitle(ytitle.c_str());

    hs->GetXaxis()->SetLabelFont(42);
    hs->GetXaxis()->SetLabelSize(0.05);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetXaxis()->SetTitleOffset(1.04);
    hs->GetXaxis()->SetTitleFont(42);
    hs->GetYaxis()->SetLabelFont(42);
    hs->GetYaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetTitleFont(42);
    hs->GetZaxis()->SetLabelFont(42);
    hs->GetZaxis()->SetLabelSize(0.035);
    hs->GetZaxis()->SetTitleSize(0.035);
    hs->GetZaxis()->SetTitleFont(42);

    THStack* h_stack = new THStack("THStcak", "");
    h_stack->SetHistogram(hs);

    h_stack->Add(h_b);
    h_stack->Add(h_t);
    h_stack->Add(h_ttbkg);
    h_stack->Add(h_tt);

    if( logy ){
        if( h_all->GetMinimum() > 0 ) 
            h_stack->SetMinimum(h_all->GetMinimum()/10);
        else
            h_stack->SetMinimum(10);
    }

    TCanvas* c1;
    if( logy ) 
        c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",261,48,1179,808);
    else
        c1 = new TCanvas( ("C_Linear_"+hName).c_str(), "",261,48,1179,808);
    c1->Range(-2.617021,-1.006188,2.382979,5.162451);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1234043);
    c1->SetRightMargin(0.07659575);
    c1->SetTopMargin(0.06209987);
    c1->SetBottomMargin(0.1376441);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0); 

    if( logy ) c1->SetLogy(1);
    else c1->SetLogy(0);
    TLegend *leg;
    leg = new TLegend(0.6561702,0.7067862,0.972766,0.9116517,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_tt, "t#bar{t}+jet (lepton+jet)", "f");
    leg->AddEntry(h_ttbkg, "t#bar{t}+jet (Other)", "f");
    leg->AddEntry(h_t, "Single top", "f");
    leg->AddEntry(h_b, "Z/#gamma*/W/WW/WZ/ZZ", "f");
    leg->AddEntry(h_bkg, "1#sigma non t#bar{t}+jet stat.", "f");
    leg->AddEntry(h_all, "1#sigma Total stat.", "f");

    TPaveText* t_title;
    t_title = new TPaveText(0.0893617,0.9443022,0.7191489,0.9891165,"brNDC");
    t_title->AddText("CMS Simulation, L = 19.7/fb, #sqrt{s} = 8TeV");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);

    h_stack->Draw("HIST");
    h_all->Draw("SAMEE2");
    h_bkg->Draw("SAMEE2");
    leg->Draw();
    t_title->Draw();

    if( logy )
        c1->SaveAs((output+"/Stack_Log_"+hName+".pdf").c_str());
    else
        c1->SaveAs((output+"/Stack_Linear_"+hName+".pdf").c_str());
}
