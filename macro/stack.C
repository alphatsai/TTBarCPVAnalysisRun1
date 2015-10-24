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
    h_bkg    = (TH1D*)((TH1D*)f->Get(("BkgMC__"+hName).c_str()))->Clone("BkgUnc");
    //h_bkg = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded_"+hName).c_str()))->Clone("BkgUnc");
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
    //hs->GetXaxis()->SetTitleOffset(1.04);
    hs->GetXaxis()->SetTitleOffset(0.98);
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
void drawStackWithData( TFile* f, std::string hName, std::string xtitle="", std::string ytitle="Events", std::string output=".", int rebin=1, bool logy=false, bool unity=false ){

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
    h_bkg    = (TH1D*)((TH1D*)f->Get(("BkgMC__"+hName).c_str()))->Clone("BkgUnc");
    h_all    = (TH1D*)((TH1D*)f->Get(("MC__"+hName).c_str()))->Clone("AllUnc");

    TH1D *h_data;
    if( hName.find("_El") != std::string::npos && hName.find("_Mu") == std::string::npos ){
        h_data = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName).c_str()))->Clone("DATA");
    }else if( hName.find("_El") == std::string::npos && hName.find("_Mu") != std::string::npos ){
        h_data = (TH1D*)((TH1D*)f->Get(("DATA_Muon__"+hName).c_str()))->Clone("DATA");
    }else{
        h_data =  (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA");
        TH1D* h_data_mu = (TH1D*)f->Get(("DATA_Muon__"+hName+"_Mu").c_str());
        h_data.Add(h_data_mu);
        std::cout<<">> [INFO] Combined two channel for data"<<std::endl;
        delete h_data_mu;
    }

    fix(h_tt);    h_tt->Rebin(rebin);
    fix(h_ttbkg); h_ttbkg->Rebin(rebin);
    fix(h_t);     h_t->Rebin(rebin);
    fix(h_b);     h_b->Rebin(rebin);
    fix(h_bkg);   h_bkg->Rebin(rebin);
    fix(h_all);   h_all->Rebin(rebin);
    fix(h_data);  h_data->Rebin(rebin);

    int bins = h_all->GetXaxis()->GetLast();
    int xMin = h_all->GetXaxis()->GetBinLowEdge(1);
    int xMax = h_all->GetXaxis()->GetBinUpEdge(bins);

    TH1D* h_mcunc = (TH1D*)h_all->Clone("mc_unc");
    h_mcunc->Divide(h_all);
    for( int b=1; b<=bins; b++ )
    {  
        float error=h_mcunc->GetBinError(b); 
        h_mcunc->SetBinContent(b, 0.);
        h_mcunc->SetBinError(b, error);
    }

    if( unity )
    {
        float data_s = 1/(h_data->Integral());
        float mc_s   = 1/(h_all->Integral());
        h_data->Scale(data_s);
        h_tt->Scale(mc_s);
        h_ttbkg->Scale(mc_s);
        h_t->Scale(mc_s);
        h_b->Scale(mc_s);
        h_bkg->Scale(mc_s);
        h_all->Scale(mc_s);
    }    

    TH1D* h_ratio = (TH1D*)h_data->Clone("ratio"); 
    h_ratio->Add(h_all, -1);  
    h_ratio->Divide(h_all);  

    h_data->SetLineWidth(lineWidth);
    h_data->SetLineColor(1);
    h_data->SetMarkerColor(1);
    h_data->SetMarkerStyle(20);

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

    hs->GetXaxis()->SetTitle(xtitle.c_str());
    hs->GetYaxis()->SetTitle(ytitle.c_str());

    hs->GetXaxis()->SetLabelFont(42);
    hs->GetXaxis()->SetLabelSize(0.05);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetXaxis()->SetTitleOffset(0.98);
    hs->GetXaxis()->SetTitleFont(42);
    hs->GetYaxis()->SetLabelFont(42);
    hs->GetYaxis()->SetLabelSize(0.06);
    hs->GetYaxis()->SetTitleOffset(0.7);
    //hs->GetYaxis()->SetTitleSize(0.07);
    hs->GetYaxis()->SetTitleSize(0.079);
    hs->GetYaxis()->SetTitleFont(42);

    THStack* h_stack = new THStack("THStcak", "");
    h_stack->SetHistogram(hs);

    h_stack->Add(h_b);
    h_stack->Add(h_t);
    h_stack->Add(h_ttbkg);
    h_stack->Add(h_tt);

    if( logy ){
        if( h_all->GetMaximum() < h_data->GetMaximum()) h_stack->SetMaximum(h_data->GetMaximum()*10);
        if( h_all->GetMinimum() > 0 )
            h_stack->SetMinimum(h_all->GetMinimum()/10);
        else
            h_stack->SetMinimum(10);
    }else{
        if( h_all->GetMaximum() < h_data->GetMaximum() ) h_stack->SetMaximum(h_data->GetMaximum()+h_data->GetMaximum()/10);
    }

    TCanvas* c1;
    TPad *p1, *p2; 
    if( logy ) 
    {
        c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",59,67,1076,824);
        //c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",59,72,1083,826);
        p1 = new TPad( ("P1_Log_"+hName).c_str(), "", 0.01725746,0.3165195,0.9976679,0.9735183);
        p2 = new TPad( ("P2_Log_"+hName).c_str(), "", 0.01399254,0.03773585,0.9981343,0.318239);
    }
    else
    {
        c1 = new TCanvas( ("C_Linear_"+hName).c_str(), "",59,67,1076,824);
        //c1 = new TCanvas( ("C_Log_"+hName).c_str(), "",59,72,1083,826);
        p1 = new TPad( ("P1_Linear_"+hName).c_str(), "", 0.01725746,0.3165195,0.9976679,0.9735183);
        p2 = new TPad( ("P2_Linear_"+hName).c_str(), "", 0.01399254,0.03773585,0.9981343,0.318239);
    }
    c1->Range(-2.617021,-1.006188,2.382979,5.162451);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1234043);
    c1->SetRightMargin(0.07659575);
    c1->SetTopMargin(0.06209987);
    c1->SetBottomMargin(0.1376441);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    p1->Draw();
    p1->cd();
    //p1->Range(-5.804078,-0.0003963955,53.70895,0.07552421);
    p1->Range(-69.93206,-11.83204,525.1982,2254.328);
    p1->SetFillColor(0);
    p1->SetBorderMode(0);
    p1->SetBorderSize(2);
    //p1->SetRightMargin(0.06232159);
    //p1->SetTopMargin(0.0532708);
    //p1->SetBottomMargin(0.005221184);
    p1->SetLeftMargin(0.1175071);
    p1->SetRightMargin(0.04234062);
    p1->SetTopMargin(0.0532708);
    p1->SetBottomMargin(0.005221184);
    p1->SetFrameBorderMode(0);
    p1->SetFrameBorderMode(0);

    if( logy ) p1->SetLogy(1);
    else p1->SetLogy(0);

    TLegend *leg;
    //leg = new TLegend(0.6998097,0.6230769,0.9995243,0.9096154,NULL,"brNDC");
    leg = new TLegend(0.6979068,0.63162,0.9976213,0.9180829,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "Data", "lp");
    leg->AddEntry(h_tt, "t#bar{t}+jet (lepton+jet)", "f");
    leg->AddEntry(h_ttbkg, "t#bar{t}+jet (Other)", "f");
    leg->AddEntry(h_t, "Single top", "f");
    leg->AddEntry(h_b, "Z/#gamma*/W/WW/WZ/ZZ", "f");
    leg->AddEntry(h_bkg, "1#sigma non t#bar{t}+jet stat.", "f");
    leg->AddEntry(h_all, "1#sigma Total stat.", "f");

    TPaveText* t_title;
    t_title = new TPaveText(0.08705995,0.9123537,0.7169363,0.9963828,"brNDC");
    //t_title = new TPaveText(0.07088487,0.9153846,0.7007612,1,"brNDC");
    t_title->AddText("CMS, L = 19.7/fb, #sqrt{s} = 8TeV");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(11);
    t_title->SetTextSize(0.04805273);

    h_stack->Draw("HIST");
    h_all->Draw("SAMEE2");
    h_bkg->Draw("SAMEE2");
    h_data->Draw("ESAME");
    leg->Draw();
    t_title->Draw();

    p1->Modified();
    c1->cd();

    p2->Draw();
    p2->cd();
    //p2->Range(-5.977337,-1.471186,53.79603,1.930508);
    p2->Range(-71.06455,-0.4118644,526.3307,0.156511);
    p2->SetFillColor(0);
    p2->SetBorderMode(0);
    p2->SetBorderSize(0);
    p2->SetTickx(1);
    //p2->SetRightMargin(0.06350711);
    //p2->SetTopMargin(0.00896861);
    //p2->SetBottomMargin(0.4618834);
    p2->SetLeftMargin(0.1189573);
    p2->SetRightMargin(0.0440758);
    p2->SetTopMargin(0.01145553);
    p2->SetBottomMargin(0.4607244);
    p2->SetFrameBorderMode(0);
    p2->SetFrameBorderMode(0);

    float x0 = h_ratio->GetXaxis()->GetBinLowEdge(h_ratio->GetXaxis()->GetFirst());
    float x1 = h_ratio->GetXaxis()->GetBinUpEdge(h_ratio->GetXaxis()->GetLast()); 
    TLine* l = new TLine( x0, 0., x1, 0.);
    //TLine* l = new TLine( x0, 1., x1, 1.);
    l->SetLineColor(1);
    l->SetLineWidth(2);
 
    h_ratio->SetLineWidth(3);
    h_mcunc->SetMarkerStyle(0);
    h_mcunc->SetMarkerSize(0);
    h_mcunc->SetLineWidth(0);
    h_mcunc->SetFillStyle(1001);
    h_mcunc->SetFillColor(kYellow);
    h_mcunc->GetXaxis()->SetTitle(xtitle.c_str());
    h_mcunc->GetXaxis()->SetLabelFont(42);
    h_mcunc->GetXaxis()->SetLabelOffset(0.021);
    h_mcunc->GetXaxis()->SetLabelSize(0.17);
    h_mcunc->GetXaxis()->SetTitleSize(0.19);
    h_mcunc->GetXaxis()->SetTitleOffset(1.06);
    h_mcunc->GetXaxis()->SetTitleFont(42);
    h_mcunc->GetYaxis()->SetTitle("#frac{Data-MC}{MC}");
    //h_mcunc->GetYaxis()->SetTitle("(Data-MC)/MC");
    //h_mcunc->GetYaxis()->SetTitle("Data/MC");
    h_mcunc->GetYaxis()->SetNdivisions(505);
    h_mcunc->GetYaxis()->SetLabelFont(42);
    h_mcunc->GetYaxis()->SetLabelSize(0.14);
    h_mcunc->GetYaxis()->SetTitleSize(0.15);
    h_mcunc->GetYaxis()->SetTitleOffset(0.32);
    h_mcunc->GetYaxis()->SetTitleFont(42);
    h_mcunc->GetYaxis()->SetRangeUser(-0.15,0.15);
    //h_mcunc->GetYaxis()->SetRangeUser(-0.9,0.9);
    //h_mcunc->GetYaxis()->SetRangeUser(0.1,1.9);
    //h_mcunc->GetYaxis()->SetRangeUser(0.85,1.15);
    //h_mcunc->GetYaxis()->SetRangeUser(0.45,1.55);
    h_mcunc->Draw("E2");
    h_ratio->Draw("esame");
    l->Draw("same");
    p2->Modified();
    c1->cd();
    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);

    if( logy )
        c1->SaveAs((output+"/StackData_Log_"+hName+".pdf").c_str());
    else
        c1->SaveAs((output+"/StackData_Linear_"+hName+".pdf").c_str());
}
