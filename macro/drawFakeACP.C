#include "caculate.C"
#include <fstream>
void drawACPLepJet( TFile* f, 
        std::string evtcat="Evt",
        std::string analysis="SemiLeptanic",
        std::string output=".",
        std::string ytitle="ACP",
        double wrt=1,
        int legX=1 )
{
    bool murmur=true;
    TH1D *hlj_o2, *hlj_o2_mu, *hlj_o2_el, *hlj_o7, *hlj_o7_mu, *hlj_o7_el;
    if( analysis.compare("") != 0 )
    { 
        hlj_o2    = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym_LepZ").c_str()); 
        hlj_o2_mu = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym_LepZ_Mu").c_str()); 
        hlj_o2_el = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym_LepZ_El").c_str()); 
        hlj_o7    = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym_JetZ").c_str()); 
        hlj_o7_mu = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym_JetZ_Mu").c_str()); 
        hlj_o7_el = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym_JetZ_El").c_str()); 
    }else{
        hlj_o2    = (TH1D*)f->Get((evtcat+"_O2Asym_LepZ").c_str()); 
        hlj_o2_mu = (TH1D*)f->Get((evtcat+"_O2Asym_LepZ_Mu").c_str()); 
        hlj_o2_el = (TH1D*)f->Get((evtcat+"_O2Asym_LepZ_El").c_str()); 
        hlj_o7    = (TH1D*)f->Get((evtcat+"_O2Asym_JetZ").c_str()); 
        hlj_o7_mu = (TH1D*)f->Get((evtcat+"_O2Asym_JetZ_Mu").c_str()); 
        hlj_o7_el = (TH1D*)f->Get((evtcat+"_O2Asym_JetZ_El").c_str()); 
    }

    const int allh=6;

    fstream out;
    out.open((output+"/FakeACPCounting_"+evtcat+".txt").c_str(),ios_base::out);
    out<<"O2 LepZ: "<<endl;
    out<<caculateACPDetail( hlj_o2 )<<endl;
    out<<"O2 LepZ muon: "<<endl;
    out<<caculateACPDetail( hlj_o2_mu )<<endl;
    out<<"O2 LepZ electron: "<<endl;
    out<<caculateACPDetail( hlj_o2_el)<<endl;
    out<<"O2 JetZ: "<<endl;
    out<<caculateACPDetail( hlj_o7)<<endl;
    out<<"O2 JetZ muon: "<<endl;
    out<<caculateACPDetail( hlj_o7_mu )<<endl;
    out<<"O2 JetZ electron: "<<endl;
    out<<caculateACPDetail( hlj_o7_el)<<endl;
    out.close();

    TH1D* h = new TH1D(("all"+evtcat).c_str(), "", allh, 0, allh);

    h->GetXaxis()->SetBinLabel(1, "O_{2,lepZ}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(2, "O_{2,lepZ}^{#mu}");
    h->GetXaxis()->SetBinLabel(3, "O_{2,lepZ}^{e}");
    h->GetXaxis()->SetBinLabel(4, "O_{2,JetZ}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(5, "O_{2,JetZ}^{#mu}");
    h->GetXaxis()->SetBinLabel(6, "O_{2,JetZ}^{e}");

    printf("O_{2,lepZ}^{e+#mu}:\n");
    h->Fill(0., caculateACP( hlj_o2, murmur ));
    //h->SetBinError(1, caculateACPerrorWrt( hlj_o2 ));
    h->SetBinError(1, caculateACPerror( hlj_o2 ));

    printf("O_{2,lepZ}^{#mu}:\n");
    h->Fill(1.,     caculateACP( hlj_o2_mu, murmur ));
    //h->SetBinError(2, caculateACPerrorWrt( hlj_o2_mu ));
    h->SetBinError(2, caculateACPerror( hlj_o2_mu ));

    printf("O_{2,lepZ}^{e}:\n");
    h->Fill(2.,   caculateACP( hlj_o2_el, murmur ));
    //h->SetBinError(3, caculateACPerrorWrt( hlj_o2_el ));
    h->SetBinError(3, caculateACPerror( hlj_o2_el ));

    printf("O_{2,JetZ}^{e+#mu}:\n");
    h->Fill(3., caculateACP( hlj_o7, murmur ));
    //h->SetBinError(4, caculateACPerrorWrt( hlj_o7 ));
    h->SetBinError(4, caculateACPerror( hlj_o7 ));

    printf("O_{2,JetZ}^{#mu}:\n");
    h->Fill(4.,     caculateACP( hlj_o7_mu, murmur ));
    //h->SetBinError(5, caculateACPerrorWrt( hlj_o7_mu ));
    h->SetBinError(5, caculateACPerror( hlj_o7_mu ));

    printf("O_{2,JetZ}^{e}:\n");
    h->Fill(5.,   caculateACP( hlj_o7_el, murmur ));
    //h->SetBinError(6, caculateACPerrorWrt( hlj_o7_el ));
    h->SetBinError(6, caculateACPerror( hlj_o7_el ));

    //h->Scale(wrt);

    TCanvas *c1 = new TCanvas(("CfakeAcp_"+evtcat).c_str(), "c1",41,89,1213,664);
    gStyle->SetOptStat(0);
    c1->Range(-1.887097,-0.07984252,12.7379,0.07015748);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1290323);
    c1->SetRightMargin(0.05045492);
    c1->SetTopMargin(0.06771654);
    c1->SetBottomMargin(0.1322835);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    TH1D* h0 = (TH1D*)h->Clone("Original");
    h->SetMaximum(0.1*wrt);
    h->SetMinimum(-0.1*wrt);
    h->SetFillColor(kBlue-9);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetLabelSize(0.07);
    h->GetXaxis()->SetLabelFont(62);
    h->GetXaxis()->SetTitleSize(0.035);
    h->GetYaxis()->SetTitle(ytitle.c_str());
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.83);
    h->GetYaxis()->SetTitleFont(42);
    h->GetZaxis()->SetLabelSize(0.035);
    h->GetZaxis()->SetTitleSize(0.035);
    h->Draw("E2");

    h0->SetMarkerStyle(21);
    h0->SetMarkerSize(2);
    h0->Draw("psame");
    TLine* line = new TLine(0,0,allh,0);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();

    TLegend *leg;
    if( legX == 0 ) //Left 
        leg = new TLegend(0.173516,0.6726768,0.4310502,0.8363384,NULL,"brNDC");
    else //right
        leg = new TLegend(0.7237386,0.7739403,0.9652605,0.8854003,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h,"1#sigma stat. error","f");
    leg->Draw();

    TPaveText* t_title;
    t_title = new TPaveText(0.09842845,0.9387755,0.7278743,0.9843014,"brNDC");
    t_title->AddText("CMS Simulation, L = 19.7/fb, #sqrt{s} = 8TeV");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);
    t_title->Draw();

    c1->SaveAs((output+"/FakeACP_"+evtcat+".pdf").c_str());
}

