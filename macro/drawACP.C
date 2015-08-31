#include "caculate.C"
#include <fstream>
void drawACP( TFile* f,
            bool wrtError=true,
            std::string analysis="SemiLeptanic",
            std::string evtcat="Evt",
            std::string obs = "O2",
            std::string output=".",
            std::string ytitle="ACP",
            std::string xtitle="O_{2}",
            double wrt=1,
            int legX=1
            )
{
    bool murmur=true;
    TH1D *hObs, *hObs_mu, *hObs_el;
    if( analysis.compare("") != 0 )
    { 
        printf((analysis+"/"+evtcat+"_"+obs+"Asym").c_str());
        hObs    = (TH1D*)f->Get((analysis+"/"+evtcat+"_"+obs+"Asym").c_str()); 
        hObs_mu = (TH1D*)f->Get((analysis+"/"+evtcat+"_"+obs+"Asym_Mu").c_str()); 
        hObs_el = (TH1D*)f->Get((analysis+"/"+evtcat+"_"+obs+"Asym_El").c_str()); 
    }else{
        printf((analysis+"/"+evtcat+"_"+obs+"Asym").c_str());
        hObs    = (TH1D*)f->Get((evtcat+"_"+obs+"Asym").c_str()); 
        hObs_mu = (TH1D*)f->Get((evtcat+"_"+obs+"Asym_Mu").c_str()); 
        hObs_el = (TH1D*)f->Get((evtcat+"_"+obs+"Asym_El").c_str()); 
    }

    const int allh=3;

    fstream out;
    out.open((output+"/ACPCounting_"+evtcat+"_"+obs+".txt").c_str(),ios_base::out);
    out<<obs<<": "<<endl;
    out<<caculateACPDetail( hObs, wrtError )<<endl;
    out<<obs<<" muon: "<<endl;
    out<<caculateACPDetail( hObs_mu, wrtError )<<endl;
    out<<obs<<" electron: "<<endl;
    out<<caculateACPDetail( hObs_el, wrtError)<<endl;
    out.close();

    TH1D* h1s = new TH1D(("all"+evtcat+"_"+obs).c_str(), "", allh, 0, allh);

    h1s->GetXaxis()->SetBinLabel(1, (xtitle+"^{e+#mu}").c_str());
    h1s->GetXaxis()->SetBinLabel(2, (xtitle+"^{#mu}").c_str()  );
    h1s->GetXaxis()->SetBinLabel(3, (xtitle+"^{e}").c_str()    );

    double eObs(0), eObsMu(0), eObsEl(0); 
    printf("%s Mean Combined, Muon, Electron channel:\n", obs.c_str());
        h1s->Fill( 0., caculateACP( hObs,    murmur ));
        h1s->Fill( 1., caculateACP( hObs_mu, murmur ));
        h1s->Fill( 2., caculateACP( hObs_el, murmur ));

    printf("%s Error Combined, Muon, Electron channel:\n", obs.c_str());
        eObs   = (wrtError)? caculateACPerrorWrt( hObs    ):caculateACPerror( hObs    );
        eObsMu = (wrtError)? caculateACPerrorWrt( hObs_mu ):caculateACPerror( hObs_mu );
        eObsEl = (wrtError)? caculateACPerrorWrt( hObs_el ):caculateACPerror( hObs_el );

    h1s->SetBinError( 1, eObs   );
    h1s->SetBinError( 2, eObsMu );
    h1s->SetBinError( 3, eObsEl );

    TH1D* h2s = (TH1D*)h1s->Clone("h_2sigma");
    h2s->SetBinError( 1, 2*eObs   );
    h2s->SetBinError( 2, 2*eObsMu );
    h2s->SetBinError( 3, 2*eObsEl );

    TCanvas *c1 = new TCanvas(("CAcp_"+evtcat+obs).c_str(), "c1", 41,111,859,637); c1->Clear();
    gStyle->SetOptStat(0);
    c1->Range(-0.4717744,-0.1330709,3.184476,0.1169291);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1290323);
    c1->SetRightMargin(0.05045492);
    c1->SetTopMargin(0.06771654);
    c1->SetBottomMargin(0.1322835);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    h2s->SetMaximum(0.1*wrt);
    h2s->SetMinimum(-0.1*wrt);
    //h->SetFillColor(kBlue-9);
    h2s->SetFillColor(kYellow);
    h2s->GetXaxis()->SetLabelOffset(0.01);
    h2s->GetXaxis()->SetLabelSize(0.07);
    h2s->GetXaxis()->SetLabelFont(62);
    h2s->GetXaxis()->SetTitleSize(0.035);
    h2s->GetYaxis()->SetTitle(ytitle.c_str());
    h2s->GetYaxis()->SetLabelOffset(0.01);
    h2s->GetYaxis()->SetLabelSize(0.05);
    h2s->GetYaxis()->SetTitleSize(0.06);
    h2s->GetYaxis()->SetTitleOffset(1.05);
    h2s->GetYaxis()->SetTitleFont(42);
    h2s->GetZaxis()->SetLabelSize(0.035);
    h2s->GetZaxis()->SetTitleSize(0.035);
    h2s->Draw("E2");

    h1s->SetFillColor(kGreen);
    h1s->Draw("E2SAME");

    TH1D* h0 = (TH1D*)h1s->Clone("Original");
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
    leg->AddEntry(h1s,"1#sigma stat. error","f");
    leg->AddEntry(h2s,"2#sigma stat. error","f");
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

    c1->SaveAs((output+"/ACP_"+evtcat+"_"+obs+".pdf").c_str());
}
void drawACPLepJet( TFile* f, 
        bool wrtError=true,
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
        hlj_o2    = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym").c_str()); 
        hlj_o2_mu = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym_Mu").c_str()); 
        hlj_o2_el = (TH1D*)f->Get((analysis+"/"+evtcat+"_O2Asym_El").c_str()); 
        hlj_o7    = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym").c_str()); 
        hlj_o7_mu = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym_Mu").c_str()); 
        hlj_o7_el = (TH1D*)f->Get((analysis+"/"+evtcat+"_O7Asym_El").c_str()); 
    }else{
        hlj_o2    = (TH1D*)f->Get((evtcat+"_O2Asym").c_str()); 
        hlj_o2_mu = (TH1D*)f->Get((evtcat+"_O2Asym_Mu").c_str()); 
        hlj_o2_el = (TH1D*)f->Get((evtcat+"_O2Asym_El").c_str()); 
        hlj_o7    = (TH1D*)f->Get((evtcat+"_O7Asym").c_str()); 
        hlj_o7_mu = (TH1D*)f->Get((evtcat+"_O7Asym_Mu").c_str()); 
        hlj_o7_el = (TH1D*)f->Get((evtcat+"_O7Asym_El").c_str()); 
    }

    const int allh=6;

    fstream out;
    out.open((output+"/ACPCounting_"+evtcat+".txt").c_str(),ios_base::out);
    out<<"O2: "<<endl;
    out<<caculateACPDetail( hlj_o2, wrtError )<<endl;
    out<<"O2 muon: "<<endl;
    out<<caculateACPDetail( hlj_o2_mu, wrtError )<<endl;
    out<<"O2 electron: "<<endl;
    out<<caculateACPDetail( hlj_o2_el, wrtError)<<endl;
    out<<"O7: "<<endl;
    out<<caculateACPDetail( hlj_o7, wrtError)<<endl;
    out<<"O7 muon: "<<endl;
    out<<caculateACPDetail( hlj_o7_mu, wrtError )<<endl;
    out<<"O7 electron: "<<endl;
    out<<caculateACPDetail( hlj_o7_el, wrtError)<<endl;
    out.close();

    TH1D* h = new TH1D(("all"+evtcat).c_str(), "", allh, 0, allh);

    h->GetXaxis()->SetBinLabel(1, "O_{2}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(2, "O_{2}^{#mu}");
    h->GetXaxis()->SetBinLabel(3, "O_{2}^{e}");
    h->GetXaxis()->SetBinLabel(4, "O_{7}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(5, "O_{7}^{#mu}");
    h->GetXaxis()->SetBinLabel(6, "O_{7}^{e}");

    printf("O_{2}^{e+#mu}:\n");
    h->Fill(0., caculateACP( hlj_o2, murmur ));
    if( wrtError ) h->SetBinError(1, caculateACPerrorWrt( hlj_o2 ));
    else h->SetBinError(1, caculateACPerror( hlj_o2 ));

    printf("O_{2}^{#mu}:\n");
    h->Fill(1.,     caculateACP( hlj_o2_mu, murmur ));
    if( wrtError ) h->SetBinError(2, caculateACPerrorWrt( hlj_o2_mu ));
    else h->SetBinError(2, caculateACPerror( hlj_o2_mu ));

    printf("O_{2}^{e}:\n");
    h->Fill(2.,   caculateACP( hlj_o2_el, murmur ));
    if( wrtError ) h->SetBinError(3, caculateACPerrorWrt( hlj_o2_el ));
    else h->SetBinError(3, caculateACPerror( hlj_o2_el ));

    printf("O_{7}^{e+#mu}:\n");
    h->Fill(3., caculateACP( hlj_o7, murmur ));
    if( wrtError ) h->SetBinError(4, caculateACPerrorWrt( hlj_o7 ));
    else h->SetBinError(4, caculateACPerror( hlj_o7 ));

    printf("O_{7}^{#mu}:\n");
    h->Fill(4.,     caculateACP( hlj_o7_mu, murmur ));
    if( wrtError ) h->SetBinError(5, caculateACPerrorWrt( hlj_o7_mu ));
    else h->SetBinError(5, caculateACPerror( hlj_o7_mu ));

    printf("O_{7}^{e}:\n");
    h->Fill(5.,   caculateACP( hlj_o7_el, murmur ));
    if( wrtError ) h->SetBinError(6, caculateACPerrorWrt( hlj_o7_el ));
    else h->SetBinError(6, caculateACPerror( hlj_o7_el ));

    TCanvas *c1 = new TCanvas(("CAcp_"+evtcat).c_str(), "c1",41,89,1213,664);
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

    c1->SaveAs((output+"/ACP_"+evtcat+".pdf").c_str());
}
void drawACP2Channel( std::string pathLepJets="../results/TTtoLepJet/result.root", 
        std::string pathMultiJets="../results/TTtoMultiJets/result.root", 
        std::string evtcat="Evt2b",
        std::string output=".",
        std::string ytitle="ACP",
        double wrt=1,
        int legX=1 )
{
    TFile* flj = new TFile(pathLepJets.c_str());
    TFile* fmj = new TFile(pathMultiJets.c_str());

    TH1D*  hlj_o2    = (TH1D*)flj->Get((evtcat+"_O2Asym").c_str()); 
    TH1D*  hlj_o2_mu = (TH1D*)flj->Get((evtcat+"_O2Asym_Mu").c_str()); 
    TH1D*  hlj_o2_el = (TH1D*)flj->Get((evtcat+"_O2Asym_El").c_str()); 
    TH1D*  hlj_o7    = (TH1D*)flj->Get((evtcat+"_O7Asym").c_str()); 
    TH1D*  hlj_o7_mu = (TH1D*)flj->Get((evtcat+"_O7Asym_Mu").c_str()); 
    TH1D*  hlj_o7_el = (TH1D*)flj->Get((evtcat+"_O7Asym_El").c_str()); 

    TH1D*  hmj_o5_p4 = (TH1D*)fmj->Get((evtcat+"_O5Asym").c_str()); 
    TH1D*  hmj_o5_p5 = (TH1D*)fmj->Get((evtcat+"_O5Asym_PT50").c_str()); 
    TH1D*  hmj_o5_p6 = (TH1D*)fmj->Get((evtcat+"_O5Asym_PT60").c_str()); 
    TH1D*  hmj_o7_p4 = (TH1D*)fmj->Get((evtcat+"_O7Asym").c_str()); 
    TH1D*  hmj_o7_p5 = (TH1D*)fmj->Get((evtcat+"_O7Asym_PT50").c_str()); 
    TH1D*  hmj_o7_p6 = (TH1D*)fmj->Get((evtcat+"_O7Asym_PT60").c_str());

    const int allh=12;

    TH1D* h = new TH1D("all", "", allh, 0, allh);

    h->GetXaxis()->SetBinLabel(1, "O_{2}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(2, "O_{2}^{#mu}");
    h->GetXaxis()->SetBinLabel(3, "O_{2}^{e}");
    h->GetXaxis()->SetBinLabel(4, "O_{7}^{e+#mu}");
    h->GetXaxis()->SetBinLabel(5, "O_{7}^{#mu}");
    h->GetXaxis()->SetBinLabel(6, "O_{7}^{e}");
    h->GetXaxis()->SetBinLabel(7, "O_{5}(p_{T}^{Bjet}#geq40)");
    h->GetXaxis()->SetBinLabel(8, "O_{5}(p_{T}^{Bjet}#geq50)");
    h->GetXaxis()->SetBinLabel(9, "O_{5}(p_{T}^{Bjet}#geq60)");
    h->GetXaxis()->SetBinLabel(10, "O_{7}(p_{T}^{Bjet}#geq40)");
    h->GetXaxis()->SetBinLabel(11, "O_{7}(p_{T}^{Bjet}#geq50)");
    h->GetXaxis()->SetBinLabel(12, "O_{7}(p_{T}^{Bjet}#geq60)");

    printf("O_{2}^{e+#mu}:\n");
    h->Fill(0., caculateACP( hlj_o2 ));
    h->SetBinError(1, caculateACPerrorWrt( hlj_o2 ));
    //h->SetBinError(1, caculateACPerror( hlj_o2 ));

    printf("O_{2}^{#mu}:\n");
    h->Fill(1.,     caculateACP( hlj_o2_mu ));
    h->SetBinError(2, caculateACPerrorWrt( hlj_o2_mu ));
    //h->SetBinError(2, caculateACPerror( hlj_o2_mu ));

    printf("O_{2}^{e}:\n");
    h->Fill(2.,   caculateACP( hlj_o2_el ));
    h->SetBinError(3, caculateACPerrorWrt( hlj_o2_el ));
    //h->SetBinError(3, caculateACPerror( hlj_o2_el ));

    printf("O_{7}^{e+#mu}:\n");
    h->Fill(3., caculateACP( hlj_o7 ));
    h->SetBinError(4, caculateACPerrorWrt( hlj_o7 ));
    //h->SetBinError(4, caculateACPerror( hlj_o7 ));

    printf("O_{7}^{#mu}:\n");
    h->Fill(4.,     caculateACP( hlj_o7_mu ));
    h->SetBinError(5, caculateACPerrorWrt( hlj_o7_mu ));
    //h->SetBinError(5, caculateACPerror( hlj_o7_mu ));

    printf("O_{7}^{e}:\n");
    h->Fill(5.,   caculateACP( hlj_o7_el ));
    h->SetBinError(6, caculateACPerrorWrt( hlj_o7_el ));
    //h->SetBinError(6, caculateACPerror( hlj_o7_el ));

    printf("O_{5}(p_{T}^{Bjet}#geq40):\n");
    h->Fill(6., caculateACP( hmj_o5_p4 ));
    h->SetBinError(7, caculateACPerrorWrt( hmj_o5_p4 ));
    //h->SetBinError(7, caculateACPerror( hmj_o5_p4 ));

    printf("O_{5}(p_{T}^{Bjet}#geq50):\n");
    h->Fill(7., caculateACP( hmj_o5_p5 ));
    h->SetBinError(8, caculateACPerrorWrt( hmj_o5_p5 ));
    //h->SetBinError(8, caculateACPerror( hmj_o5_p5 ));

    printf("O_{5}(p_{T}^{Bjet}#geq60):\n");
    h->Fill(8., caculateACP( hmj_o5_p6 ));
    h->SetBinError(9, caculateACPerrorWrt( hmj_o5_p6 ));
    //h->SetBinError(9, caculateACPerror( hmj_o5_p6 ));

    printf("O_{7}(p_{T}^{Bjet}#geq40):\n");
    h->Fill(9., caculateACP( hmj_o7_p4 ));
    h->SetBinError(10, caculateACPerrorWrt( hmj_o7_p4 ));
    //h->SetBinError(10, caculateACPerror( hmj_o7_p4 ));

    printf("O_{7}(p_{T}^{Bjet}#geq50):\n");
    h->Fill(10., caculateACP( hmj_o7_p5 ));
    h->SetBinError(11, caculateACPerrorWrt( hmj_o7_p5 ));
    //h->SetBinError(11, caculateACPerror( hmj_o7_p5 ));

    printf("O_{7}(p_{T}^{Bjet}#geq60):\n");
    h->Fill(11., caculateACP( hmj_o7_p6 ));
    h->SetBinError(12, caculateACPerrorWrt( hmj_o7_p6 ));
    //h->SetBinError(12, caculateACPerror( hmj_o7_p6 ));

    //h->Scale(wrt);

    TCanvas *c1 = new TCanvas("c1", "c1",41,89,1213,664);
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

    TLine* line = new TLine(0,0,allh,0);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();

    TH1D* h0 = (TH1D*)h->Clone("Original");
    h->SetMaximum(0.06*wrt);
    h->SetMinimum(-0.06*wrt);
    h->SetFillColor(kBlue-9);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetLabelSize(0.05);
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
    //TLine* line = new TLine(0,0,allh,0);
    //line->SetLineColor(2);
    //line->SetLineWidth(3);
    //line->Draw();

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

    c1->SaveAs((output+"/ACP_"+evtcat+".pdf").c_str());
}
