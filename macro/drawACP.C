#include "caculate.C" 
void drawACP( std::string pathLepJets="../results/TTtoLepJet/result.root", 
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
	h->Fill(0, caculateACP( hlj_o2 ));
	h->SetBinError(1, caculateACPerror( hlj_o2 ));

	printf("O_{2}^{#mu}:\n");
	h->Fill(1,     caculateACP( hlj_o2_mu ));
	h->SetBinError(2, caculateACPerror( hlj_o2_mu ));

	printf("O_{2}^{e}:\n");
	h->Fill(2,   caculateACP( hlj_o2_el ));
	h->SetBinError(3, caculateACPerror( hlj_o2_el ));

	printf("O_{7}^{e+#mu}:\n");
	h->Fill(3, caculateACP( hlj_o7 ));
	h->SetBinError(4, caculateACPerror( hlj_o7 ));

	printf("O_{7}^{#mu}:\n");
	h->Fill(4,     caculateACP( hlj_o7_mu ));
	h->SetBinError(5, caculateACPerror( hlj_o7_mu ));

	printf("O_{7}^{e}:\n");
	h->Fill(5,   caculateACP( hlj_o7_el ));
	h->SetBinError(6, caculateACPerror( hlj_o7_el ));

	printf("O_{5}(p_{T}^{Bjet}#geq40):\n");
	h->Fill(6, caculateACP( hmj_o5_p4 ));
	h->SetBinError(7, caculateACPerror( hmj_o5_p4 ));

	printf("O_{5}(p_{T}^{Bjet}#geq50):\n");
	h->Fill(7, caculateACP( hmj_o5_p5 ));
	h->SetBinError(8, caculateACPerror( hmj_o5_p5 ));

	printf("O_{5}(p_{T}^{Bjet}#geq60):\n");
	h->Fill(8, caculateACP( hmj_o5_p6 ));
	h->SetBinError(9, caculateACPerror( hmj_o5_p6 ));

	printf("O_{7}(p_{T}^{Bjet}#geq40):\n");
	h->Fill(9, caculateACP( hmj_o7_p4 ));
	h->SetBinError(10, caculateACPerror( hmj_o7_p4 ));

	printf("O_{7}(p_{T}^{Bjet}#geq50):\n");
	h->Fill(10, caculateACP( hmj_o7_p5 ));
	h->SetBinError(11, caculateACPerror( hmj_o7_p5 ));

	printf("O_{7}(p_{T}^{Bjet}#geq60):\n");
	h->Fill(11, caculateACP( hmj_o7_p6 ));
	h->SetBinError(12, caculateACPerror( hmj_o7_p6 ));

	h->Scale(wrt);

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

	TH1D* h0 = (TH1D*)h->Clone("Original");
	TColor *color; // for color definition with alpha
	ci = TColor::GetColor("#9999ff");
	h->SetMaximum(0.06*wrt);
	h->SetMinimum(-0.06*wrt);
	h->SetFillColor(ci);
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
	
	c1->SaveAs((output+"/ACP_"+evtcat+".pdf").c_str());
}
