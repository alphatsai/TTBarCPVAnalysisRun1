void drawObservableDist( TFile* f, std::string output=".", std::string histName="Evt_O", std::string unit="O/M_{top}^{3}", std::string CHname="LepJets", bool legX=1 )
{
   bool isLepJets=false;
   bool isMultiJets=false;
   if( CHname.compare("LepJets") == 0 ) isLepJets=true;
   if( CHname.compare("MultiJets") == 0 ) isMultiJets=true;

   TCanvas *c1 = new TCanvas("c1", "c1",57,97,1099,752);
   gStyle->SetOptStat(0);
   c1->Range(-2.530337,-153.6535,2.391011,1015.73);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetRightMargin(0.07945205);
   c1->SetTopMargin(0.0691563);
   c1->SetBottomMargin(0.131397);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
	
   TH1D* Oh0 = (TH1D*)f->Get(histName.c_str());
   TH1D* Oh = (TH1D*)Oh0->Clone("CopyOriginal");
   TH1D* Oh_mu, *Oh_el, *Oh_p5, *Oh_p6;
   Oh->SetLineColor(4);
   Oh->SetLineWidth(3);
   Oh->GetXaxis()->SetTitle(unit.c_str());
   Oh->GetXaxis()->SetTitleSize(0.05);
   Oh->GetXaxis()->SetTitleFont(62);
   Oh->GetXaxis()->SetLabelSize(0.05);
   Oh->GetXaxis()->SetLabelFont(62);
   Oh->GetYaxis()->SetLabelFont(62);
   Oh->GetYaxis()->SetTitle("Events");
   Oh->GetYaxis()->SetTitleSize(0.06);
   Oh->GetYaxis()->SetTitleOffset(0.77);
   Oh->GetYaxis()->SetTitleFont(42);
   Oh->GetZaxis()->SetLabelFont(62);
   Oh->GetZaxis()->SetLabelSize(0.035);
   Oh->GetZaxis()->SetTitleSize(0.035);
   Oh->GetZaxis()->SetTitleFont(62);
   Oh->Draw("HISTE");
   if( isLepJets ){
    Oh_mu = (TH1D*)f->Get((histName+"_Mu").c_str());
    Oh_el = (TH1D*)f->Get((histName+"_El").c_str());
    Oh_mu->SetLineColor(kOrange+4);
    Oh_el->SetLineColor(kGreen+3);
    Oh_mu->SetLineWidth(3);
    Oh_el->SetLineWidth(3);
    Oh_mu->Draw("SAMEHISTE");
    Oh_el->Draw("SAMEHISTE");
   }	
   if( isMultiJets ){
    Oh_p5 = (TH1D*)f->Get((histName+"_PT50").c_str());
    Oh_p6 = (TH1D*)f->Get((histName+"_PT60").c_str());
    Oh_p5->SetLineColor(kOrange+4);
    Oh_p6->SetLineColor(kGreen+3);
    Oh_p5->SetLineWidth(3);
    Oh_p6->SetLineWidth(3);
    Oh_p5->Draw("SAMEHISTE");
    Oh_p6->Draw("SAMEHISTE");
   }
   TLegend *leg;
   if( legX == 0 ) //Left 
   	//leg = new TLegend(0.173516,0.6726768,0.4310502,0.8363384,NULL,"brNDC");
   	leg = new TLegend(0.153516,0.6726768,0.4110502,0.8363384,NULL,"brNDC");
   else //right
    leg = new TLegend(0.6283105,0.7026279,0.9324201,0.8990318,NULL,"brNDC");
	
   leg->SetBorderSize(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   if( isMultiJets ){
    leg->AddEntry(Oh,"Multi-Jets, p_{T}(Bjet)>40 GeV","l");
    leg->AddEntry(Oh_p5,"Multi-Jets, p_{T}(Bjet)>50 GeV","l");
    leg->AddEntry(Oh_p6,"Multi-Jets, p_{T}(Bjet)>60 GeV","l");
   }
   if( isLepJets ){
    leg->AddEntry(Oh,"Combined channel","l");
    leg->AddEntry(Oh_mu,"Muon+Jets channel","l");
    leg->AddEntry(Oh_el,"Electron+Jets channel","l");
   }
   leg->Draw();
   c1->SaveAs((output+"/"+histName+".pdf").c_str());
}

void drawObservable( TFile* f, std::string output=".", std::string histName, std::string Oname="O", std::string CHname="", bool legX=0 )
{
   TCanvas *c1 = new TCanvas("c1", "c1",113,93,1099,750);
   gStyle->SetOptStat(0);
   c1->Range(-0.25,1247.471,2.25,1397.296);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
	
   TH1D* OSigma0 = (TH1D*)f->Get(histName.c_str());
   TH1D* OSigma = (TH1D*)OSigma0->Clone("CopyOriginal");
   TH1D *O = (TH1D*)OSigma->Clone("O");
   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9999ff");
   OSigma->SetFillColor(ci);
   OSigma->GetXaxis()->SetLabelSize(0.08);
   OSigma->GetXaxis()->SetLabelFont(62);
   OSigma->GetXaxis()->SetTitleSize(0.035);
   OSigma->GetXaxis()->SetTitleFont(62);
   OSigma->GetYaxis()->SetLabelFont(62);
   OSigma->GetYaxis()->SetTitle("Events");
   OSigma->GetYaxis()->SetTitleSize(0.06);
   OSigma->GetYaxis()->SetTitleOffset(0.77);
   OSigma->GetYaxis()->SetTitleFont(42);
   OSigma->GetZaxis()->SetLabelFont(62);
   OSigma->GetZaxis()->SetLabelSize(0.035);
   OSigma->GetZaxis()->SetTitleSize(0.035);
   OSigma->GetZaxis()->SetTitleFont(62);
   OSigma->Draw("E2");
 
   ci = TColor::GetColor("#000099");
   O->SetLineColor(ci);
   O->SetLineWidth(3);
   O->GetXaxis()->SetBinLabel(1,(Oname+"<0").c_str());
   O->GetXaxis()->SetBinLabel(2,(Oname+">0").c_str());
   O->GetXaxis()->SetLabelFont(42);
   O->GetXaxis()->SetLabelSize(0.035);
   O->GetXaxis()->SetTitleSize(0.035);
   O->GetXaxis()->SetTitleFont(42);
   O->GetYaxis()->SetLabelFont(42);
   O->GetYaxis()->SetLabelSize(0.035);
   O->GetYaxis()->SetTitleSize(0.035);
   O->GetYaxis()->SetTitleFont(42);
   O->GetZaxis()->SetLabelFont(42);
   O->GetZaxis()->SetLabelSize(0.035);
   O->GetZaxis()->SetTitleSize(0.035);
   O->GetZaxis()->SetTitleFont(42);
   O->Draw("histsametext0");

   TLegend *leg;
   if( legX == 0 ) //Left 
   	leg = new TLegend(0.173516,0.6726768,0.4310502,0.8363384,NULL,"brNDC");
   else //right
    leg = new TLegend(0.5780822,0.7040111,0.8356164,0.8672199,NULL,"brNDC");
	
   leg->SetBorderSize(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(OSigma,"1#sigma stat.","f");
   if( CHname.size() == 0 ){
   	leg->AddEntry(O,Oname.c_str(),"l");
   }else{
	std::string legs=Oname+", "+CHname+" channel";
   	leg->AddEntry(O,legs.c_str(),"l");
   }
   leg->Draw();
   c1->SaveAs((output+"/"+histName+".pdf").c_str());
}
