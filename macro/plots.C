#include "TColor.h"
#include "help.C"
// For semi-lepton case
void drawObservableData( TFile* f, std::string output=".", std::string oName="O2", std::string eName="EvtChi2", std::string headName="0 b-jet CR", std::string xtitle="O_{2}/m_{top}^{3}", std::string ytitle="Events (k)", float scaleY=0.001, float xmin=-1, float xmax=1, int rebin=1 )
{
    //TCanvas *c1 = new TCanvas(("c1_obdist_"+eName+oName).c_str(), "c1",1394,85,881,747);
    TCanvas *c1 = new TCanvas(("c1_obdist_"+eName+oName).c_str(), "c1", 1394, 85, W, H);
    gStyle->SetOptStat(0);
    //c1->Range(-1.432161,-5.570236,1.234506,28.3845);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1620603);
    c1->SetRightMargin(0.0879397);
    c1->SetTopMargin(0.07329843);
    c1->SetBottomMargin(0.1640489);
    //c1->SetLeftMargin( L/W );
    //c1->SetRightMargin( R/W );
    //c1->SetTopMargin( T/H );
    //c1->SetBottomMargin( B/H );
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    TH1D* Oh_el = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+eName+"_"+oName+"_El").c_str()))->Clone();
    TH1D* Oh_mu = (TH1D*)((TH1D*)f->Get(("DATA_Muon__"+eName+"_"+oName+"_Mu").c_str()))->Clone();
    fix(Oh_el); Oh_el->Rebin(rebin);
    fix(Oh_mu); Oh_mu->Rebin(rebin);
    
    fixNewRange(Oh_el, xmin, xmax );
    fixNewRange(Oh_mu, xmin, xmax );

    Oh_el->Scale(scaleY);
    Oh_mu->Scale(scaleY);

    TH1D* Oh = (TH1D*)Oh_el->Clone();
    Oh->Add(Oh_mu);

    Oh->SetLineColor(4);
    Oh->SetLineWidth(3);
    Oh->GetXaxis()->SetTitle(xtitle.c_str());
    Oh->GetXaxis()->SetNdivisions(509);
    Oh->GetXaxis()->SetLabelSize(0.06);
    Oh->GetXaxis()->SetTitleSize(0.07);
    Oh->GetYaxis()->SetTitle(ytitle.c_str());
    Oh->GetYaxis()->SetNdivisions(505);
    Oh->GetYaxis()->SetLabelSize(0.06);
    Oh->GetYaxis()->SetTitleSize(0.07);
    Oh->GetYaxis()->SetTitleOffset(0.79);
    Oh->GetYaxis()->SetTitleFont(42);
    Oh->Draw("HISTE");

    //Oh_mu->SetLineColor(kOrange+4);
    //Oh_el->SetLineColor(kGreen+3);
    Oh_mu->SetLineColor(kRed);
    Oh_el->SetLineColor(8);
    Oh_mu->SetLineWidth(3);
    Oh_el->SetLineWidth(3);
    Oh_mu->Draw("SAMEHISTE");
    Oh_el->Draw("SAMEHISTE");

    TLegend *leg;
    leg = new TLegend(0.6319095,0.6684119,0.9208543,0.8813264,NULL,"brNDC");
    leg->SetTextSize(0.04861111);
    leg->SetHeader(headName.c_str());
    leg->SetBorderSize(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(Oh_mu,"Muon",    "le");
    leg->AddEntry(Oh_el,"Electron","le");
    leg->AddEntry(Oh,   "Combined","le");
    leg->Draw();

    TPaveText* t_title;
    t_title = new TPaveText(0.6595477,0.938918,0.9396985,0.9842932,"brNDC");
    t_title->AddText("19.7 fb^{-1} (8TeV)");
    t_title->SetFillColor(0);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextFont(42);
    t_title->SetTextSize(0.05235602);
    t_title->Draw();

    TPaveText* t_CMS;
    t_CMS = new TPaveText(0.2022613,0.8481675,0.2952261,0.8935428,"brNDC");
    t_CMS->AddText("CMS");
    t_CMS->SetFillColor(0);
    t_CMS->SetFillStyle(0);
    t_CMS->SetBorderSize(0);
    t_CMS->SetTextAlign(12);
    t_CMS->SetTextSize(0.06108203);
    t_CMS->Draw();

    TPaveText* t_perliminary;
    t_perliminary = new TPaveText(0.2022613,0.7853403,0.2952261,0.8307155,"brNDC");
    t_perliminary->AddText("Perliminary");
    t_perliminary->SetFillColor(0);
    t_perliminary->SetFillStyle(0);
    t_perliminary->SetBorderSize(0);
    t_perliminary->SetTextAlign(12);
    t_perliminary->SetTextFont(52);
    t_perliminary->SetTextSize(0.05235602);
    t_perliminary->Draw();

    c1->SaveAs((output+"/ObsDistribution_DATA_"+eName+oName+".pdf").c_str());
}

// For delphes case
void drawObservableDist( TFile* f, std::string output=".", std::string histName="Evt_O", std::string unit="O/M_{top}^{3}", std::string CHname="LepJets", bool legX=1 )
{
    bool isLepJets=false;
    bool isMultiJets=false;
    if( CHname.compare("LepJets") == 0 ) isLepJets=true;
    if( CHname.compare("MultiJets") == 0 ) isMultiJets=true;

    TCanvas *c1 = new TCanvas(("c1_obdist_"+histName).c_str(), "c1",57,97,1099,752);
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
    double wrtOh = 1/Oh->Integral();
    Oh->Scale(wrtOh);
    Oh->SetLineColor(4);
    Oh->SetLineWidth(3);
    Oh->GetXaxis()->SetTitle(unit.c_str());
    Oh->GetXaxis()->SetTitleSize(0.05);
    Oh->GetXaxis()->SetTitleFont(62);
    Oh->GetXaxis()->SetLabelSize(0.05);
    Oh->GetXaxis()->SetLabelFont(62);
    Oh->GetYaxis()->SetLabelFont(62);
    Oh->GetYaxis()->SetTitle("Events/Total");
    Oh->GetYaxis()->SetTitleSize(0.06);
    Oh->GetYaxis()->SetTitleOffset(0.77);
    Oh->GetYaxis()->SetTitleFont(42);
    Oh->GetZaxis()->SetLabelFont(62);
    Oh->GetZaxis()->SetLabelSize(0.035);
    Oh->GetZaxis()->SetTitleSize(0.035);
    Oh->GetZaxis()->SetTitleFont(62);
    Oh->Draw("HISTE");
    if( isLepJets ){
        Oh_mu = (TH1D*)((TH1D*)f->Get((histName+"_Mu").c_str()))->Clone("muon");
        Oh_el = (TH1D*)((TH1D*)f->Get((histName+"_El").c_str()))->Clone("electron");
        //double wrtmu = 1/Oh_mu->Integral();
        //double wrtel = 1/Oh_el->Integral();
        //Oh_mu->Scale(wrtmu);
        //Oh_el->Scale(wrtel);
        Oh_mu->Scale(wrtOh);
        Oh_el->Scale(wrtOh);
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
        //double wrtp5 = 1/Oh_p5->Integral();
        //double wrtp6 = 1/Oh_p6->Integral();
        //Oh_p5->Scale(wrtp5);
        //Oh_p6->Scale(wrtp6);
        Oh_p5->Scale(wrtOh);
        Oh_p6->Scale(wrtOh);
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

    TPaveText* t_title;
    t_title = new TPaveText(0.06940639,0.9337931,0.6995434,0.9793103,"brNDC");
    t_title->AddText("CMS Simulation, L = 19.7/fb, #sqrt{s} = 8TeV");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);
    t_title->Draw();

    c1->SaveAs((output+"/"+histName+".pdf").c_str());
}
void drawObservable( TFile* f, std::string histName, std::string output=".", std::string Oname="O", std::string CHname="", float range=0.1, bool legX=0 )
{
    TCanvas *c1 = new TCanvas(("c1_ob_"+histName).c_str(), "c1",113,93,1099,750);
    gStyle->SetOptStat(0);
    c1->Range(-0.25,1247.471,2.25,1397.296);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    TH1D* h0 = (TH1D*)f->Get(histName.c_str());
    TH1D* O = (TH1D*)h0->Clone("O");
    double wrt = 1/(O->Integral());
    O->Scale(wrt);	
    TH1D* OSigma = (TH1D*)O->Clone("CopyOriginal");
    OSigma->SetFillColor(kBlue-9);
    OSigma->GetXaxis()->SetLabelSize(0.08);
    OSigma->GetXaxis()->SetLabelFont(62);
    OSigma->GetXaxis()->SetTitleSize(0.035);
    OSigma->GetXaxis()->SetTitleFont(62);
    OSigma->GetYaxis()->SetLabelFont(62);
    OSigma->GetYaxis()->SetTitle("Events/Total");
    OSigma->GetYaxis()->SetTitleSize(0.06);
    OSigma->GetYaxis()->SetTitleOffset(0.77);
    OSigma->GetYaxis()->SetTitleFont(42);
    OSigma->GetYaxis()->SetRangeUser(0.5-range,0.5+range);
    OSigma->GetZaxis()->SetLabelFont(62);
    OSigma->GetZaxis()->SetLabelSize(0.035);
    OSigma->GetZaxis()->SetTitleSize(0.035);
    OSigma->GetZaxis()->SetTitleFont(62);
    OSigma->Draw("E2");

    O->SetLineColor(kBlue+3);
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

    TPaveText* t_title;
    t_title = new TPaveText(0.07214612,0.9004149,0.7022831,0.9460581,"brNDC");
    t_title->AddText("CMS Simulation, L = 19.7/fb, #sqrt{s} = 8TeV");
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);
    t_title->Draw();

    c1->SaveAs((output+"/"+histName+".pdf").c_str());
}
