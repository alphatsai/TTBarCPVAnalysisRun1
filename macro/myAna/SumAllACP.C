#include <iostream>
#include <vector>
#include <string>
#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "/Users/Alpha/myAna/TTBarCPV/TTBarCPVAnalysisRun1/macro/CMS_lumi.C"
void copyToSumHist( TH1D* sumHist, TH1D* h, int onum, std::string label )
{ 
    int osizes=sumHist->GetNbinsX()/3; 
    int bin=onum+1;
    sumHist->GetXaxis()->SetBinLabel( bin,          label.c_str());
    sumHist->GetXaxis()->SetBinLabel( bin+osizes,   label.c_str());
    sumHist->GetXaxis()->SetBinLabel( bin+2*osizes, label.c_str());
    // Electron channel
    sumHist->SetBinContent( bin,         h->GetBinContent(2) );
    sumHist->SetBinError(   bin,         h->GetBinError(2)   );
    // Muon channel
    sumHist->SetBinContent( bin+osizes,  h->GetBinContent(3) );
    sumHist->SetBinError(   bin+osizes,  h->GetBinError(3)   );
    // Combined channel
    sumHist->SetBinContent( bin+2*osizes, h->GetBinContent(1) );
    sumHist->SetBinError(   bin+2*osizes, h->GetBinError(1)   );
}
void SumAllACP( std::vector<std::string> Onames, std::vector<std::string> Olabels, std::string inputpath=".", std::string output="." )
{
    if( Onames.size() != Olabels.size() ) return;

    const int Osize=Onames.size();
    const int Bins=Osize*3;
    TH1D* hall = new TH1D("All", "", Bins, 0, Bins ); hall->Sumw2();
    TH1D* hbkg = new TH1D("Bkg", "", Bins, 0, Bins ); hbkg->Sumw2();
    TH1D* hsig = new TH1D("Sig", "", Bins, 0, Bins ); hsig->Sumw2();
    //TGraphAsymmErrors* hsigE;
 
    for( int i=0; i<Osize; i++)
    //for( int i=0; i<2; i++)
    {
        TFile* fin = new TFile((inputpath+"/FinalACPSumHist_"+Onames[i]+".root").c_str());
        TH1D* h_ACP = (TH1D*)fin->Get("ACP");
        TH1D* h_Bkg = (TH1D*)fin->Get("ACPbkg");
        TH1D* h_Sig = (TH1D*)fin->Get("ACPsig");
        copyToSumHist( hall, h_ACP, i, Olabels[i] );
        copyToSumHist( hbkg, h_Bkg, i, Olabels[i] );
        copyToSumHist( hsig, h_Sig, i, Olabels[i] );
        fin->Close();
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

    TH1D* hbkgLine = (TH1D*)hbkg->Clone();
    float width=4;
    //float width=3;
    hbkg->SetMaximum(    width );
    hbkg->SetMinimum( -1*width );
    hbkg->GetXaxis()->SetLabelOffset(0.01);
    hbkg->GetXaxis()->SetLabelSize(0.07);
    hbkg->GetXaxis()->SetLabelFont(62);
    hbkg->GetXaxis()->SetTitleSize(0.035);
    hbkg->GetYaxis()->SetTitle("A'_{CP} [%]");
    hbkg->GetYaxis()->SetLabelOffset(0.01);
    hbkg->GetYaxis()->SetLabelSize(0.06);
    hbkg->GetYaxis()->SetTitleSize(0.07);
    hbkg->GetYaxis()->SetTitleOffset(0.84);
    hbkg->GetYaxis()->SetTitleFont(42);
    hbkg->GetZaxis()->SetLabelSize(0.035);
    hbkg->GetZaxis()->SetTitleSize(0.035);
    hbkg->SetLineWidth(2);
    hbkg->SetLineColor(kGreen+2);
    hbkg->SetFillColor(kGreen-8);
    hbkg->SetFillStyle(3008);
    hbkg->Draw("e2");

    hbkgLine->SetLineWidth(2);
    hbkgLine->SetLineColor(kGreen+2);
    hbkgLine->Draw("histsame");

    hall->SetMarkerStyle(4);
    hall->SetMarkerSize(1.5);
    hall->SetMarkerColor(2);
    hall->SetLineColor(2);
    hall->SetLineWidth(2);
    hall->Draw("pe1same");

    hsig->SetLineColor(1);
    hsig->SetLineWidth(2);
    hsig->SetMarkerColor(1);
    hsig->SetMarkerStyle(8);
    hsig->SetMarkerSize(1.5);
    hsig->Draw("PESAME");

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
    leg->AddEntry(hsig, "t#bar{t} events #pm 1#sigma (stat.+syst.)","pel");
    leg->AddEntry(hall, "Before background subtraction",            "pel");
    leg->AddEntry(hbkg, "Estimated background",                     "fl");
    leg->Draw();

    writeExtraText=false;
    CMS_lumi( c1, 2, 0, writeExtraText );
    c1->SaveAs((output+"/FinalACPSum.pdf").c_str());
}
