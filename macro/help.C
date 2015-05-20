#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMath.h>
void fix(TH1* histo) {

	Float_t val, errval;

	val=histo->GetBinContent(1)+histo->GetBinContent(0);
	errval=0;
	if(histo->GetBinContent(1)!=0.)
		errval+=TMath::Power(histo->GetBinError(1),2);
	if(histo->GetBinContent(0)!=0.)
		errval+=TMath::Power(histo->GetBinError(0),2);
	errval=TMath::Sqrt(errval);
	histo->SetBinContent(1,val);
	histo->SetBinError(1,errval);
	histo->SetBinContent(0,0);
	histo->SetBinError(0,0);

	Int_t highbin=histo->GetNbinsX();

	val=histo->GetBinContent(highbin)+histo->GetBinContent(highbin+1);
	errval=0;
	if(histo->GetBinContent(highbin)!=0.)
		errval+=TMath::Power(histo->GetBinError(highbin),2);
	if(histo->GetBinContent(highbin+1)!=0.)
		errval+=TMath::Power(histo->GetBinError(highbin+1),2);
	errval=TMath::Sqrt(errval);
	histo->SetBinContent(highbin,val);
	histo->SetBinError(highbin,errval);
	histo->SetBinContent(highbin+1,0);
	histo->SetBinError(highbin+1,0);
}

void variableRebinMTW(TH1* hin,TH1* hrebin) {

	for (Int_t ii=1;ii<=50;++ii) {
		hrebin->SetBinContent(ii,hin->GetBinContent(ii));
	}
}

void variableRebinEta(TH1* hin,TH1* hrebin) {

	for (Int_t ii=1;ii<=hin->GetNbinsX();++ii) {
		double xx = hin->GetBinCenter(ii);
		if (xx<2.5) hrebin->Fill(xx,hin->GetBinContent(ii));
	}
}

void variableRebin(TH1* hin,TH1* hrebin) {

	for (Int_t ii=1;ii<=hin->GetNbinsX();++ii) {
		hrebin->Fill(ii,hin->GetBinContent(ii));
	}
}

void beautifyNotFix(TH1* in, int color=kBlack, int style=0) {
	in->SetFillStyle(style);
	in->SetFillColor(color);
	in->SetLineColor(color);
}

void beautifyFillStyle(TH1* in, int color=kBlack, int style=0) {
	in->SetFillStyle(style);
	in->SetFillColor(color);
}

void beautify(TH1* in, int color=kBlack, int fillstyle=0, int linestyle=1, int linewidth=2) {
	in->SetFillStyle(fillstyle);
	in->SetFillColor(color);
	in->SetLineColor(color);
	in->SetLineStyle(linestyle);
	in->SetLineWidth(linewidth);
}

void beautifyAxis(TAxis* ax) {
	ax->SetLabelFont(132);
	ax->SetTitleFont(132);
	ax->SetLabelSize(0.07);
	ax->SetTitleSize(0.09);
}

void beautifyTopPad(TPad* pad) {
	pad->SetFillColor(0);
	pad->SetBorderMode(0);
	pad->SetBorderSize(2);
	pad->SetTopMargin(0.01);
	pad->SetBottomMargin(0.0);
	pad->SetFrameFillStyle(0);
	pad->SetFrameBorderMode(0);
	pad->SetFrameFillStyle(0);
	pad->SetFrameBorderMode(0);
}

void beautifyBottomPad(TPad* pad,TAxis* ax,TAxis* ay) {
	pad->SetFillColor(0);
	pad->SetBorderMode(0);
	pad->SetBorderSize(2);
	pad->SetTopMargin(0.0);
	pad->SetBottomMargin(0.52);
	pad->SetFrameFillStyle(0);
	pad->SetFrameBorderMode(0);
	pad->SetFrameFillStyle(0);
	pad->SetFrameBorderMode(0);
	pad->SetTicky();

	ax->SetLabelFont(132);
	ax->SetTitleFont(132);
	ax->SetLabelSize(0.18);
	ax->SetTitleSize(0.21);
	ax->SetTitleOffset(0.86);

	ay->SetNdivisions(3,4,0);
	ay->SetLabelFont(132);
	ay->SetTitleFont(132);
	ay->SetLabelSize(0.18);
	ay->SetTitleSize(0.20);
	ay->SetTitleOffset(0.34);

}
