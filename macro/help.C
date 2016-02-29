#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMath.h>
const int wtopx=50;
const int wtopy=50;
const int W = 800;
const int H = 600;
const float T = 0.08*H;
const float B = 0.12*H; 
const float L = 0.12*W;
const float R = 0.04*W;
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

void fixNewRange( TH1D* h, float xMin, float xMax)
{
    int b1_0 = h->GetXaxis()->GetFirst();
    int b2_0 = h->GetXaxis()->GetLast();
    float binWidth = h->GetBinWidth(1);
    float lowValue_0 = h->GetXaxis()->GetBinLowEdge(b1_0);
    float upValue_0 = h->GetXaxis()->GetBinUpEdge(b2_0);
    //if( xMin % binWidth != 0 || xMax % binWidth != 0 ) 
    //{
    //    std::cout<<">> [ERROR] Can't make range "<<xMin<<" - "<<xMax<<" with bin zise "<<binWidth<<std::endl;
    //    return;
    //}
    if( xMin > lowValue_0 )
    {
        float newcontent=0;
        float newsumw2=0;
        for( int b=b1_0; b<=b2_0; b++ )
        {
            if( h->GetXaxis()->GetBinLowEdge(b) < xMin )
            {
                newcontent += h->GetBinContent(b);
                newsumw2   += h->GetBinError(b)*h->GetBinError(b);
            }
            else
            {
                newcontent += h->GetBinContent(b);
                newsumw2   += h->GetBinError(b)*h->GetBinError(b);
                h->SetBinContent( b, newcontent    );
                h->SetBinError(   b, sqrt(newsumw2));
                break;
            }
        }
    }else
    { xMin = lowValue_0; }
    if( xMax < upValue_0 )
    {
        float newcontent=0;
        float newsumw2=0;
        for( int b=b2_0; b>=b1_0; b-- )
        {
            if( h->GetXaxis()->GetBinUpEdge(b) > xMax )
            {
                newcontent += h->GetBinContent(b);
                newsumw2   += h->GetBinError(b)*h->GetBinError(b);
            }
            else
            {
                newcontent += h->GetBinContent(b);
                newsumw2   += h->GetBinError(b)*h->GetBinError(b);
                h->SetBinContent( b, newcontent    );
                h->SetBinError(   b, sqrt(newsumw2));
                break;
            }
        }
    }else
    { xMax = upValue_0; }
    h->GetXaxis()->SetRangeUser(xMin, xMax); 
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
