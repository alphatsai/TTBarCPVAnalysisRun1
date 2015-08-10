#include "caculate.C"
#include <fstream>
void drawObs( TFile* f, 
        std::string hname="",
        std::string analysis="SemiLeptanic",
        std::string output=".",
        std::string ytitle="O",
        std::string xtitle="",
        int         rebinX=1,
        int         rebinY=20 )
{
    bool murmur=true;
    TH2D *h;
    if( analysis.compare("") != 0 )
    { 
        h = (TH2D*)f->Get((analysis+"/"+hname).c_str()); 
    }else{
        h = (TH2D*)f->Get(hname.c_str()); 
    }

    const int allh=6;

    h->Rebin2D( rebinX, rebinY );

    TCanvas *c1 = new TCanvas(("C2D_"+hname).c_str(), "c1",41,89,1213,664);
    gStyle->SetOptStat(0);
    c1->Range(-1.887097,-0.07984252,12.7379,0.07015748);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1290323);
    //c1->SetRightMargin(0.05045492);
    c1->SetTopMargin(0.06771654);
    c1->SetBottomMargin(0.1322835);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetLabelSize(0.05);
    //h->GetXaxis()->SetLabelSize(0.07);
    //h->GetXaxis()->SetLabelFont(62);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.83);
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitle(ytitle.c_str());
    h->Draw("COLZTEXT");

    double x0 = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
    double x1 = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast());

    TLine* line = new TLine( x0, 0, x1, 0 );
    //line->SetLineColor(2);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->Draw();

    TPaveText* t_title;
    t_title = new TPaveText(0.09842845,0.9387755,0.7278743,0.9843014,"brNDC");
    t_title->AddText(("CMS Simulation, L = 19.7/fb, #sqrt{s} = 8TeV ["+hname+"]").c_str());
    t_title->SetTextColor(kBlack);
    t_title->SetFillColor(kWhite);
    t_title->SetFillStyle(0);
    t_title->SetBorderSize(0);
    t_title->SetTextAlign(12);
    t_title->SetTextSize(0.04);
    t_title->Draw();

    c1->SaveAs((output+"/2D_"+hname+".pdf").c_str());
}
