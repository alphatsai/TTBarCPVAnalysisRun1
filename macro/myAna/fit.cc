#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooFitResult.h>

using namespace RooFit;

void plotDressing(TCanvas *canvas, RooPlot *frame)
{
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(2);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.035);
    canvas->SetTopMargin(0.07);
    canvas->SetBottomMargin(0.15);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderMode(0);
    
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitleSize(0.06);
    frame->GetXaxis()->SetTitleOffset(1.09);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetLabelSize(0.055);
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelOffset(0.01);
    frame->GetYaxis()->SetTitleOffset(1.20);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelSize(0.055);
}

RooFitResult *doDataFit(bool doDisplay = kTRUE)
{
    TFile *fin = new TFile("TemplateSyst_EvtChi2_Top_Leptonic_Mbl.root");
    TH1D* DATA_El = (TH1D*)fin->Get("DATA_El");
    TH1D* SigMC_El = (TH1D*)fin->Get("SigMC_El");
    TH1D* BkgMC_El = (TH1D*)fin->Get("BkgMC_El");
    
    RooRealVar mass("mass","",0.,500.);
    
    RooDataHist dh_DATA_El("dh_DATA_El","",mass,DATA_El);
    RooDataHist dh_SigMC_El("dh_SigMC_El","",mass,SigMC_El);
    RooDataHist dh_BkgMC_El("dh_BkgMC_El","",mass,BkgMC_El);
    
    RooHistPdf pdf_SigMC_El("pdf_SigMC_El","",mass,dh_SigMC_El);
    RooHistPdf pdf_BkgMC_El("pdf_BkgMC_El","",mass,dh_BkgMC_El);
    
    RooRealVar nsig("nsig","",30000.,0.,6E4);
    RooRealVar nbkg("nbkg","",2000.,0.,1E4);
    
    RooAddPdf model("model","",RooArgList(pdf_SigMC_El,pdf_BkgMC_El),RooArgList(nsig,nbkg));
    
    RooFitResult *res = model.fitTo(dh_DATA_El,Minos(kTRUE),Save(kTRUE),Extended(kTRUE));
    
    if (doDisplay) {
        RooPlot* frame_m = mass.frame();
        
        dh_DATA_El.plotOn(frame_m);
        model.plotOn(frame_m,LineColor(kBlack),LineWidth(2),LineStyle(kSolid));
        model.plotOn(frame_m,Components(pdf_SigMC_El),LineColor(kBlue),LineWidth(2),LineStyle(kDashed));
        model.plotOn(frame_m,Components(pdf_BkgMC_El),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed),DrawOption("F"));
        
        TCanvas *c1 = new TCanvas("c1","c1",768,568);

        plotDressing(c1,frame_m);
        frame_m->GetXaxis()->SetTitle("M(lb) [GeV]");
        frame_m->GetYaxis()->SetTitle("Entries / 10 GeV");
        frame_m->Draw();
    }
    
    return res;
}

RooFitResult *doToyFit(RooFitResult *res_mu)
{
    TFile *fin = new TFile("TemplateSyst_EvtChi2_Top_Leptonic_Mbl.root");
    TH1D* SigMC_El = (TH1D*)fin->Get("SigMC_El");
    TH1D* BkgMC_El = (TH1D*)fin->Get("BkgMC_El");
    
    RooRealVar mass("mass","",0.,500.);
    
    RooDataHist dh_SigMC_El("dh_SigMC_El","",mass,SigMC_El);
    RooDataHist dh_BkgMC_El("dh_BkgMC_El","",mass,BkgMC_El);
    
    RooHistPdf pdf_SigMC_El("pdf_SigMC_El","",mass,dh_SigMC_El);
    RooHistPdf pdf_BkgMC_El("pdf_BkgMC_El","",mass,dh_BkgMC_El);
    
    RooRealVar *nsig_mu = (RooRealVar*)res_mu->floatParsFinal().find("nsig");
    RooRealVar *nbkg_mu = (RooRealVar*)res_mu->floatParsFinal().find("nbkg");
    
    RooRealVar nsig("nsig","",nsig_mu->getVal(),0.,6E4);
    RooRealVar nbkg("nbkg","",nbkg_mu->getVal(),0.,1E4);
    
    RooAddPdf model("model","",RooArgList(pdf_SigMC_El,pdf_BkgMC_El),RooArgList(nsig,nbkg));
    
    RooDataHist *toy = model.generateBinned(RooArgSet(mass),Extended());
    
    RooFitResult *res = model.fitTo(*toy,Minos(kFALSE),Save(kTRUE),Extended(kTRUE));
    
    fin->Close();
    delete fin;
    
    return res;
}

void fit()
{
    RooFitResult *res_mu = doDataFit(kTRUE);
    
    TNtupleD *nt_toy = new TNtupleD("nt_toy","","pull");
    
    for(int idx=0;idx<1000;idx++) {
        RooFitResult *res = doToyFit(res_mu);
        
        RooRealVar *nsig_mu = (RooRealVar*)res_mu->floatParsFinal().find("nsig");
        RooRealVar *nsig = (RooRealVar*)res->floatParsFinal().find("nsig");
        
        double pull = (nsig->getVal()-nsig_mu->getVal())/nsig->getError();
        nt_toy->Fill(&pull);
        
        delete res;
    }
    
    RooRealVar pull("pull","",-5.,5.);
    
    RooRealVar mean("mean","", 0.,-1.,1.);
    RooRealVar width("width","",1.,0.8,1.2);
    RooGaussian gaussian("gaussian","",pull,mean,width);
    
    RooDataSet *toy = new RooDataSet("toy","",nt_toy,RooArgSet(pull));
    RooFitResult *res_pull = gaussian.fitTo(*toy,Minos(kTRUE),Save(kTRUE));
    
    RooPlot* frame_p = pull.frame();
    toy->plotOn(frame_p,Binning(50));
    gaussian.plotOn(frame_p,LineWidth(3));
    
    TCanvas *c2 = new TCanvas("c2","c2",768,568);
    plotDressing(c2,frame_p);
    frame_p->GetXaxis()->SetTitle("Pull [#sigma]");
    frame_p->GetYaxis()->SetTitle("Entries");
    frame_p->Draw();
    
    //print out the results
    res_mu->Print("v");
    res_pull->Print("v");
}