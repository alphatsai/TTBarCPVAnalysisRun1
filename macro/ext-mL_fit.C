#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <TMinuit.h>
#include "vector.h"
#include <TMath.h>
#include "TVirtualFitter.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TLine.h"
#include "TRandom2.h"


using namespace std;

#define NPAR 2

TF1* fchi2;
//TH1F* h1;
TF1* bsFin;

// Remember to change the values of xh and xl in all_tree.C too!!!
//Double_t xh = 6.5;
//Double_t xl = 4.5;
Double_t xl = -1.;
Double_t xh = 1.;
Double_t bin_size = 0.01;

const double _two_pi = 2.0 * TMath::Pi();
Double_t fit_lo_edge = -1.;
Double_t fit_hi_edge = 1.;

Int_t toymc_sig=13084;
Int_t toymc_bkg=12504;

vector<Double_t> dataColl;
vector<Double_t> sigColl;
vector<Double_t> bkgColl;

vector<Double_t> totalColl;
vector<Double_t> ctauColl;

vector<Double_t> info;
vector<Double_t> info_err;

//par[0] = fs Jpsi signal fraction;
//-----mass part-----------------------------------------------------
//par[1] = g norm; par[2] g1 mean; par[3] g1 width; 
//par[4] g2 ratio;  par[5] g2 mean; par[6] g2 width;
//par[7] bg norm; par[8] bg slope;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  Double_t Nevt=0.;
  Double_t fs = par[0];
  Double_t fb = par[1];

  for ( int i=0; i<dataColl.size(); i++ ) {
    Nevt += dataColl[i];
    //PDF for signal and background
    Double_t Ls = sigColl[i];
    Double_t Lb = bkgColl[i];	
    for (int data=0; data<dataColl[i]; data++) {	
      //Get Log Likelihood
      //if(Ls!=0. ||  Lb!=0.) Lsum += TMath::Log( (fs*Ls + fb*Lb) / (fs+fb) );
      if(Ls!=0. &&  Lb!=0.) Lsum += TMath::Log( (fs*Ls + fb*Lb) / (fs+fb) );
    }
  }
  f=2*( -1*Lsum + (fs+fb)  - Nevt*TMath::Log(fs+fb) );
  
}


//___________________________________________________________________________
Double_t* Ifit( TFile* fin, std::string name, std::string output=".", std::string xTitle="", std::string yTitle="Events", int rebin=1, int channel=0, int fit_data=1 )
{
    Double_t* fitted = new Double_t[8];

    TCanvas *c1 = new TCanvas("HF1", "Histos1", 258,92,748,702);
    c1->Range(-104.4905,-2560.33,537.9965,11563.46);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1626344);
    c1->SetRightMargin(0.05913978);
    c1->SetTopMargin(0.05349183);
    c1->SetBottomMargin(0.1812779);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

    double count=0;
    dataColl.clear();
    sigColl.clear();
    bkgColl.clear();

    totalColl.clear();
    ctauColl.clear();

    //Get data from looping tree

    //TFile *fin = new TFile("results/15Dec_LepJet_MCDATA/Template_EvtChi2_Top_Hadronic_Mass.root");
    //TFile *fin = new TFile("results/15Dec_LepJet_MCDATA/Template_EvtChi2_Top_Leptonic_Mbl.root");


    //   TH1D *hsig = new TH1D();
    //   TH1D *hbkg = new TH1D();
    TH1D *hsig_toymc = new TH1D();
    TH1D *hbkg_toymc = new TH1D();

    char hname[30];
    std::string ch;
    if( channel == 1 )
        ch="_El";
    else if( channel == 2)
        ch="_Mu";
    else
        ch="";

    TH1D * hsig = (TH1D*)((TH1D*)fin->Get(("SigMC"+ch).c_str()))->Clone();
    TH1D * hbkg = (TH1D*)((TH1D*)fin->Get(("BkgMC"+ch).c_str()))->Clone();
    TH1D *hEGdata;
    hEGdata = (TH1D*)((TH1D*)fin->Get(("DATA"+ch).c_str()))->Clone();

    // hsig->Sumw2();
    //hbkg->Sumw2();

    // hsig->Rebin(2);
    // hbkg->Rebin(2);
    // hEGdata->Rebin(2);
    // hbkg_toymc->Rebin(2);
    // hsig_toymc->Rebin(2);

    hsig->Rebin(rebin);
    hbkg->Rebin(rebin);
    hEGdata->Rebin(rebin);

    // normalize template
    hsig->Scale(1./hsig->Integral());
    hbkg->Scale(1./hbkg->Integral());  
    if(fit_data==0){
        hsig_toymc->Scale(1./hsig_toymc->Integral());
        hbkg_toymc->Scale(1./hbkg_toymc->Integral());
    }

  TH1D *hsum = new TH1D();
  int ntemplate = 1000.;
  float sigfrac = 0.5;
  TH1D *hsum_norm = new TH1D();
  TH1D *hdata = new TH1D();

  int ndata=0;
  if ( fit_data>0 ) {
    hdata = (TH1D*)hEGdata->Clone();
    ndata = hdata->GetEntries();
  }else { //generate toymc
    hsum = (TH1D*)hsig_toymc->Clone();
    hsum->Scale(toymc_sig);
    hsum->Add(hbkg_toymc,toymc_bkg);
    
    hsum_norm = (TH1D*)hsum->Clone();  
    hsum_norm->Scale(1./hsum->Integral());
    hdata = (TH1D*)hsum_norm->Clone();
    //ndata = (int) gRandom->Poisson(hsum->Integral());
    ndata=toymc_sig+toymc_bkg;
    hdata->FillRandom(hsum_norm, ndata);
  }
  if(ndata==0) {
    printf(" ---  no events in the fit \n");
    fitted[0] = 0.;
    fitted[1] = 0.;
    fitted[2] = 0.;
    fitted[3] = 0.;
    fitted[4] = 0.;
    fitted[5] = 0.;
    fitted[6] = 0.;
    fitted[7] = 0.;
    fin_data->Close();
    fin->Close();
    fin_gjet6000->Close();

    return fitted;
  }    

  printf(" --------- before the fit ------------- \n");
  printf("Nsig %2.3f, Nbg %2.3f, Ntemplate %3.3f \n", hsig->Integral(), hbkg->Integral(), ntemplate);
  printf("Purity %2.3f, init size %4.3f,  test sample size %4d\n", hsig->Integral()/hsum->Integral(), hsum->Integral(), ndata);
  printf(" -------------------------------------- \n");

  int nbins = hdata->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++) {
    dataColl.push_back(hdata->GetBinContent(ibin));
    sigColl.push_back(hsig->GetBinContent(ibin));
    bkgColl.push_back(hbkg->GetBinContent(ibin));    
  }
  printf( " -----  Got %d, %d, %d events for fit ----- \n ", dataColl.size(), sigColl.size(), bkgColl.size() );  
  if ( dataColl.size() != sigColl.size() || sigColl.size()!=bkgColl.size() ) {
    printf(" error ...  inconsistent hit collection size \n");
    fin_data->Close();
    fin->Close();
    fin_gjet6000->Close();

    return fitted;
  }

  //--------------------------------------------------
  //init parameters for fit
  Double_t vstart[10] = {1., 1.};
  vstart[0] = sigfrac*ndata;
  vstart[1] = (1-sigfrac)*ndata;
 
  TMinuit *gMinuit = new TMinuit(NPAR);  
  gMinuit->Command("SET STR 1");
  gMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  arglist[0] = 1;
  gMinuit->mnexcm("SET PRINT", arglist ,1,ierflg);

  Double_t step[] = { 0.1, 0.1,};

  gMinuit->mnparm(0,  "Signal yield"  , vstart[0],  step[0], 0., ndata*2.  , ierflg);
  gMinuit->mnparm(1,  "background yield"  , vstart[1],  step[1], 0., ndata*2. , ierflg);
  
  printf(" --------------------------------------------------------- \n");
  printf(" Now ready for minimization step \n --------------------------------------------------------- \n");
  
  arglist[0] = 2000; // number of iteration
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  printf (" -------------------------------------------- \n");
  printf("Finished.  ierr = %2.2f \n", ierflg);

  info.clear();
  info_err.clear();

  double para[NPAR+1],errpara[NPAR+1];
  if ( ierflg == 0 ) 
    {
      for(int j=0; j<=NPAR-1;j++) {
        gMinuit->GetParameter(j, para[j],errpara[j]);
        para[NPAR] = dataColl.size();
        info.push_back(para[j]);
        info_err.push_back(errpara[j]);
        printf("Parameter (yeild) %d = %f +- %f\n",j,para[j],errpara[j]);
	
      }
      printf(" fitted yield %2.3f \n", (para[0]+para[1])/ndata );

      info.push_back(sigColl.size());

      //do minos if fit sucessed.
//       printf("         ---------------------------------------------------------\n");
//       printf("          Now call for minos step \n");
//       printf("         ---------------------------------------------------------\n");
      
//       arglist[0] = 200; // number of iteration
//       arglist[1] = 1;
//       gMinuit->mnexcm("MINOS", arglist ,2,ierflg);
//       printf("         --------------------------------------------------------- \n");
//       printf("         Done Minos.  ierr = %d \n", ierflg);
//       Double_t amin;
//       gMinuit->mnprin(1,amin);
    }
  else {
    printf(" *********** Fit failed! ************\n");
    gMinuit->GetParameter(0, para[0],errpara[0]);
    gMinuit->GetParameter(1, para[1],errpara[1]);
    para[0]=0.; errpara[0]=0.;
  }

  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(1,amin);  
  gMinuit->mnmatu(1);
  printf(" ========= happy ending !? =========================== \n");
  
  printf("FCN =  %3.3f \n", amin);

  double yerr[100];
  for(int i=0;i<100;i++){
    yerr[i] = 0.;
  }

  hsig->Scale(para[0]);
  hbkg->Scale(para[1]);
  TH1D *hfit = (TH1D*)hbkg->Clone();
  hfit->Add(hsig);

  hsig->SetLineColor(4);
  hsig->SetLineWidth(2);
//   hsig->SetFillColor(5);
//   hsig->SetFillStyle(3001);

//   hbkg->SetLineWidth(2);
  // plot
  c1->Draw();  
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0); 
  //gPad->SetLogy();
  hdata->SetLineColor(1);

  hdata->SetXTitle(xTitle.c_str());
  hdata->SetYTitle(yTitle.c_str());
  hdata->SetTitle("");
  hdata->SetMarkerStyle(8);
  hdata->SetMinimum(0.);
  hdata->GetXaxis()->SetNdivisions(505);
  hdata->GetXaxis()->SetLabelFont(42);
  hdata->GetXaxis()->SetLabelSize(0.05);
  hdata->GetXaxis()->SetTitleSize(0.06);
  hdata->GetXaxis()->SetTitleOffset(1.15);
  hdata->GetXaxis()->SetTitleFont(42);
  hdata->GetYaxis()->SetNdivisions(505);
  hdata->GetYaxis()->SetLabelFont(42);
  hdata->GetYaxis()->SetLabelSize(0.035);
  hdata->GetYaxis()->SetTitleSize(0.06);
  hdata->GetYaxis()->SetTitleOffset(1.21);
  hdata->GetYaxis()->SetTitleFont(42);

  float ymax = hdata->GetMaximum();
  if ( hfit->GetMaximum() > hdata->GetMaximum() ) ymax = hfit->GetMaximum();
  if ( hdata->GetMaximum() < 15 ) ymax = 15;
  hdata->SetMaximum(ymax*1.4);
  hfit->SetMaximum(ymax*1.4);
  hsig->SetMaximum(ymax*1.4);
  hbkg->SetMaximum(ymax*1.4);

  hdata->Draw("p e");

  hbkg->SetMarkerStyle(0);
  hbkg->SetFillColor(2);
  hbkg->SetLineWidth(1);
  hbkg->SetLineColor(2);
  hbkg->SetFillStyle(3005);
  hbkg->SetError(yerr);
  hbkg->Draw("h same");

  hsig->SetMarkerStyle(0);
  hsig->SetError(yerr);
  hsig->Draw("h same");

  hfit->SetMarkerStyle(0);
  hfit->SetLineColor(1);
  hfit->SetLineWidth(2);
  hfit->SetError(yerr);
  //printf("nbins hfit %d \n", hfit->GetNbinsX());
  hfit->Draw("h same");
  hdata->Draw("p e same");
  
  TLegend *tleg = new TLegend(0.5241935,0.6344725,0.8682796,0.9331352,NULL,"brNDC");
  char text[50];
  sprintf(text,"Top Mass");
  //tleg->SetHeader(text);
  tleg->SetBorderSize(0);
  tleg->SetTextSize(0.03120357);
  tleg->SetLineColor(1);
  tleg->SetLineStyle(1);
  tleg->SetLineWidth(1);
  tleg->SetFillColor(0);
  tleg->SetFillStyle(0);
  sprintf(text,"Data %5.1f events",hdata->Integral());
  tleg->AddEntry(hdata,text,"pl");
  sprintf(text,"Fitted %5.1f events",hfit->Integral());
  
  tleg->AddEntry(hfit,text,"l");
  sprintf(text,"SIG %5.1f #pm %5.1f events",para[0], errpara[0]);
  tleg->AddEntry(hsig,text,"f");
  sprintf(text,"BKG %5.1f #pm %5.1f events",para[1], errpara[1]);
  
  tleg->AddEntry(hbkg,text,"f");
  tleg->Draw();
  
  TLatex *tlx = new TLatex(6.247421e-06,9218.143,"CMS #sqrt{s} = 8TeV, L=19.7 fb^{-1}");
  tlx->SetTextSize(0.035);
  tlx->SetLineWidth(2);
  tlx->Draw();

  //gPad->RedrawAxis();
  
  if(fit_data>0) hdata->Chi2Test(hfit,"P");

  c1->SaveAs((output+"/FittingResults_"+name+ch+".pdf").c_str());
  return fitted;

//   float sig_part = hsig->Integral(ibin1,hfit->GetNbinsX());
//   float sig_part_err = hsig->Integral(ibin1,hfit->GetNbinsX())*errpara[0]/para[0];
//   float bkg_part = hbkg->Integral(ibin1,hfit->GetNbinsX());
//   float bkg_part_err = hbkg->Integral(ibin1,hfit->GetNbinsX())*errpara[1]/para[1];
//   printf("%s Data %5.1f events, fitted %5.1f\n", EBEE, hdata->Integral(), hfit->Integral());
//   printf("%s Data %5.1f, and fitted (in 5GeV) %5.1f events \n", EBEE, hdata->Integral(ibin1,hfit->GetNbinsX()), hfit->Integral(ibin1,hfit->GetNbinsX()));
//   printf("%s SIG %5.1f #pm %5.1f events \n", EBEE, para[0], errpara[0]);
//   printf("%s SIG (in 5GeV) %5.1f #pm %5.1f events \n", EBEE, sig_part, sig_part_err);
//   printf("%s BKG %5.1f #pm %5.1f events \n", EBEE, para[1], errpara[1]);
//   printf("%s BKG (in 5GeV) %5.1f #pm %5.1f events \n", EBEE, bkg_part, bkg_part_err);
   

//   char fname[30];
//   sprintf(fname,"plots/test_Ifit%s_%d_%d.pdf",EBEE, jetbin, ptbin);
  
//   printf("----- fit results with signal projection   ----------- \n");
//   if(fit_data>0) hdata->Chi2Test(hfit,"P");
//   //ftemplate->Close();

//   fitted[0] = para[0];
//   fitted[1] = errpara[0]/TMath::Sqrt(2);
//   fitted[2] = para[1];
//   if (fit_data==2 ) fitted[2] += hdata->GetBinContent(hdata->GetNbinsX()+1);
//   fitted[3] = errpara[1]/TMath::Sqrt(2);
//   fitted[4] = sig_part;
//   fitted[5] = sig_part_err/TMath::Sqrt(2);
//   fitted[6] = bkg_part;
//   fitted[7] = bkg_part_err/TMath::Sqrt(2);
  
//   if(fit_data==0){
//     fin_filter->Close();
//     fin_data->Close();
//     fin->Close();
//     fin_gjet6000->Close();
//     fin_DYMC->Close();
//     fin_DYData->Close();
//     fin_WJetMC->Close();
//     fin_WJetData->Close();
//     fin_WJetTemplate->Close();
//     fin_WJetTemplate_alt->Close();
//   }

  return fitted;
}

void pulltest(int ptbin=13, char EBEE[10]="EB", int jetbin=0, int sig=400, int bkg=400){

  toymc_sig = sig;
  toymc_bkg = bkg;

  double *fitted;
  TH1F *h1 = new TH1F("h1","",100,-10., 10.);
  TH1F *h2 = new TH1F("h2","",40000, 0., 40000);

  h1->SetNdivisions(505,"XY");
  h2->SetNdivisions(505,"XY");

  int nexp=200;
  Double_t Nevt=0.;

  for (int i=0; i<nexp; i++) {
    fitted = Ifit(ptbin,EBEE, 0, jetbin);
     Nevt=sig+bkg;
    // for ( int ii=0; ii<dataColl.size(); ii++ ) {
    //   Nevt += dataColl[ii];
    // }
    // printf("fit purity %2.2f +- %2.2f err with %d events. \n", info[0], info_err[0], Nevt);
    // h1->Fill((info[0]/Nevt-input)/(info_err[0]/Nevt));
    printf("input %d, %d, fitted %.2f %.2f \n", toymc_sig, toymc_bkg, fitted[0], fitted[1]);
    h1->Fill((fitted[0]-toymc_sig)/fitted[1]);
    h2->Fill(fitted[0]);
  }    

  TCanvas *c2 = new TCanvas("c2","",1000,500);
  c2->Divide(2,1);
  c2->cd(1);
  char txt[100];
  sprintf(txt, "(fitted-gen)/error");
  h1->SetXTitle(txt);
  h1->Fit("gaus");
  //h1->GetXaxis()->SetRangeUser(-5.,5.);
  h1->Draw();
  c2->cd(2);
  sprintf(txt, "fitted signal (input %d)", toymc_sig);
  h2->SetXTitle(txt);
  h2->Fit("gaus");
  h2->GetXaxis()->SetRangeUser(toymc_sig*0.5, toymc_sig*1.5);
  //if ( input >0.8 )  h2->GetXaxis()->SetRangeUser(0., Nevt*1.4);
  h2->Draw();  
  sprintf(txt, "plots/extmLfit_pull_%s_%d_pt%d.pdf", EBEE, jetbin,  ptbin);
  c2->SaveAs(txt);
  sprintf(txt, "plots/extmLfit_pull_%s_%d_pt%d.C", EBEE, jetbin,  ptbin);
  c2->SaveAs(txt);

  TF1 *func = (TF1*)h2->GetFunction("gaus");
  printf(" TOYTOY pt %d, mean %.2f \n", ptbin, func->GetParameter(1));
  printf(" %5.2f \\\% \\pm %5.2f \\\% \n", 100.-(func->GetParameter(1)/sig*100.), TMath::Sqrt(sig)/sig*100.);
  
}

void Draw_yield_treeeff(char EBEE[20]="EB", int jetbin=0){

  int ptbin_int=0;
  float ptcut[30] = {22, 30, 36, 50, 75, 90, 105,  120, 135, 150, 165, 175, 185,
		     190, 200, 220, 250, 300, 350, 400, 500, 750, 1000, 1500, 2000, 3000, 10000}; //22 bins
  //                  13   14   15   16   17   18   19   20   21    22    23    24    25     26

  int nbin=22;

  TH1F *h_yield = new TH1F("h_yield","",nbin, ptcut);
  TH1F *h_purity = new TH1F("h_purity","",nbin, ptcut);
  TH1F *h_purity_tight = new TH1F("h_purity_tight","",nbin, ptcut);

  TH1F *h_sig_yield = new TH1F("h_sig_yield","",nbin, ptcut);
  TH1F *h_sig_yield_tight = new TH1F("h_sig_yield_tight","",nbin, ptcut);
  TH1F *h_bkg_yield = new TH1F("h_bkg_yield","",nbin, ptcut);
  TH1F *h_bkg_yield_tight = new TH1F("h_bkg_yield_tight","",nbin, ptcut);


  TH1F *h_xs = new TH1F("h_xs","",nbin,ptcut);
  double *fitted;
  float lumi = 2568.83;
  float lumi_err = lumi*0.046;
  float deta = 1.4442*2.; 
  float template_sys = TMath::Sqrt(3.2*3.2+3.*3.)*0.01;
  int ebeebin=0;
  if(strcmp(EBEE,"EE")==0) {
    deta = (2.5-1.566)*2.;
    ebeebin=1;
  }

  TTree *tt = new TTree();
  tt->ReadFile("eff.dat");
  Int_t   ptbin_;
  Int_t   ebee_;
  Int_t   jetbin_;
  Float_t trigeff;
  Float_t trigeff_err;
  Float_t recoeff;
  Float_t recoeff_err;
  Float_t preseleff;
  Float_t preseleff_err;
  Float_t SF;
  Float_t SF_err;       

  tt->SetBranchAddress("ptbin", &ptbin_);
  tt->SetBranchAddress("EBEE", &ebee_);
  tt->SetBranchAddress("jetbin", &jetbin_);
  tt->SetBranchAddress("trigeff", &trigeff);
  tt->SetBranchAddress("trigeff_err", &trigeff_err);
  tt->SetBranchAddress("recoeff", &recoeff);
  tt->SetBranchAddress("recoeff_err", &recoeff_err);
  tt->SetBranchAddress("preseleff", &preseleff);
  tt->SetBranchAddress("preseleff_err", &preseleff_err);
  tt->SetBranchAddress("SF", &SF);
  tt->SetBranchAddress("SF_err", &SF_err);

  char txt[100];
  for(int ii=13; ii<22; ii++){
    //perform fit for yield
    fitted=Ifit(ii, EBEE,1, jetbin);
    if(fitted[0]>0.) {

      h_yield->SetBinContent(ii+1, fitted[0]/h_yield->GetBinWidth(ii+1));
      h_yield->SetBinError(ii+1, fitted[1]/h_yield->GetBinWidth(ii+1));

      Long64_t jentry = ii-13 + 9*3*ebeebin + 9*jetbin;
      tt->GetEntry(jentry);
      printf("bin %d, %d, %d \n", ptbin_, ebee_, jetbin_);
    
      float djet_eta=1.5*2.;
      if(jetbin==1) djet_eta = 0.9*2.;
      if(jetbin==2) djet_eta = 1;
      printf("bin %d, ptcut %.1f,  fit %.0f , bkg %.0f , eff %.2f, %2.f, binwidth %.1f, deta %.2f \n", ii, ptcut[ii], fitted[0], fitted[2], recoeff, preseleff, h_yield->GetBinWidth(ii+1), deta );
      printf(" %f %f %f %f %f %f %f %f %f \n",  fitted[0] , lumi , trigeff, recoeff , preseleff , h_yield->GetBinWidth(ii+1) ,SF ,deta , djet_eta);
      float xs = fitted[0] / lumi / trigeff/ recoeff / preseleff / h_yield->GetBinWidth(ii+1) /SF /deta / djet_eta ; //xs per GeV
      float xs_err = (fitted[1]/fitted[0])*(fitted[1]/fitted[0]) + (lumi_err/lumi)*(lumi_err/lumi) + 
	(preseleff_err/preseleff)*(preseleff_err/preseleff) + 
	(recoeff_err/recoeff)*(recoeff_err/recoeff) + 
	(trigeff_err/trigeff)*(trigeff_err/trigeff) +
	(SF_err/SF)*(SF_err/SF) +
	(template_sys*template_sys);
      xs_err = TMath::Sqrt(xs_err)*xs;
      printf("xs %f , xs_err %f \n", xs, xs_err);
      h_xs->SetBinContent(ii+1, xs);
      h_xs->SetBinError(ii+1, xs_err);

      h_sig_yield->SetBinContent(ii+1, fitted[0]);
      h_bkg_yield->SetBinContent(ii+1, fitted[2]);

      h_sig_yield_tight->SetBinContent(ii+1, fitted[4]);
      h_bkg_yield_tight->SetBinContent(ii+1, fitted[6]);

      // //perform fit for purity in WP90 region
      // h_purity_tight->SetBinContent(ii+1, fitted[4]/(fitted[4]+fitted[6]));
      // float err = TMath::Sqrt(fitted[4]*fitted[6]/(fitted[4]+fitted[6])/(fitted[4]+fitted[6])/(fitted[4]+fitted[6]));
      // h_purity_tight->SetBinError(ii+1, err);    
      // h_purity->SetBinContent(ii+1, fitted[0]/(fitted[0]+fitted[2]));
      // err = TMath::Sqrt(fitted[0]*fitted[2]/(fitted[0]+fitted[2])/(fitted[0]+fitted[2])/(fitted[0]+fitted[2]));
      // h_purity->SetBinError(ii+1, err);    
    }
  }
  for(int ii=13; ii<nbin; ii++){
    //   printf("ptbin %d, %.2f, xs %f , xs_err %f \n", ii, ptcut[ii], h_xs->GetBinContent(ii+1), h_xs->GetBinError(ii+1));
    printf("ptbin %d, %.2f, yield, %.2f, %.2f \n", ii, ptcut[ii], h_sig_yield_tight->GetBinContent(ii+1), h_bkg_yield_tight->GetBinContent(ii+1));
  }

  TCanvas *c10 = new TCanvas("c10","",600,600);
  c10->Draw();
  gPad->SetLogy();
  h_yield->SetNdivisions(505,"XY");
  h_yield->SetXTitle("p_{T} (GeV)");
  h_yield->SetYTitle("Entries / GeV");
  h_yield->SetMarkerStyle(8);
  h_yield->GetXaxis()->SetRangeUser(150.,1000.);
  h_yield->Draw("pe");

  char pho_text[100];
  char jet_text[100];
  if(ebeebin==0) sprintf(pho_text,"|#eta_{#gamma}|<1.4442");
  else sprintf(pho_text,"1.566<|#eta_{#gamma}|<2.5");
  if(jetbin==0) sprintf(jet_text,"|#eta_{jet}|<1.5");
  else sprintf(jet_text,"1.5<|#eta_{jet}|<2.4");


  TLegend *tleg = new TLegend(0.4, 0.65, 0.85, 0.85);
  char text[50];
  sprintf(text,"CMS 13TeV, %.0f pb^{-1}",lumi);
  tleg->SetHeader(text);
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);  
  sprintf(text,"%s, %s",pho_text, jet_text);
  if(jetbin==2)   sprintf(text,"%s",pho_text);
  tleg->AddEntry(h_yield,text,"pl");
  tleg->Draw();
  

  TLatex *tlx = new TLatex();
  tlx->SetTextSize(0.04);
  tlx->DrawLatex(100, h_yield->GetMaximum()*1.3, "CMS Preliminary #sqrt{s} = 13TeV");

  // //purity plot
  // TCanvas *c11 = new TCanvas("c11","",600,600);
  // h_purity->SetNdivisions(505,"XY");
  // h_purity->SetXTitle("P_{T} (GeV)");
  // h_purity->SetYTitle("Purity");
  // h_purity->SetMarkerStyle(8);
  // h_purity->GetXaxis()->SetRangeUser(150.,1000.);
  // h_purity->GetYaxis()->SetRangeUser(0., 1.01);
  // h_purity->Draw();
  // // h_purity_tight->SetMarkerStyle(8);
  // // h_purity_tight->SetMarkerColor(4);
  // // h_purity_tight->SetLineColor(4);
  // // h_purity_tight->Draw("pe same");

  
  // TGraphAsymmErrors *tgrs = new TGraphAsymmErrors();
  // h_bkg_yield->Add(h_sig_yield);

  // tgrs->BayesDivide(h_sig_yield, h_bkg_yield);
  // tgrs->SetMarkerStyle(8);
  // tgrs->Draw("pe same");

  // TGraphAsymmErrors *tgrs_tight = new TGraphAsymmErrors();
  // h_bkg_yield_tight->Add(h_sig_yield_tight);
  // for(int ii=13; ii<nbin; ii++){
  //   //   printf("ptbin %d, %.2f, xs %f , xs_err %f \n", ii, ptcut[ii], h_xs->GetBinContent(ii+1), h_xs->GetBinError(ii+1));
  //   printf("ptbin %d, %.2f, yield, %.2f, %.2f \n", ii, ptcut[ii], h_sig_yield_tight->GetBinContent(ii+1), h_bkg_yield_tight->GetBinContent(ii+1));
  // }

  // tgrs_tight->BayesDivide(h_sig_yield_tight, h_bkg_yield_tight,"v");
  // tgrs_tight->SetMarkerStyle(8);
  // tgrs_tight->SetMarkerColor(4);
  // tgrs_tight->SetLineColor(4);
  // tgrs_tight->Draw("pe same");

  // tleg = new TLegend(0.3, 0.15, 0.85, 0.35);
  // char text[50];
  // sprintf(text,"Single Photon PD %.0f pb^{-1}", lumi);
  // tleg->SetHeader(text);
  // tleg->SetFillColor(0);
  // tleg->SetShadowColor(0);
  // tleg->SetBorderSize(0);
  // sprintf(text,"%s, %s",pho_text, jet_text);
  // tleg->AddEntry(tgrs,text,"pl");
  // sprintf(text,"%s, %s, BDT>0.37",pho_text, jet_text);
  // tleg->AddEntry(tgrs_tight,text,"pl");
  // tleg->Draw();

  // tlx->SetTextSize(0.04);
  // tlx->DrawLatex(200, 1.01, "CMS Preliminary #sqrt{s} = 13TeV");
  

  //Draw XS plot
  TCanvas *c12 = new TCanvas("c12","",600,700);  
  c12->Draw();

  TPad* pad1 = new TPad("pad1","",0.02, 0.25, 0.99, 0.99);
  TPad* pad2 = new TPad("pad1","",0.02, 0.02, 0.99, 0.25);


  // pad1->SetLeftMargin(0.02);
  pad1->SetRightMargin(0.035);
  // pad1->SetTopMargin(0.02);
  pad1->SetBottomMargin(0.0);
  pad1->SetFrameBorderMode(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);

  // pad2->SetLeftMargin(0.02);
  pad2->SetRightMargin(0.035);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.275);
  pad2->SetFrameBorderMode(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);

  pad1->Draw();
  pad2->Draw();
  
  pad1->cd();

  gPad->SetLogy();

  TH2F *hh2 = new TH2F("hh2","",100, 150, 1000, 100, 2e-7, 1.);
  if(ebeebin==0) hh2 = new TH2F("hh2","",100, 150, 1000, 100, 2e-7, 1.);

  hh2->SetNdivisions(505,"XY");
  hh2->SetXTitle("P_{T} (GeV)");
  hh2->SetYTitle("d^{3}#sigma / d^{#gamma}p_{T}d^{#gamma}#eta d^{jet}#eta [pb/GeV]");
  if(jetbin==2)   hh2->SetYTitle("d^{2}#sigma / d^{#gamma}p_{T}d^{#gamma}#eta [pb/GeV]");
  hh2->GetXaxis()->SetLabelSize(0.035);
  hh2->GetYaxis()->SetLabelSize(0.035);
  hh2->GetYaxis()->SetTitleSize(0.04);
  hh2->GetYaxis()->SetTitleOffset(1.1);
  hh2->Draw();
  
  h_xs->SetMarkerStyle(8);
  h_xs->Draw("pel same");
  //Draw JetPhox value
  TTree *t_th = new TTree();
  t_th->ReadFile("jetphox.dat");
  t_th->Print();

  Int_t ptbin_1;
  Int_t ebee_1;
  Int_t jetbin_1;
  Double_t xs_th;
  Double_t xs_th_err;
  
  
  t_th->SetBranchAddress("ptbin", &ptbin_1);
  t_th->SetBranchAddress("EBEE", &ebee_1);
  t_th->SetBranchAddress("jetbin", &jetbin_1);
  t_th->SetBranchAddress("xs", &xs_th);
  t_th->SetBranchAddress("xs_err", &xs_th_err);
  
  
  TH1F *h_th_xs = new TH1F("h_th_xs","",nbin, ptcut);
  TH1F *h_th_xs_err = new TH1F("h_th_xs_err","",nbin, ptcut);

  for(int ii=13; ii<22; ii++){
    Long64_t jentry = ii-13 + 9*3*ebeebin + 9*jetbin;
    t_th->GetEntry(jentry);
    //printf("bin %d, %d, %d \n", ptbin_1, ebee_1, jetbin_1);
    h_th_xs->SetBinContent(ii+1, xs_th);
    h_th_xs->SetBinError(ii+1, xs_th_err);
    h_th_xs_err->SetBinContent(ii+1, xs_th_err);

    if(jetbin==2){
      if(ebeebin==0){
	h_th_xs->SetBinContent(ii+1, xs_th/deta);
	h_th_xs->SetBinError(ii+1, xs_th_err/deta);
	h_th_xs_err->SetBinContent(ii+1, xs_th_err/deta);
      }else{
	h_th_xs->SetBinContent(ii+1, xs_th/deta*2.);
	h_th_xs->SetBinError(ii+1, xs_th_err/deta*2.);
	h_th_xs_err->SetBinContent(ii+1, xs_th_err/deta*2.);
      }
    }
  }

  h_th_xs->SetMarkerStyle(25);
  h_th_xs->SetMarkerColor(2);
  h_th_xs->SetLineColor(2);
  h_th_xs->Draw("same ple");

  TLegend *tleg = new TLegend(0.4, 0.65, 0.85, 0.85);
  char text[50];
  sprintf(text,"CMS 13TeV, %.0f pb^{-1}",lumi);
  tleg->SetHeader(text);
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);  
  sprintf(text,"%s, %s",pho_text, jet_text);
  if(jetbin==2)   sprintf(text,"%s",pho_text);
  tleg->AddEntry(h_xs,text,"pl");
  tleg->AddEntry(h_th_xs,"JETPHOX","pl");
  tleg->Draw();
  
  tlx->DrawLatex(150, 1., "CMS Preliminary #sqrt{s} = 13TeV");

  
  //Draw data/th XS ratio
  // TCanvas *c15 = new TCanvas("c15","",600,600);
  // c15->Draw();
  pad2->cd();
  gPad->SetGridy();
  TH1F *h_xs_ratio = (TH1F*)h_xs->Clone();
  h_xs_ratio->Draw();


  h_xs_ratio->Divide(h_th_xs);
  // TGraphAsymmErrors *tgrs = new TGraphAsymmErrors();
  // tgrs->Divide(h_xs, h_th_xs);
  // TH1F *h_xs_ratio = new TH1F("h_xs_ratio","",nbin,ptcut);
  h_xs_ratio->SetMinimum(0.);
  h_xs_ratio->SetMaximum(2.);
  h_xs_ratio->SetNdivisions(505,"XY");
  h_xs_ratio->SetXTitle("P_{T} (GeV)");
  h_xs_ratio->SetYTitle("Data / MC");
  h_xs_ratio->GetXaxis()->SetLabelSize(0.12);
  h_xs_ratio->GetYaxis()->SetLabelSize(0.1);
  h_xs_ratio->GetXaxis()->SetTitleSize(0.12);
  h_xs_ratio->GetYaxis()->SetTitleSize(0.12);
  h_xs_ratio->GetXaxis()->SetTitleOffset(1.);
  h_xs_ratio->GetYaxis()->SetTitleOffset(0.4);

  h_xs_ratio->GetXaxis()->SetRangeUser(150,1000);
  h_xs_ratio->GetYaxis()->SetRangeUser(0.55, 1.45);
  if(ebeebin==1) h_xs_ratio->GetYaxis()->SetRangeUser(0., 2.2);
  if(jetbin==1) h_xs_ratio->GetYaxis()->SetRangeUser(0., 2.2);

  //h_xs_ratio->GetYaxis()->SetNdivisions(102,"Y");
  h_xs_ratio->Draw("pe");
  // tgrs->SetMarkerStyle(8);
  // tgrs->Draw("pe same");
  
}
