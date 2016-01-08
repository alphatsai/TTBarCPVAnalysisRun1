#include <string>
#include "help.C"
#define NPAR 2
vector<Double_t> info;
vector<Double_t> info_err;
void mkTemplate( TFile* f, std::string hName, std::string output=".", int rebin=1 )
{
    TFile* fout = new TFile((output+"/Template_"+hName+".root").c_str(), "RECREATE");
    TH1D *h_data, *h_sigMC, *h_bkgMC;
    TH1D *h_data_El, *h_sigMC_El, *h_bkgMC_El;
    TH1D *h_data_Mu, *h_sigMC_Mu, *h_bkgMC_Mu;

    h_data  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA");
    h_sigMC = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_El").c_str()))->Clone("SigMC");
    h_bkgMC = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_El").c_str()))->Clone("BkgMC");

    h_data_El  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA_El");
    h_sigMC_El = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_El").c_str()))->Clone("SigMC_El");
    h_bkgMC_El = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_El").c_str()))->Clone("BkgMC_El");

    h_data_Mu  = (TH1D*)((TH1D*)f->Get(("DATA_Muon__"+hName+"_Mu").c_str()))->Clone("DATA_Mu");
    h_sigMC_Mu = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_Mu").c_str()))->Clone("SigMC_Mu");
    h_bkgMC_Mu = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_Mu").c_str()))->Clone("BkgMC_Mu");

    h_data->Rebin(rebin);  h_data_El->Rebin(rebin);  h_data_Mu->Rebin(rebin);
    h_sigMC->Rebin(rebin); h_sigMC_El->Rebin(rebin); h_sigMC_Mu->Rebin(rebin);
    h_bkgMC->Rebin(rebin); h_bkgMC_El->Rebin(rebin); h_bkgMC_Mu->Rebin(rebin);
    
    fix(h_data);  fix(h_data_El);  fix(h_data_Mu);
    fix(h_sigMC); fix(h_sigMC_El); fix(h_sigMC_Mu);
    fix(h_bkgMC); fix(h_bkgMC_El); fix(h_bkgMC_Mu);

    h_data->Add(h_data_Mu);
    h_sigMC->Add(h_sigMC_Mu);
    h_bkgMC->Add(h_bkgMC_Mu);

    getHistStatUnc( h_data,     "DATA_Stat");
    getHistStatUnc( h_data_El,  "DATA_El_Stat");
    getHistStatUnc( h_data_Mu,  "DATA_Mu_Stat");
    getHistStatUnc( h_sigMC,    "SigMC_Stat");
    getHistStatUnc( h_sigMC_El, "SigMC_El_Stat");
    getHistStatUnc( h_sigMC_Mu, "SigMC_Mu_Stat");
    getHistStatUnc( h_bkgMC,    "BkgMC_Stat");
    getHistStatUnc( h_bkgMC_El, "BkgMC_El_Stat");
    getHistStatUnc( h_bkgMC_Mu, "BkgMC_Mu_Stat");

    fout->Write(); 
}
void getHistStatUnc( TH1D* h, std::string name )
{
    TH1D* h_u = (TH1D*)h->Clone((name+"Up").c_str());
    TH1D* h_d = (TH1D*)h->Clone((name+"Down").c_str());
    int bins=h->GetNbinsX();
    for( int i=1; i<=bins; i++)
    {
        h_u->SetBinContent(i, h->GetBinContent(i)+h->GetBinError(i));
        h_d->SetBinContent(i, h->GetBinContent(i)-h->GetBinError(i));
    }
}
void mkFittedTemplate( TFile* f, std::string output="", std::string output="." )
{
    Ifit(f,)


}

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
Double_t* Ifit( TFile* fin, std::string output=".", int rebin=1, int channel=0, int fit_data=1 )
{
  Double_t* fitted = new Double_t[8];

  TCanvas *c1 = new TCanvas("HF1", "Histos1", 0, 0, 600, 600);
  double count=0;
  dataColl.clear();
  sigColl.clear();
  bkgColl.clear();

  totalColl.clear();
  ctauColl.clear();

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
  c1->Draw();  
  hdata->SetLineColor(1);

  hdata->SetNdivisions(505,"XY");
  hdata->SetXTitle("t\bar{t} mass (GeV)");
  hdata->SetYTitle("Entries");
  hdata->SetTitle("");
  hdata->SetMarkerStyle(8);
  hdata->SetMinimum(0.);
  
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
  
  TLegend *tleg = new TLegend(0.5, 0.65, 0.9, 0.85);
  char text[50];
  sprintf(text,"Top Mass");
  tleg->SetHeader(text);
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  sprintf(text,"Data %5.1f events",hdata->Integral());
  tleg->AddEntry(hdata,text,"pl");
  sprintf(text,"Fitted %5.1f events",hfit->Integral());
  
  tleg->AddEntry(hfit,text,"l");
  sprintf(text,"SIG %5.1f #pm %5.1f events",para[0], errpara[0]);
  tleg->AddEntry(hsig,text,"f");
  sprintf(text,"BKG %5.1f #pm %5.1f events",para[1], errpara[1]);
  
  tleg->AddEntry(hbkg,text,"f");
  tleg->Draw();
  
  TLatex *tlx = new TLatex();
  tlx->SetTextSize(0.035);
  tlx->DrawLatex(-0.5, hdata->GetMaximum(), "CMS Preliminary #sqrt{s} = 8TeV, L=21 fb^{-1}");

  gPad->RedrawAxis();
  
  if(fit_data>0) hdata->Chi2Test(hfit,"P");

  c1->SaveAs(("FittingResults_"+output+ch+".pdf").c_str());
  return fitted;
}
