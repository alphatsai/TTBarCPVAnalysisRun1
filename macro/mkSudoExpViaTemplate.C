#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TLine.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
const std::string fin_s  ="Final_histograms_SemiLeptanic.root";
const std::string fout_s ="Final_SudoExpViaTemplate.root";
const std::string tout_s ="tree";

const int NENTRY=12;
const int NOBS=2; // O2=0, O7=1
const int NCH=3;  // Combined=2, Electron=0, Muon=1
const int CH_Electron=0;
const int CH_Muon=1;
const int CH_Combined=2;
const int OBS_O2=0;
const int OBS_O7=1;

void mkSudoExpViaTemplate()
{
    TFile* fin = new TFile(fin_s.c_str());  

    TH1D *h_all[NOBS][NCH], *h_sudoAll[NOBS][NCH], *h_sudoAllAsym[NOBS][NCH], *h_sudoSig[NOBS][NCH], *h_sudoSigAsym[NOBS][NCH]; 
    TH1D *h_bkg[NOBS][NCH];
    TH1D *h_bkgAsym[NOBS][NCH];

    h_all[OBS_O2][CH_Muon]     = (TH1D*)fin->Get("MC__Evt_O2_Mu");
    h_all[OBS_O2][CH_Electron] = (TH1D*)fin->Get("MC__Evt_O2_El");
    h_all[OBS_O2][CH_Combined] = (TH1D*)fin->Get("MC__Evt_O2");
    h_all[OBS_O7][CH_Muon]     = (TH1D*)fin->Get("MC__Evt_O7_Mu");
    h_all[OBS_O7][CH_Electron] = (TH1D*)fin->Get("MC__Evt_O7_El");
    h_all[OBS_O7][CH_Combined] = (TH1D*)fin->Get("MC__Evt_O7");

    h_bkg[OBS_O2][CH_Electron]     = (TH1D*)fin->Get("BkgMC__Evt_O2_El"    );
    h_bkg[OBS_O2][CH_Muon]         = (TH1D*)fin->Get("BkgMC__Evt_O2_Mu"    );
    h_bkg[OBS_O2][CH_Combined]     = (TH1D*)fin->Get("BkgMC__Evt_O2"       );
    h_bkg[OBS_O7][CH_Electron]     = (TH1D*)fin->Get("BkgMC__Evt_O7_El"    );
    h_bkg[OBS_O7][CH_Muon]         = (TH1D*)fin->Get("BkgMC__Evt_O7_Mu"    );
    h_bkg[OBS_O7][CH_Combined]     = (TH1D*)fin->Get("BkgMC__Evt_O7"       );
    h_bkgAsym[OBS_O2][CH_Electron] = (TH1D*)fin->Get("BkgMC__Evt_O2Asym_El");
    h_bkgAsym[OBS_O2][CH_Muon]     = (TH1D*)fin->Get("BkgMC__Evt_O2Asym_Mu");
    h_bkgAsym[OBS_O2][CH_Combined] = (TH1D*)fin->Get("BkgMC__Evt_O2Asym"   );
    h_bkgAsym[OBS_O7][CH_Electron] = (TH1D*)fin->Get("BkgMC__Evt_O7Asym_El");
    h_bkgAsym[OBS_O7][CH_Muon]     = (TH1D*)fin->Get("BkgMC__Evt_O7Asym_Mu");
    h_bkgAsym[OBS_O7][CH_Combined] = (TH1D*)fin->Get("BkgMC__Evt_O7Asym"   );

    TFile* fout = new TFile(fout_s.c_str(), "RECREATE");

    h_sudoAll[OBS_O2][CH_Electron]     = new TH1D("SudoExp_O2_El",        "", 100, -5, 5);
    h_sudoAll[OBS_O2][CH_Muon]         = new TH1D("SudoExp_O2_Mu",        "", 100, -5, 5);
    h_sudoAll[OBS_O2][CH_Combined]     = new TH1D("SudoExp_O2",           "", 100, -5, 5);
    h_sudoAll[OBS_O7][CH_Electron]     = new TH1D("SudoExp_O7_El",        "", 100, -5, 5);
    h_sudoAll[OBS_O7][CH_Muon]         = new TH1D("SudoExp_O7_Mu",        "", 100, -5, 5);
    h_sudoAll[OBS_O7][CH_Combined]     = new TH1D("SudoExp_O7",           "", 100, -5, 5);
    h_sudoAllAsym[OBS_O2][CH_Electron] = new TH1D("SudoExp_O2Asym_El",    "",   2,  0, 2);
    h_sudoAllAsym[OBS_O2][CH_Muon]     = new TH1D("SudoExp_O2Asym_Mu",    "",   2,  0, 2);
    h_sudoAllAsym[OBS_O2][CH_Combined] = new TH1D("SudoExp_O2Asym",       "",   2,  0, 2);
    h_sudoAllAsym[OBS_O7][CH_Electron] = new TH1D("SudoExp_O7Asym_El",    "",   2,  0, 2);
    h_sudoAllAsym[OBS_O7][CH_Muon]     = new TH1D("SudoExp_O7Asym_Mu",    "",   2,  0, 2);
    h_sudoAllAsym[OBS_O7][CH_Combined] = new TH1D("SudoExp_O7Asym",       "",   2,  0, 2);
    h_sudoSig[OBS_O2][CH_Electron]     = new TH1D("SudoExpSig_O2_El",     "", 100, -5, 5);
    h_sudoSig[OBS_O2][CH_Muon]         = new TH1D("SudoExpSig_O2_Mu",     "", 100, -5, 5);
    h_sudoSig[OBS_O2][CH_Combined]     = new TH1D("SudoExpSig_O2",        "", 100, -5, 5);
    h_sudoSig[OBS_O7][CH_Electron]     = new TH1D("SudoExpSig_O7_El",     "", 100, -5, 5);
    h_sudoSig[OBS_O7][CH_Muon]         = new TH1D("SudoExpSig_O7_Mu",     "", 100, -5, 5);
    h_sudoSig[OBS_O7][CH_Combined]     = new TH1D("SudoExpSig_O7",        "", 100, -5, 5);
    h_sudoSigAsym[OBS_O2][CH_Electron] = new TH1D("SudoExpSig_O2Asym_El", "",   2,  0, 2);
    h_sudoSigAsym[OBS_O2][CH_Muon]     = new TH1D("SudoExpSig_O2Asym_Mu", "",   2,  0, 2);
    h_sudoSigAsym[OBS_O2][CH_Combined] = new TH1D("SudoExpSig_O2Asym",    "",   2,  0, 2);
    h_sudoSigAsym[OBS_O7][CH_Electron] = new TH1D("SudoExpSig_O7Asym_El", "",   2,  0, 2);
    h_sudoSigAsym[OBS_O7][CH_Muon]     = new TH1D("SudoExpSig_O7Asym_Mu", "",   2,  0, 2);
    h_sudoSigAsym[OBS_O7][CH_Combined] = new TH1D("SudoExpSig_O7Asym",    "",   2,  0, 2);

    for( int iobs=0; iobs<NOBS; iobs++)
    {
        h_sudoAll[iobs][CH_Combined]    ->Sumw2();
        h_sudoAllAsym[iobs][CH_Combined]->Sumw2();
        h_sudoSig[iobs][CH_Combined]    ->Sumw2();
        h_sudoSigAsym[iobs][CH_Combined]->Sumw2();
        for( int ich=0; ich<CH_Combined; ich++)
        {
            double evts = h_all[iobs][ich]->Integral();
            h_all[iobs][ich]->Scale((1./evts)); 
 
            h_sudoAll[iobs][ich]    ->Sumw2();
            h_sudoAllAsym[iobs][ich]->Sumw2();
            h_sudoSig[iobs][ich]    ->Sumw2();
            h_sudoSigAsym[iobs][ich]->Sumw2();

            // Get sudo exp.
            // GetRandom used Rndm()
            //for( int evt=0; evt<=int(evts); evt++)
            //{
            //    double rd = h_all[iobs][ich]->GetRandom();
            //    h_sudoSig[iobs][ich]        ->Fill(rd);
            //    h_sudoSig[iobs][CH_Combined]->Fill(rd);
            //    h_sudoAll[iobs][ich]        ->Fill(rd);
            //    h_sudoAll[iobs][CH_Combined]->Fill(rd);
            //    if( rd < 0 ){ 
            //        h_sudoAllAsym[iobs][ich]        ->Fill(0); 
            //        h_sudoAllAsym[iobs][CH_Combined]->Fill(0); 
            //    }else{         
            //        h_sudoAllAsym[iobs][ich]        ->Fill(1); 
            //        h_sudoAllAsym[iobs][CH_Combined]->Fill(1); 
            //    }
            //}
            //
            // FillRandom is Poisson in same binning histo
            h_sudoSig[iobs][ich]->FillRandom( h_all[iobs][ich], evts );
            h_sudoAll[iobs][ich]->FillRandom( h_all[iobs][ich], evts );

            // Subtract bkg MC
            h_sudoSig[iobs][ich]->Add( h_bkg[iobs][ich], -1);

        }// ich end
        h_sudoAll[iobs][CH_Combined]->Add( h_sudoAll[iobs][CH_Electron]);
        h_sudoAll[iobs][CH_Combined]->Add( h_sudoAll[iobs][CH_Muon]    );
        h_sudoSig[iobs][CH_Combined]->Add( h_sudoSig[iobs][CH_Electron]);
        h_sudoSig[iobs][CH_Combined]->Add( h_sudoSig[iobs][CH_Muon]    );
        // Subtract bkg MC
        h_sudoSig[iobs][CH_Combined]->Add( h_bkg[iobs][CH_Electron], -1);
        h_sudoSig[iobs][CH_Combined]->Add( h_bkg[iobs][CH_Muon],     -1);

        // Caculate asym for sudo all
        for( int ich=0; ich<NCH; ich++)
        {
            double e1(0), e2(0), v1(0), v2(0);
            for( int b=1; b<=h_sudoAll[iobs][ich]->GetXaxis()->GetLast(); b++)
            {
                double v = h_sudoAll[iobs][ich]->GetBinContent(b);
                double e = h_sudoAll[iobs][ich]->GetBinError(b);
                if( h_sudoAll[iobs][ich]->GetBinLowEdge(b) < 0 ){ 
                    v1 += v; 
                    e1 = sqrt( e*e + e1*e1 ); 
                }else{
                    v2 += v; 
                    e2 = sqrt( e*e + e2*e2 ); 
                }
                if( v < 0 ) cout<<">> [WARNING] Sudo all has nagtive entry [obs:ch:b:v] "<<iobs<<", "<<ich<<", "<<b<<", "<<v<<endl; 
            }
            h_sudoAllAsym[iobs][ich]->SetBinContent(1, v1);
            h_sudoAllAsym[iobs][ich]->SetBinContent(2, v2);
            h_sudoAllAsym[iobs][ich]->SetBinError(1, e1);
            h_sudoAllAsym[iobs][ich]->SetBinError(2, e2);

        }// ich end

        // Caculate asym for sudo signal
        for( int ich=0; ich<NCH; ich++)
        {
            double e1(0), e2(0), v1(0), v2(0);
            for( int b=1; b<=h_sudoSig[iobs][ich]->GetXaxis()->GetLast(); b++)
            {
                double v = h_sudoSig[iobs][ich]->GetBinContent(b);
                double e = h_sudoSig[iobs][ich]->GetBinError(b);
                //if( v < 0 ) continue;
                if( h_sudoSig[iobs][ich]->GetBinLowEdge(b) < 0 ){ 
                    v1 += v; 
                    e1 = sqrt( e*e + e1*e1 ); 
                }else{
                    v2 += v; 
                    e2 = sqrt( e*e + e2*e2 ); 
                }
                if( v < 0 ) cout<<">> [WARNING] Sudo signal nagtive entry [obs:ch:b:v] "<<iobs<<", "<<ich<<", "<<b<<", "<<v<<endl; 
            }
            h_sudoSigAsym[iobs][ich]->SetBinContent(1, v1);
            h_sudoSigAsym[iobs][ich]->SetBinContent(2, v2);
            h_sudoSigAsym[iobs][ich]->SetBinError(1, e1);
            h_sudoSigAsym[iobs][ich]->SetBinError(2, e2);

        }// ich end
    }// iobs end

    fout->Write(); 
    fout->Close(); 
}
