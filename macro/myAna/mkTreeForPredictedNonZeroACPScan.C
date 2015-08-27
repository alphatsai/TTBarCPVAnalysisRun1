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
const std::string  out_s ="./";
const std::string fout_s = out_s + "Final_PredictionTree.root";
const std::string tout_s ="tree";

const int NEXP=1000;
const int NENTRY=13;
const int NOBS=2; // O2=0, O7=1
const int NCH=3;  // Combined=2, Electron=0, Muon=1
const int CH_Electron=0;
const int CH_Muon=1;
const int CH_Combined=2;
const int OBS_O2=0;
const int OBS_O7=1;

const std::string obs_s[NOBS] = {"O2", "O7"};
const std::string  ch_s[NCH]  = {"_El", "_Mu", "" };

float getValue( float per, float evnts, bool isPos);
float getACPMean( float Op, float Om );
float getACPUncs( float Op, float Om );
float getACPUncs( float Op, float Om, float Ope, float Ome );

void mkTreeForPredictedNonZeroACPScanViaSudoExp()
{
    float  assumedACP[NENTRY]  ={  -30,   -25,   -20,   -15,   -10,   -5,   0,    5,    10,    15,    20,    25,    30 };
    string assumedACP_s[NENTRY]={ "m30", "m25", "m20", "m15", "m10", "m5", "0", "p5", "p10", "p15", "p20", "p25", "p30"};

    TFile* fin = new TFile(fin_s.c_str());  

    //TH1D *h_sudoAllAsym[NOBS][NCH], *h_sudoSigAsym[NOBS][NCH]; 
    TH1D *h_mcCutFlow[2]; 
    TH1D *h_sigCutFlow[2]; 
    TH1D *h_SigAsym[NOBS][NCH];
    TH1D *h_bkg[NOBS][NCH];
    TH1D *h_bkgAsym[NOBS][NCH];

    h_mcCutFlow[CH_Muon]           = (TH1D*)fin->Get("MC__Evt_CutFlow_Mu");
    h_mcCutFlow[CH_Electron]       = (TH1D*)fin->Get("MC__Evt_CutFlow_El");
    h_sigCutFlow[CH_Muon]          = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_CutFlow_Mu");
    h_sigCutFlow[CH_Electron]      = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_CutFlow_El");
    h_SigAsym[OBS_O2][CH_Electron] = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O2Asym_El");
    h_SigAsym[OBS_O2][CH_Muon]     = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O2Asym_Mu");
    h_SigAsym[OBS_O2][CH_Combined] = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O2Asym"   );
    h_SigAsym[OBS_O7][CH_Electron] = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O7Asym_El");
    h_SigAsym[OBS_O7][CH_Muon]     = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O7Asym_Mu");
    h_SigAsym[OBS_O7][CH_Combined] = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_O7Asym"   );

    h_bkgAsym[OBS_O2][CH_Electron] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym_El");
    h_bkgAsym[OBS_O2][CH_Muon]     = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym_Mu");
    h_bkgAsym[OBS_O2][CH_Combined] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym"   );
    h_bkgAsym[OBS_O7][CH_Electron] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym_El");
    h_bkgAsym[OBS_O7][CH_Muon]     = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym_Mu");
    h_bkgAsym[OBS_O7][CH_Combined] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym"   );

    TFile* fout = new TFile(fout_s.c_str(), "RECREATE");
    TTree* tout = new TTree(tout_s.c_str(), "");

    TH1D *h_template[NENTRY][NOBS][NCH];
    TH1D *h_sudoExpACPmean[NENTRY][NOBS][NCH];
    TH1D *h_sudoExpACPuncs[NENTRY][NOBS][NCH];
    TH1D *h_sudoSigACPmean[NENTRY][NOBS][NCH];
    TH1D *h_sudoSigACPuncs[NENTRY][NOBS][NCH];

    int lastBin = h_sigCutFlow[CH_Electron]->GetXaxis()->GetLast();
    float nSig[NCH], eSig[NCH], nMC[NCH];
    nSig[CH_Muon]     = h_sigCutFlow[CH_Muon]    ->GetBinContent(lastBin);
    nSig[CH_Electron] = h_sigCutFlow[CH_Electron]->GetBinContent(lastBin);
    nSig[CH_Combined] = nSig[CH_Muon]+nSig[CH_Electron];
    eSig[CH_Muon]     = h_sigCutFlow[CH_Muon]    ->GetBinError(lastBin);
    eSig[CH_Electron] = h_sigCutFlow[CH_Electron]->GetBinError(lastBin);
    eSig[CH_Combined] = sqrt( eSig[CH_Muon]*eSig[CH_Muon] + eSig[CH_Electron]*eSig[CH_Electron] );

    nMC[CH_Muon]     = h_mcCutFlow[CH_Muon]->GetBinContent(lastBin);
    nMC[CH_Electron] = h_mcCutFlow[CH_Electron]->GetBinContent(lastBin);
    nMC[CH_Combined] = nMC[CH_Muon] + nMC[CH_Electron];

    float assumedACPMean;
    float assumedACPUncs[NOBS][NCH];
    float nSigP[NOBS][NCH];
    float nSigM[NOBS][NCH];
    float eSigP[NOBS][NCH];
    float eSigM[NOBS][NCH];
    float nBkgP[NOBS][NCH];
    float eBkgP[NOBS][NCH];
    float nBkgM[NOBS][NCH];
    float eBkgM[NOBS][NCH];
    float MC_nEventsPredictedObsP[NOBS][NCH]; 
    float MC_nEventsPredictedObsM[NOBS][NCH]; 
    float MC_eEventsPredictedObsP[NOBS][NCH]; 
    float MC_eEventsPredictedObsM[NOBS][NCH];
    float MC_acpMean[NOBS][NCH]; 
    float MC_acpUncs[NOBS][NCH]; 
    float SudoExp_acpMean[NOBS][NCH]; 
    float SudoExp_acpUncs[NOBS][NCH]; 
    float SudoSig_acpMean[NOBS][NCH]; 
    float SudoSig_acpUncs[NOBS][NCH]; 

    tout->Branch("assumedACPMean",               &assumedACPMean,                     "assumedACPMean/F"                    );
    tout->Branch("assumedACPUncs",               &assumedACPUncs[0][0],               "assumedACPUncs[2][3]/F"              );
    tout->Branch("nSigP",                        &nSigP[0][0],                        "nSigP[2][3]/F"                       );
    tout->Branch("eSigP",                        &eSigP[0][0],                        "eSigP[2][3]/F"                       );
    tout->Branch("nSigM",                        &nSigM[0][0],                        "nSigM[2][3]/F"                       );
    tout->Branch("eSigM",                        &eSigM[0][0],                        "eSigM[2][3]/F"                       );
    tout->Branch("nBkgP",                        &nBkgP[0][0],                        "nBkgP[2][3]/F"                       );
    tout->Branch("eBkgP",                        &eBkgP[0][0],                        "eBkgP[2][3]/F"                       );
    tout->Branch("nBkgM",                        &nBkgM[0][0],                        "nBkgM[2][3]/F"                       );
    tout->Branch("eBkgM",                        &eBkgM[0][0],                        "eBkgM[2][3]/F"                       );
    tout->Branch("MC_nEventsPredictedObsP",      &MC_nEventsPredictedObsP[0][0],      "MC_nEventsPredictedObsP[2][3]/F"     );
    tout->Branch("MC_nEventsPredictedObsM",      &MC_nEventsPredictedObsM[0][0],      "MC_nEventsPredictedObsM[2][3]/F"     );
    tout->Branch("MC_eEventsPredictedObsP",      &MC_eEventsPredictedObsP[0][0],      "MC_eEventsPredictedObsP[2][3]/F"     );
    tout->Branch("MC_eEventsPredictedObsM",      &MC_eEventsPredictedObsM[0][0],      "MC_eEventsPredictedObsM[2][3]/F"     );
    tout->Branch("MC_acpMean",                   &MC_acpMean[0][0],                   "MC_acpMean[2][3]/F"                  );
    tout->Branch("MC_acpUncs",                   &MC_acpUncs[0][0],                   "MC_acpUncs[2][3]/F"                  );
    tout->Branch("SudoExp_acpMean",              &SudoExp_acpMean[0][0],              "SudoExp_acpMean[2][3]/F"             );
    tout->Branch("SudoExp_acpUncs",              &SudoExp_acpUncs[0][0],              "SudoExp_acpUncs[2][3]/F"             );
    tout->Branch("SudoSig_acpMean",              &SudoSig_acpMean[0][0],              "SudoSig_acpMean[2][3]/F"             );
    tout->Branch("SudoSig_acpUncs",              &SudoSig_acpUncs[0][0],              "SudoSig_acpUncs[2][3]/F"             );

    printf("=============================================================\n");
    printf("Total MC Combined %6.2f, Electron %6.2f, Muon %6.2f\n", nMC[CH_Combined], nMC[CH_Electron], nMC[CH_Muon]);
    for( int entry=0; entry<NENTRY; entry++)
    {
        assumedACPMean = assumedACP[entry]/100;
        // Checking bkg MC effect 
        for( int ich=0; ich<NCH; ich++)
        {
            for( int iobs=0; iobs<NOBS; iobs++)
            {
                h_sudoExpACPmean[entry][iobs][ich] = new TH1D(("SudoExp_"+assumedACP_s[entry]+"ACPmean_"+obs_s[iobs]+ch_s[ich]).c_str(), "", 2000, -100, 100);
                h_sudoExpACPuncs[entry][iobs][ich] = new TH1D(("SudoExp_"+assumedACP_s[entry]+"ACPunus_"+obs_s[iobs]+ch_s[ich]).c_str(), "", 1000,    0,   1);
                h_sudoSigACPmean[entry][iobs][ich] = new TH1D(("SudoSig_"+assumedACP_s[entry]+"ACPmean_"+obs_s[iobs]+ch_s[ich]).c_str(), "", 2000, -100, 100);
                h_sudoSigACPuncs[entry][iobs][ich] = new TH1D(("SudoSig_"+assumedACP_s[entry]+"ACPunus_"+obs_s[iobs]+ch_s[ich]).c_str(), "", 1000,    0,   1);
                h_sudoExpACPmean[entry][iobs][ich]->Sumw2();
                h_sudoExpACPuncs[entry][iobs][ich]->Sumw2();
                h_sudoSigACPmean[entry][iobs][ich]->Sumw2();
                h_sudoSigACPuncs[entry][iobs][ich]->Sumw2();

                h_template[entry][iobs][ich] = new TH1D(("MC_"+assumedACP_s[entry]+"ACP_"+obs_s[iobs]+"Asym"+ch_s[ich]).c_str(), "", 2, 0, 2);
                h_template[entry][iobs][ich] ->Sumw2();

                // Get signal MC or give assumed signal MC (non-0 acp)
                if( assumedACPMean == 0. )
                {
                    nSigP[iobs][ich] = h_SigAsym[iobs][ich]->GetBinContent(2);
                    nSigM[iobs][ich] = h_SigAsym[iobs][ich]->GetBinContent(1);
                }else{
                    nSigP[iobs][ich] = getValue( assumedACPMean, nSig[ich], true  );
                    nSigM[iobs][ich] = getValue( assumedACPMean, nSig[ich], false );
                }
                eSigP[iobs][ich] = sqrt(nSigP[iobs][ich]);
                eSigM[iobs][ich] = sqrt(nSigM[iobs][ich]);
                assumedACPUncs[iobs][ich] = getACPUncs( nSigP[iobs][ich], nSigM[iobs][ich] );
    
                // Get background yieds
                nBkgP[iobs][ich] = h_bkgAsym[iobs][ich]->GetBinContent(2);
                nBkgM[iobs][ich] = h_bkgAsym[iobs][ich]->GetBinContent(1);
                eBkgP[iobs][ich] = h_bkgAsym[iobs][ich]->GetBinError(2);
                eBkgM[iobs][ich] = h_bkgAsym[iobs][ich]->GetBinError(1);

                // Get assumed signal + background, ACP, and fill template  
                MC_nEventsPredictedObsP[iobs][ich] = nBkgP[iobs][ich] + nSigP[iobs][ich];
                MC_nEventsPredictedObsM[iobs][ich] = nBkgM[iobs][ich] + nSigM[iobs][ich];
                MC_eEventsPredictedObsP[iobs][ich] = sqrt( eBkgP[iobs][ich]*eBkgP[iobs][ich] + eSigP[iobs][ich]*eSigP[iobs][ich] );
                MC_eEventsPredictedObsM[iobs][ich] = sqrt( eBkgM[iobs][ich]*eBkgM[iobs][ich] + eSigM[iobs][ich]*eSigM[iobs][ich] );
        
                MC_acpMean[iobs][ich] = getACPMean( MC_nEventsPredictedObsP[iobs][ich], MC_nEventsPredictedObsM[iobs][ich]);
                MC_acpUncs[iobs][ich] = getACPUncs( MC_nEventsPredictedObsP[iobs][ich], MC_nEventsPredictedObsM[iobs][ich], MC_eEventsPredictedObsP[iobs][ich], MC_eEventsPredictedObsM[iobs][ich] );

                h_template[entry][iobs][ich]->SetBinContent(1, MC_nEventsPredictedObsM[iobs][ich]);
                h_template[entry][iobs][ich]->SetBinContent(2, MC_nEventsPredictedObsP[iobs][ich]);
                h_template[entry][iobs][ich]->SetBinError  (1, MC_eEventsPredictedObsM[iobs][ich]);
                h_template[entry][iobs][ich]->SetBinError  (2, MC_eEventsPredictedObsP[iobs][ich]);
            }
        }
            
        // Sudo expriment
        for( int iexp=0; iexp<NEXP; iexp++)
        {
            for( int iobs=0; iobs<NOBS; iobs++ )
            {
                for( int ich=0; ich<NCH; ich++)
                {
                    float nEvtsSudoExp_M, nEvtsSudoSig_M;
                    float nEvtsSudoExp_P, nEvtsSudoSig_P;
                    float eEvtsSudoExp_M, eEvtsSudoSig_M;
                    float eEvtsSudoExp_P, eEvtsSudoSig_P;

                    // Do toy MC
                    TH1D *h_sudoExp = new TH1D("sudoExp_Asym", "", 2, 0, 2); 
                    h_sudoExp->Sumw2();
                    h_sudoExp->FillRandom( h_template[entry][iobs][ich], int(nMC[ich]));
                    nEvtsSudoExp_M = h_sudoExp->GetBinContent(1);
                    nEvtsSudoExp_P = h_sudoExp->GetBinContent(2);
                    eEvtsSudoExp_M = h_sudoExp->GetBinError(1);
                    eEvtsSudoExp_P = h_sudoExp->GetBinError(2);

                    // Subtract bkg, get toy signal 
                    TH1D* h_sudoSig = (TH1D*)h_sudoExp->Clone("sudoSig_Asym");
                    h_sudoSig->Add( h_bkgAsym[iobs][ich], -1);
                    nEvtsSudoSig_M = h_sudoSig->GetBinContent(1);
                    nEvtsSudoSig_P = h_sudoSig->GetBinContent(2);
                    eEvtsSudoSig_M = h_sudoSig->GetBinError(1);
                    eEvtsSudoSig_P = h_sudoSig->GetBinError(2);

                    // Fill pull
                    h_sudoExpACPmean[entry][iobs][ich]->Fill( getACPMean(nEvtsSudoExp_P, nEvtsSudoExp_M)*100. );
                    h_sudoSigACPmean[entry][iobs][ich]->Fill( getACPMean(nEvtsSudoSig_P, nEvtsSudoSig_M)*100. );
                    h_sudoExpACPuncs[entry][iobs][ich]->Fill( getACPUncs(nEvtsSudoExp_P, nEvtsSudoExp_M)*100. );
                    h_sudoSigACPuncs[entry][iobs][ich]->Fill( getACPUncs(nEvtsSudoSig_P, nEvtsSudoSig_M, eEvtsSudoSig_P, eEvtsSudoSig_M )*100. );
        
                    delete h_sudoExp;
                    delete h_sudoSig;
                }
            }
        }
        for( int ich=0; ich<NCH; ich++)
        {
            for( int iobs=0; iobs<NOBS; iobs++ )
            {
                TF1* gaus = new TF1("gaus1_","gaus", assumedACP[entry]-10, assumedACP[entry]+10);
                gaus->SetLineColor(2);

                h_sudoExpACPmean[entry][iobs][ich]->Fit( gaus, "WR" );
                SudoExp_acpMean[iobs][ich] = gaus->GetParameter(1)/100;
                SudoExp_acpUncs[iobs][ich] = gaus->GetParameter(2)/100;

                h_sudoSigACPmean[entry][iobs][ich]->Fit( gaus, "WR" );
                SudoSig_acpMean[iobs][ich] = gaus->GetParameter(1)/100;
                SudoSig_acpUncs[iobs][ich] = gaus->GetParameter(2)/100;

                delete gaus;
            }
        }
        // Debug
        printf("=============================================================\n");
        printf("Assumed ACP %.0f%s\n", assumedACPMean*100, "%");
        printf("O2 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", assumedACPUncs[OBS_O2][CH_Combined]*100, assumedACPUncs[OBS_O2][CH_Electron]*100, assumedACPUncs[OBS_O2][CH_Muon]*100);
        printf("O7 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", assumedACPUncs[OBS_O7][CH_Combined]*100, assumedACPUncs[OBS_O7][CH_Electron]*100, assumedACPUncs[OBS_O7][CH_Muon]*100);
        printf("O2+ Evt Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigP[OBS_O2][CH_Combined], nSigP[OBS_O2][CH_Electron], nSigP[OBS_O2][CH_Muon]);
        printf("O2- Evt Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigM[OBS_O2][CH_Combined], nSigM[OBS_O2][CH_Electron], nSigM[OBS_O2][CH_Muon]);
        printf("O7+ Evt Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigP[OBS_O7][CH_Combined], nSigP[OBS_O7][CH_Electron], nSigP[OBS_O7][CH_Muon]);
        printf("O7- Evt Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigM[OBS_O7][CH_Combined], nSigM[OBS_O7][CH_Electron], nSigM[OBS_O7][CH_Muon]);
        printf("-------------------------------------------------------------\n");
        printf("O2+ Bkg Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nBkgP[OBS_O2][CH_Combined], nBkgP[OBS_O2][CH_Electron], nBkgP[OBS_O2][CH_Muon]);
        printf("O2+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eBkgP[OBS_O2][CH_Combined], eBkgP[OBS_O2][CH_Electron], eBkgP[OBS_O2][CH_Muon]);
        printf("O2- Bkg Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nBkgM[OBS_O2][CH_Combined], nBkgM[OBS_O2][CH_Electron], nBkgM[OBS_O2][CH_Muon]);
        printf("O2- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eBkgM[OBS_O2][CH_Combined], eBkgM[OBS_O2][CH_Electron], eBkgM[OBS_O2][CH_Muon]);
        printf("O7+ Bkg Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nBkgP[OBS_O7][CH_Combined], nBkgP[OBS_O7][CH_Electron], nBkgP[OBS_O7][CH_Muon]);
        printf("O7+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eBkgP[OBS_O7][CH_Combined], eBkgP[OBS_O7][CH_Electron], eBkgP[OBS_O7][CH_Muon]);
        printf("O7- Bkg Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nBkgM[OBS_O7][CH_Combined], nBkgM[OBS_O7][CH_Electron], nBkgM[OBS_O7][CH_Muon]);
        printf("O7- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eBkgM[OBS_O7][CH_Combined], eBkgM[OBS_O7][CH_Electron], eBkgM[OBS_O7][CH_Muon]);
        printf("-------------------------------------------------------------\n");
        printf("O2+ All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_nEventsPredictedObsP[OBS_O2][CH_Combined], MC_nEventsPredictedObsP[OBS_O2][CH_Electron], MC_nEventsPredictedObsP[OBS_O2][CH_Muon]);
        printf("O2+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_eEventsPredictedObsP[OBS_O2][CH_Combined], MC_eEventsPredictedObsP[OBS_O2][CH_Electron], MC_eEventsPredictedObsP[OBS_O2][CH_Muon]);
        printf("O2- All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_nEventsPredictedObsM[OBS_O2][CH_Combined], MC_nEventsPredictedObsM[OBS_O2][CH_Electron], MC_nEventsPredictedObsM[OBS_O2][CH_Muon]);
        printf("O2- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_eEventsPredictedObsM[OBS_O2][CH_Combined], MC_eEventsPredictedObsM[OBS_O2][CH_Electron], MC_eEventsPredictedObsM[OBS_O2][CH_Muon]);
        printf("O7+ All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_nEventsPredictedObsP[OBS_O7][CH_Combined], MC_nEventsPredictedObsP[OBS_O7][CH_Electron], MC_nEventsPredictedObsP[OBS_O7][CH_Muon]);
        printf("O7+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_eEventsPredictedObsP[OBS_O7][CH_Combined], MC_eEventsPredictedObsP[OBS_O7][CH_Electron], MC_eEventsPredictedObsP[OBS_O7][CH_Muon]);
        printf("O7- All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_nEventsPredictedObsM[OBS_O7][CH_Combined], MC_nEventsPredictedObsM[OBS_O7][CH_Electron], MC_nEventsPredictedObsM[OBS_O7][CH_Muon]);
        printf("O7- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", MC_eEventsPredictedObsM[OBS_O7][CH_Combined], MC_eEventsPredictedObsM[OBS_O7][CH_Electron], MC_eEventsPredictedObsM[OBS_O7][CH_Muon]);
        printf("-------------------------------------------------------------\n");
        printf("O2 Mean Combined %6.2f, Electron %6.2f, Muon %6.2f\n", MC_acpMean[OBS_O2][CH_Combined]*100, MC_acpMean[OBS_O2][CH_Electron]*100, MC_acpMean[OBS_O2][CH_Muon]*100);
        printf("O2 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", MC_acpUncs[OBS_O2][CH_Combined]*100, MC_acpUncs[OBS_O2][CH_Electron]*100, MC_acpUncs[OBS_O2][CH_Muon]*100);
        printf("O7 Mean Combined %6.2f, Electron %6.2f, Muon %6.2f\n", MC_acpMean[OBS_O7][CH_Combined]*100, MC_acpMean[OBS_O7][CH_Electron]*100, MC_acpMean[OBS_O7][CH_Muon]*100);
        printf("O7 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", MC_acpUncs[OBS_O7][CH_Combined]*100, MC_acpUncs[OBS_O7][CH_Electron]*100, MC_acpUncs[OBS_O7][CH_Muon]*100);
        printf("\n");

        tout->Fill();
    }//Entry End

    fout->Write(); 
    fout->Close(); 
    cout<<">> Wrote to file: "<<fout_s.c_str()<<endl;
}

float getValue( float per, float evnts, bool isPos)
{
    float x = (1-per)*evnts/2;
    float out; 
    if( isPos ) out=evnts-x; 
    else out=x;
    return out; 
}
float getACPMean( float Op, float Om )
{
    return (Op-Om)/(Op+Om); 
}
float getACPUncs( float Op, float Om )
{
    float acpUnc = sqrt( 4*Op*Om/pow(Op+Om,3.) );
    return acpUnc; 
}
float getACPUncs( float Op, float Om, float Ope, float Ome )
{
    float acpUnc = sqrt( 4*(Op*Op*Ome*Ome+Om*Om*Ope*Ope)/pow(Op+Om,4.) );
    return acpUnc;
} 
