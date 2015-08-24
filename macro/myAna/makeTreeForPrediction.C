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
const std::string fout_s ="Final_PredictionTree.root";
const std::string tout_s ="tree";

const int NENTRY=12;
const int NOBS=2; // O2=0, O7=1
const int NCH=3;  // Combined=2, Electron=0, Muon=1
const int CH_Electron=0;
const int CH_Muon=1;
const int CH_Combined=2;
const int OBS_O2=0;
const int OBS_O7=1;

float getValue( float per, float evnts, bool isPos);
float getACPMean( float Op, float Om );
float getACPUncs( float Op, float Om );
float getACPUncs( float Op, float Om, float Ope, float Ome );

void makeTreeForPrediction()
{
    float nonZeroACP[NENTRY]={ -30, -25, -20, -15, -10, -5, 5, 10, 15, 20, 25, 30};

    TFile* fin = new TFile(fin_s.c_str());  

    TH1D *h_sig[2], *h_bkg[NOBS][NCH];
    h_sig[CH_Muon]             = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_CutFlow_Mu");
    h_sig[CH_Electron]         = (TH1D*)fin->Get("TTJets_SemiLeptMGDecays__Evt_CutFlow_El");
    h_bkg[OBS_O2][CH_Electron] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym_El");
    h_bkg[OBS_O2][CH_Muon]     = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym_Mu");
    h_bkg[OBS_O2][CH_Combined] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O2Asym"   );
    h_bkg[OBS_O7][CH_Electron] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym_El");
    h_bkg[OBS_O7][CH_Muon]     = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym_Mu");
    h_bkg[OBS_O7][CH_Combined] = (TH1D*)fin->Get("BkgMC_TTJetsNonSemiLeptMGDecaysIncluded__Evt_O7Asym"   );

    int lastBin = h_sig[CH_Electron]->GetXaxis()->GetLast();
    float nSig[NCH], eSig[NCH];
    nSig[CH_Muon]     = h_sig[CH_Muon]    ->GetBinContent(lastBin);
    nSig[CH_Electron] = h_sig[CH_Electron]->GetBinContent(lastBin);
    nSig[CH_Combined] = nSig[CH_Muon]+nSig[CH_Electron];
    eSig[CH_Muon]     = h_sig[CH_Muon]    ->GetBinError(lastBin);
    eSig[CH_Electron] = h_sig[CH_Electron]->GetBinError(lastBin);
    eSig[CH_Combined] = sqrt( eSig[CH_Muon]*eSig[CH_Muon] + eSig[CH_Electron]*eSig[CH_Electron] );

    TFile* fout = new TFile(fout_s.c_str(), "RECREATE");
    TTree* tout = new TTree(tout_s.c_str(), "");

    float predictNonZeroACPMean;
    float predictNonZeroACPUncs[NCH];
    float nSigP[NCH];
    float nSigM[NCH];
    float eSigP[NCH];
    float eSigM[NCH];
    float nBkgP[NOBS][NCH];
    float eBkgP[NOBS][NCH];
    float nBkgM[NOBS][NCH];
    float eBkgM[NOBS][NCH];
    float nEventsPredictedObsP[NOBS][NCH]; 
    float nEventsPredictedObsM[NOBS][NCH]; 
    float eEventsPredictedObsP[NOBS][NCH]; 
    float eEventsPredictedObsM[NOBS][NCH];
    float acpMean[NOBS][NCH]; 
    float acpUncs[NOBS][NCH]; 

    tout->Branch("predictNonZeroACPMean", &predictNonZeroACPMean,      "predictNonZeroACPMean/F");
    tout->Branch("predictNonZeroACPUncs", &predictNonZeroACPUncs[0],   "predictNonZeroACPUncs[3]/F");
    tout->Branch("nSigP",                 &nSigP[0],                   "nSigP[3]/F");
    tout->Branch("eSigP",                 &eSigP[0],                   "eSigP[3]/F");
    tout->Branch("nSigM",                 &nSigM[0],                   "nSigM[3]/F");
    tout->Branch("eSigM",                 &eSigM[0],                   "eSigM[3]/F");
    tout->Branch("nBkgP",                 &nBkgP[0][0],                "nBkgP[2][3]/F");
    tout->Branch("eBkgP",                 &eBkgP[0][0],                "eBkgP[2][3]/F");
    tout->Branch("nBkgM",                 &nBkgM[0][0],                "nBkgM[2][3]/F");
    tout->Branch("eBkgM",                 &eBkgM[0][0],                "eBkgM[2][3]/F");
    tout->Branch("nEventsPredictedObsP",  &nEventsPredictedObsP[0][0], "nEventsPredictedObsP[2][3]/F");
    tout->Branch("nEventsPredictedObsM",  &nEventsPredictedObsM[0][0], "nEventsPredictedObsM[2][3]/F");
    tout->Branch("eEventsPredictedObsP",  &eEventsPredictedObsP[0][0], "eEventsPredictedObsP[2][3]/F");
    tout->Branch("eEventsPredictedObsM",  &eEventsPredictedObsM[0][0], "eEventsPredictedObsM[2][3]/F");
    tout->Branch("acpMean",               &acpMean[0][0],              "acpMean[2][3]/F");
    tout->Branch("acpUncs",               &acpUncs[0][0],              "acpUncs[2][3]/F");

    for( int entry=0; entry<NENTRY; entry++)
    {
        predictNonZeroACPMean = nonZeroACP[entry]/100;
        for( int ich=0; ich<NCH; ich++)
        {
            nSigP[ich] = getValue( predictNonZeroACPMean, nSig[ich], true  );
            nSigM[ich] = getValue( predictNonZeroACPMean, nSig[ich], false );
            eSigP[ich] = sqrt(nSigP[ich]);
            eSigM[ich] = sqrt(nSigM[ich]);
            predictNonZeroACPUncs[ich] = getACPUncs( nSigP[ich], nSigM[ich] );

            for( int iobs=0; iobs<NOBS; iobs++)
            {
                nBkgP[iobs][ich] = h_bkg[iobs][ich]->GetBinContent(2);
                nBkgM[iobs][ich] = h_bkg[iobs][ich]->GetBinContent(1);
                eBkgP[iobs][ich] = h_bkg[iobs][ich]->GetBinError(2);
                eBkgM[iobs][ich] = h_bkg[iobs][ich]->GetBinError(1);

                nEventsPredictedObsP[iobs][ich] = nBkgP[iobs][ich] + nSigP[ich];
                nEventsPredictedObsM[iobs][ich] = nBkgM[iobs][ich] + nSigM[ich];
                eEventsPredictedObsP[iobs][ich] = sqrt( eBkgP[iobs][ich]*eBkgP[iobs][ich] + eSigP[ich]*eSigP[ich] );
                eEventsPredictedObsM[iobs][ich] = sqrt( eBkgM[iobs][ich]*eBkgM[iobs][ich] + eSigM[ich]*eSigM[ich] );

                acpMean[iobs][ich] = getACPMean( nEventsPredictedObsP[iobs][ich], nEventsPredictedObsM[iobs][ich]);
                acpUncs[iobs][ich] = getACPUncs( nEventsPredictedObsP[iobs][ich], nEventsPredictedObsM[iobs][ich], eEventsPredictedObsP[iobs][ich], eEventsPredictedObsM[iobs][ich] );
            }
        }

        // Debug
        printf("=============================================================\n");
        printf("Non ACP %.0f%s\n", predictNonZeroACPMean*100, "%");
        printf("Acp Unc Combined %6.2f, Electron %6.2f, Muon %6.2f\n", predictNonZeroACPUncs[CH_Combined]*100, predictNonZeroACPUncs[CH_Electron]*100, predictNonZeroACPUncs[CH_Muon]*100);
        printf("Event + Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigP[CH_Combined], nSigP[CH_Electron], nSigP[CH_Muon]);
        printf("Event - Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nSigM[CH_Combined], nSigM[CH_Electron], nSigM[CH_Muon]);
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
        printf("O2+ All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nEventsPredictedObsP[OBS_O2][CH_Combined], nEventsPredictedObsP[OBS_O2][CH_Electron], nEventsPredictedObsP[OBS_O2][CH_Muon]);
        printf("O2+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eEventsPredictedObsP[OBS_O2][CH_Combined], eEventsPredictedObsP[OBS_O2][CH_Electron], eEventsPredictedObsP[OBS_O2][CH_Muon]);
        printf("O2- All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nEventsPredictedObsM[OBS_O2][CH_Combined], nEventsPredictedObsM[OBS_O2][CH_Electron], nEventsPredictedObsM[OBS_O2][CH_Muon]);
        printf("O2- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eEventsPredictedObsM[OBS_O2][CH_Combined], eEventsPredictedObsM[OBS_O2][CH_Electron], eEventsPredictedObsM[OBS_O2][CH_Muon]);
        printf("O7+ All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nEventsPredictedObsP[OBS_O7][CH_Combined], nEventsPredictedObsP[OBS_O7][CH_Electron], nEventsPredictedObsP[OBS_O7][CH_Muon]);
        printf("O7+ Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eEventsPredictedObsP[OBS_O7][CH_Combined], eEventsPredictedObsP[OBS_O7][CH_Electron], eEventsPredictedObsP[OBS_O7][CH_Muon]);
        printf("O7- All Combined %6.0f, Electron %6.0f, Muon %6.0f\n", nEventsPredictedObsM[OBS_O7][CH_Combined], nEventsPredictedObsM[OBS_O7][CH_Electron], nEventsPredictedObsM[OBS_O7][CH_Muon]);
        printf("O7- Unc Combined %6.0f, Electron %6.0f, Muon %6.0f\n", eEventsPredictedObsM[OBS_O7][CH_Combined], eEventsPredictedObsM[OBS_O7][CH_Electron], eEventsPredictedObsM[OBS_O7][CH_Muon]);
        printf("-------------------------------------------------------------\n");
        printf("O2 Mean Combined %6.2f, Electron %6.2f, Muon %6.2f\n", acpMean[OBS_O2][CH_Combined]*100, acpMean[OBS_O2][CH_Electron]*100, acpMean[OBS_O2][CH_Muon]*100);
        printf("O2 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", acpUncs[OBS_O2][CH_Combined]*100, acpUncs[OBS_O2][CH_Electron]*100, acpUncs[OBS_O2][CH_Muon]*100);
        printf("O7 Mean Combined %6.2f, Electron %6.2f, Muon %6.2f\n", acpMean[OBS_O7][CH_Combined]*100, acpMean[OBS_O7][CH_Electron]*100, acpMean[OBS_O7][CH_Muon]*100);
        printf("O7 Uncs Combined %6.2f, Electron %6.2f, Muon %6.2f\n", acpUncs[OBS_O7][CH_Combined]*100, acpUncs[OBS_O7][CH_Electron]*100, acpUncs[OBS_O7][CH_Muon]*100);
        printf("\n");

        tout->Fill();
    }//Entry End
    fout->Write(); 
    fout->Close(); 
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
