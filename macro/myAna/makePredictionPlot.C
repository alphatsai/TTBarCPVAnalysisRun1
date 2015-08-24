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
#include <TPolyLine.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
void makePredictionPlot()
{
    const int NOBS=2; // O2=0, O7=1
    const int NCH=3;  // Combined=2, Electron=0, Muon=1
    const int CH_Electron=0;
    const int CH_Muon=1;
    const int CH_Combined=2;
    const int OBS_O2=0;
    const int OBS_O7=1;
    const int UP=0;
    const int LOW=1;

    const int OBS=OBS_O7;
    const int CH=CH_Combined;

    TFile* fin = new TFile("Final_PredictionTree.root");
    TTree* tin =(TTree*)fin->Get("tree");

    int NPOINTS = tin->GetEntries();
    float predictNonZeroACPMean;
    float predictNonZeroACPUncs[NCH];
    float acpMean[NOBS][NCH]; 
    float acpUncs[NOBS][NCH];
    tin->SetBranchAddress("predictNonZeroACPMean", &predictNonZeroACPMean  );
    tin->SetBranchAddress("predictNonZeroACPUncs", &predictNonZeroACPUncs[0]);
    tin->SetBranchAddress("acpMean", &acpMean[0][0]);
    tin->SetBranchAddress("acpUncs", &acpUncs[0][0]);

    float nonZeroACP[NPOINTS];
    float nonZeroACPUnc[NPOINTS][2];
    float acp[NPOINTS]; 
    float unc_1s[NPOINTS][2];
    float unc_2s[NPOINTS][2];

    for(int idx=0; idx<NPOINTS; idx++){
        tin->GetEntry(idx);

        nonZeroACP[idx]    = predictNonZeroACPMean*100;
        nonZeroACPUnc[idx][UP]  = nonZeroACP[idx] + predictNonZeroACPUncs[CH]*100;
        nonZeroACPUnc[idx][LOW] = nonZeroACP[idx] - predictNonZeroACPUncs[CH]*100;

        acp[idx] = acpMean[OBS][CH]*100;

        unc_1s[idx][UP]  = acp[idx] + acpUncs[OBS][CH]*100;
        unc_2s[idx][UP]  = acp[idx] + acpUncs[OBS][CH]*100*2;
        unc_1s[idx][LOW] = acp[idx] - acpUncs[OBS][CH]*100;
        unc_2s[idx][LOW] = acp[idx] - acpUncs[OBS][CH]*100*2;
        printf("%d\n", idx);
        printf("+ 1s %.2f, 2s %.2f\n", unc_1s[idx][UP], unc_2s[idx][UP]);
        printf("- 1s %.2f, 2s %.2f\n", unc_1s[idx][LOW], unc_2s[idx][LOW]);
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH2F *frame = new TH2F("frame","frame", NPOINTS, nonZeroACP[0], nonZeroACP[NPOINTS-1], NPOINTS, nonZeroACP[0], nonZeroACP[NPOINTS-1]);

    frame->SetStats(kFALSE);
    frame->SetTitle("O_{2} Electron+Muon channel");
    //frame->SetTitle("O_{7} Electron+Muon channel");
    frame->SetXTitle("Ideal non-Zero ACP [%]");
    frame->SetYTitle("Predicted non-Zero ACP [%]");
    frame->Draw();

    //c1->SetLogy();
    c1->SetGridx();
    c1->SetGridy();

    TPolyLine *pl_2s = new TPolyLine(NPOINTS*2);
    TPolyLine *pl_1s = new TPolyLine(NPOINTS*2);
    TPolyLine *o_1s  = new TPolyLine(NPOINTS*2);

    TGraph *pl_med = new TGraph(NPOINTS);
    TGraph *pl_obs = new TGraph(NPOINTS);

    for(int i=0; i<NPOINTS; i++) {
        pl_2s->SetNextPoint( nonZeroACP[i], unc_2s[i][UP] );
        pl_1s->SetNextPoint( nonZeroACP[i], unc_1s[i][UP] );
        o_1s->SetNextPoint(  nonZeroACP[i], nonZeroACPUnc[i][UP]);

        pl_med->SetPoint( i, nonZeroACP[i], nonZeroACP[i]);
        pl_obs->SetPoint( i, nonZeroACP[i], acp[i]       );

    }
    for(int i=NPOINTS-1; i>=0; i--) {
        pl_2s->SetNextPoint(  nonZeroACP[i], unc_2s[i][LOW] );
        pl_1s->SetNextPoint(  nonZeroACP[i], unc_1s[i][LOW] );
        o_1s->SetNextPoint(   nonZeroACP[i], nonZeroACPUnc[i][LOW]);
    }

    pl_2s->SetFillColor(46);
    pl_2s->Draw("f");
    pl_1s->SetFillColor(kGreen);
    pl_1s->Draw("f");
    o_1s->SetFillStyle(3244);
    o_1s->SetFillColor(13);
    o_1s->Draw("f");
    pl_med->SetLineStyle(7);
    pl_med->SetLineWidth(2);
    pl_med->Draw("l");
    pl_obs->Draw("*lsame");

}
