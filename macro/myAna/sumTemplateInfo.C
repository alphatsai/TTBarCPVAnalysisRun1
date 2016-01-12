#include <TROOT.h>
#include <TFile.h>
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
#include "getTemplateInfo.C"
using namespace std;
void sumTemplateInfo( TFile* f, string output=".", string xElTitle="", string xMuTitel="", string yTitle="Events", bool logy=true )
{
    const int nElSyst=6, nMuSyst=7;
    float *mcEl[nElSyst], *sigEl[nElSyst], *bkgEl[nElSyst];
    float *mcMu[nMuSyst], *sigMu[nMuSyst], *bkgMu[nMuSyst];
    float sumw2El[2], sumw2sigEl[2], sumw2bkgEl[2];
    float sumw2Mu[2], sumw2sigMu[2], sumw2bkgMu[2];
    string systNameEl[nElSyst]={"Stat", "PU", "JER", "BTagSF", "TopPT", "elID"};
    string systNameMu[nMuSyst]={"Stat", "PU", "JER", "BTagSF", "TopPT", "muID", "muISO"};
    FILE* outTxt;

    outTxt = fopen((output+"/SystUncertainties.txt").c_str(),"w");
    fprintf( outTxt, "Electron channel\n", );
    fprintf( outTxt, "%8s", " ");
    for( int i=0; i<nElSyst; i++ )
    {
        mcEl[i]  = drawSyst( f, "MC_El",    systNameEl[i], output, xElTitle, yTitle, logy);
        sigEl[i] = drawSyst( f, "SigMC_El", systNameEl[i], output, xElTitle, yTitle, logy);
        bkgEl[i] = drawSyst( f, "BkgMC_El", systNameEl[i], output, xElTitle, yTitle, logy);
        fprintf( outTxt, "%8s", systNameEl[i].c_str() );
    }
    fprintf( outTxt, "\nMC:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (mcEl[i][1]-mcEl[i][0])/mcEl[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (mcEl[i][2]-mcEl[i][0])/mcEl[i][0]*100 ); }   
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (sigEl[i][1]-sigEl[i][0])/sigEl[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (sigEl[i][2]-sigEl[i][0])/sigEl[i][0]*100 ); }   
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (bkgEl[i][1]-bkgEl[i][0])/bkgEl[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+8.2f", (bkgEl[i][2]-bkgEl[i][0])/bkgEl[i][0]*100 ); }   
    fprintf( outTxt, "\n\n");

    fprintf( outTxt, "Muon channel\n" );
    fprintf( outTxt, "%8s", " ");
    for( int i=0; i<nMuSyst; i++ )
    {
        mcMu[i]  = drawSyst( f, "MC_Mu",    systNameMu[i], output, xMuTitel, yTitle, logy);
        sigMu[i] = drawSyst( f, "SigMC_Mu", systNameMu[i], output, xMuTitel, yTitle, logy);
        bkgMu[i] = drawSyst( f, "BkgMC_Mu", systNameMu[i], output, xMuTitel, yTitle, logy);
        fprintf( outTxt, "%8s", systNameMu[i].c_str() );
    }
    fprintf( outTxt, "%8s", "Syst" );
    fprintf( outTxt, "\nMC:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(mcMu[i][1]-mcMu[i][0])/mcMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2Mu[0]+=v*v; }  
    fprintf( outTxt, "%+8.2f", sqrt(sumw2Mu[0]) );
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(mcMu[i][2]-mcMu[i][0])/mcMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2Mu[1]+=v*v; }   
    fprintf( outTxt, "%+8.2f", -sqrt(sumw2Mu[1]) );
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(sigMu[i][1]-sigMu[i][0])/sigMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2sigMu[0]+=v*v; }   
    fprintf( outTxt, "%+8.2f", sqrt(sumw2sigMu[0]) );
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(sigMu[i][2]-sigMu[i][0])/sigMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2sigMu[1]+=v*v; }   
    fprintf( outTxt, "%+8.2f", -sqrt(sumw2sigMu[1]) );
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(bkgMu[i][1]-bkgMu[i][0])/bkgMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2bkgMu[0]+=v*v; }   
    fprintf( outTxt, "%+8.2f", sqrt(sumw2bkgMu[0]) );
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(bkgMu[i][2]-bkgMu[i][0])/bkgMu[i][0]*100; fprintf( outTxt, "%+8.2f", v ); sumw2bkgMu[1]+=v*v; }   
    fprintf( outTxt, "%+8.2f", -sqrt(sumw2bkgMu[1]) );
    fprintf( outTxt, "\n\n");
    fclose( outTxt );
}
