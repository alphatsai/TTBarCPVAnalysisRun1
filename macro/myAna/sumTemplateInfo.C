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
    const int nElSyst=5, nMuSyst=6;
    float *mcEl[nElSyst], *sigEl[nElSyst], *bkgEl[nElSyst];
    float *mcMu[nMuSyst], *sigMu[nMuSyst], *bkgMu[nMuSyst];
    string systNameEl[nElSyst]={"PU", "JER", "BTagSF", "TopPT", "elID"};
    string systNameMu[nMuSyst]={"PU", "JER", "BTagSF", "TopPT", "muID", "muISO"};
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
    fprintf( outTxt, "\nMC:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (mcMu[i][1]-mcMu[i][0])/mcMu[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (mcMu[i][2]-mcMu[i][0])/mcMu[i][0]*100 ); }   
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (sigMu[i][1]-sigMu[i][0])/sigMu[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (sigMu[i][2]-sigMu[i][0])/sigMu[i][0]*100 ); }   
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%8s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (bkgMu[i][1]-bkgMu[i][0])/bkgMu[i][0]*100 ); }   
    fprintf( outTxt, "\n%8s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ fprintf( outTxt, "%+8.2f", (bkgMu[i][2]-bkgMu[i][0])/bkgMu[i][0]*100 ); }   
    fprintf( outTxt, "\n\n");
    fclose( outTxt );
}
