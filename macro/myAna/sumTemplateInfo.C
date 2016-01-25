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
void sumTemplateInfo( TFile* f, string name, string output=".", string xElTitle="", string xMuTitle="", string yTitle="Events", bool logy=true )
{
    //const int nElSyst=9, nMuSyst=10;
    //string systNameEl[nElSyst]={"Stat", "TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID"};
    //string systNameMu[nMuSyst]={"Stat", "TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "muID", "muISO"};
    //const int nElSyst=8, nMuSyst=9;
    //string systNameEl[nElSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID"};
    //string systNameMu[nMuSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "muID", "muISO"};
    const int nElSyst=10, nMuSyst=10, nSyst=10;
    string systNameEl[nElSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};
    string systNameMu[nMuSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};
    string systName[nSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};

    float *mcEl[nElSyst], *sigEl[nElSyst], *bkgEl[nElSyst];
    float *mcMu[nMuSyst], *sigMu[nMuSyst], *bkgMu[nMuSyst];
    float sumw2El[2]={0,0}, sumw2sigEl[2]={0,0}, sumw2bkgEl[2]={0,0};
    float sumw2Mu[2]={0,0}, sumw2sigMu[2]={0,0}, sumw2bkgMu[2]={0,0};
  
    FILE* outTxt;
    outTxt = fopen((output+"/SystUncertainties_MC.txt").c_str(),"w");
    fprintf( outTxt, "Electron channel\n" );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<nElSyst; i++ )
    {
        mcEl[i]  = drawSyst( f, "MC_El",    systNameEl[i], output, xElTitle, yTitle, logy);
        sigEl[i] = drawSyst( f, "SigMC_El", systNameEl[i], output, xElTitle, yTitle, logy);
        bkgEl[i] = drawSyst( f, "BkgMC_El", systNameEl[i], output, xElTitle, yTitle, logy);
        fprintf( outTxt, "%9s", systNameEl[i].c_str() );
    }
    fprintf( outTxt, "%9s", "Syst" );
    fprintf( outTxt, "\nMC:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(mcEl[i][1]-mcEl[i][0])/mcEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2El[0]+=v*v;} }  
    fprintf( outTxt, "%+9.2f", sqrt(sumw2El[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(mcEl[i][2]-mcEl[i][0])/mcEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2El[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2El[1]) );
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(sigEl[i][1]-sigEl[i][0])/sigEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sigEl[0]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2sigEl[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(sigEl[i][2]-sigEl[i][0])/sigEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sigEl[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2sigEl[1]) );
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(bkgEl[i][1]-bkgEl[i][0])/bkgEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkgEl[0]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2bkgEl[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nElSyst; i++ ){ float v=(bkgEl[i][2]-bkgEl[i][0])/bkgEl[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkgEl[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2bkgEl[1]) );
    fprintf( outTxt, "\n\n");
    //fprintf( outTxt, "\nMC:");
    //fprintf( outTxt, "\n%9s", "+1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (mcEl[i][1]-mcEl[i][0])/mcEl[i][0]*100 ); }   
    //fprintf( outTxt, "\n%9s", "-1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (mcEl[i][2]-mcEl[i][0])/mcEl[i][0]*100 ); }   
    //fprintf( outTxt, "\nSig:");
    //fprintf( outTxt, "\n%9s", "+1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (sigEl[i][1]-sigEl[i][0])/sigEl[i][0]*100 ); }   
    //fprintf( outTxt, "\n%9s", "-1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (sigEl[i][2]-sigEl[i][0])/sigEl[i][0]*100 ); }   
    //fprintf( outTxt, "\nBkg:");
    //fprintf( outTxt, "\n%9s", "+1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (bkgEl[i][1]-bkgEl[i][0])/bkgEl[i][0]*100 ); }   
    //fprintf( outTxt, "\n%9s", "-1sigma" );
    //for( int i=0; i<nElSyst; i++ ){ fprintf( outTxt, "%+9.2f", (bkgEl[i][2]-bkgEl[i][0])/bkgEl[i][0]*100 ); }   
    //fprintf( outTxt, "\n\n");

    fprintf( outTxt, "Muon channel\n" );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<nMuSyst; i++ )
    {
        mcMu[i]  = drawSyst( f, "MC_Mu",    systNameMu[i], output, xMuTitle, yTitle, logy);
        sigMu[i] = drawSyst( f, "SigMC_Mu", systNameMu[i], output, xMuTitle, yTitle, logy);
        bkgMu[i] = drawSyst( f, "BkgMC_Mu", systNameMu[i], output, xMuTitle, yTitle, logy);
        fprintf( outTxt, "%9s", systNameMu[i].c_str() );
    }
    fprintf( outTxt, "%9s", "Syst" );
    fprintf( outTxt, "\nMC:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(mcMu[i][1]-mcMu[i][0])/mcMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2Mu[0]+=v*v;} }  
    fprintf( outTxt, "%+9.2f", sqrt(sumw2Mu[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(mcMu[i][2]-mcMu[i][0])/mcMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2Mu[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2Mu[1]) );
    fprintf( outTxt, "\nSig:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(sigMu[i][1]-sigMu[i][0])/sigMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sigMu[0]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2sigMu[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(sigMu[i][2]-sigMu[i][0])/sigMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2sigMu[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2sigMu[1]) );
    fprintf( outTxt, "\nBkg:");
    fprintf( outTxt, "\n%9s", "+1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(bkgMu[i][1]-bkgMu[i][0])/bkgMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkgMu[0]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", sqrt(sumw2bkgMu[0]) );
    fprintf( outTxt, "\n%9s", "-1sigma" );
    for( int i=0; i<nMuSyst; i++ ){ float v=(bkgMu[i][2]-bkgMu[i][0])/bkgMu[i][0]*100; fprintf( outTxt, "%+9.2f", v ); if( i!=0 ){ sumw2bkgMu[1]+=v*v;} }   
    fprintf( outTxt, "%+9.2f", -sqrt(sumw2bkgMu[1]) );
    fprintf( outTxt, "\n\n");
    fclose( outTxt );
    
    drawFittedStack( f, name, systNameEl, nElSyst, 1, output, xElTitle, yTitle );
    drawFittedStack( f, name, systNameMu, nMuSyst, 2, output, xMuTitle, yTitle );

    bool unBlind=false;

    FILE* outTxt1;
    outTxt1 = fopen((output+"/FinalAcpResults.txt").c_str(),"w");
    getSubtractBkgResults( f, outTxt1, "O2", nElSyst, systNameEl, 0 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O3", nElSyst, systNameEl, 0 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O4", nElSyst, systNameEl, 0 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O7", nElSyst, systNameEl, 0 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O2", nMuSyst, systNameMu, 1 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O3", nMuSyst, systNameMu, 1 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O4", nMuSyst, systNameMu, 1 , unBlind );
    getSubtractBkgResults( f, outTxt1, "O7", nMuSyst, systNameMu, 1 , unBlind );
    fclose( outTxt1 );

    FILE* outTxt2;
    outTxt2 = fopen((output+"/FinalAcpResultsCombined.txt").c_str(),"w");
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, output, "O2", "O_{2}", unBlind );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, output, "O3", "O_{3}", unBlind );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, output, "O4", "O_{4}", unBlind );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, output, "O7", "O_{7}", unBlind );
    fclose( outTxt2 );
}
