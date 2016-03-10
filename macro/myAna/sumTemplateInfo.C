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
//void sumTemplateInfo( TFile* f, string name, string output=".", string xElTitle="", string xMuTitle="", string yTitle="Events" )
void sumTemplateInfo( TFile* f, string name, string output=".", string xElTitle="", string xMuTitle="", string yTitle="Events", TFile* fCR=NULL, bool unBlind=false )
{
    bool logySyst=true;
    bool doTopMassUncRescale=true; 
    int nPDFup=34;
    int nPDFdown=35; 
    //const int nElSyst=9, nMuSyst=10;
    //string systNameEl[nElSyst]={"Stat", "TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID"};
    //string systNameMu[nMuSyst]={"Stat", "TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "muID", "muISO"};
    //const int nElSyst=8, nMuSyst=9;
    //string systNameEl[nElSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID"};
    //string systNameMu[nMuSyst]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "muID", "muISO"};
    const int nElSyst=12, nMuSyst=12, nSyst=12, nPDF=50;
    string systNameEl[nElSyst]={"TopMatch", "TopScale", "TopPT", "TopGen", "TopMass", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};
    string systNameMu[nMuSyst]={"TopMatch", "TopScale", "TopPT", "TopGen", "TopMass", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};
    string systName[nSyst]={"TopMatch", "TopScale", "TopPT", "TopGen", "TopMass", "PU", "JER", "JES", "BTagSF", "elID", "muID", "muISO"};

    float *mcEl[nElSyst], *sigEl[nElSyst], *bkgEl[nElSyst];
    float *mcMu[nMuSyst], *sigMu[nMuSyst], *bkgMu[nMuSyst];
    float *mcElPDF[nPDF], *sigElPDF[nPDF], *bkgElPDF[nPDF];
    float *mcMuPDF[nPDF], *sigMuPDF[nPDF], *bkgMuPDF[nPDF];
    float sumw2El[2]={0,0}, sumw2sigEl[2]={0,0}, sumw2bkgEl[2]={0,0};
    float sumw2Mu[2]={0,0}, sumw2sigMu[2]={0,0}, sumw2bkgMu[2]={0,0};
  
    FILE* outTxt;
    outTxt = fopen((output+"/SystUncertainties_MC.txt").c_str(),"w");
    fprintf( outTxt, "!!!! Without TopMass variation re-scale !!!!\n" );
    fprintf( outTxt, "Electron channel\n" );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<nElSyst; i++ )
    {
        mcEl[i]  = drawSyst( f, "MC_El",    systNameEl[i], output, xElTitle, yTitle, logySyst);
        sigEl[i] = drawSyst( f, "SigMC_El", systNameEl[i], output, xElTitle, yTitle, logySyst);
        bkgEl[i] = drawSyst( f, "BkgMC_El", systNameEl[i], output, xElTitle, yTitle, logySyst);
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

    // Print PDF
    int iminElmc(0), iminElsig(0), iminElbkg(0);
    int imaxElmc(0), imaxElsig(0), imaxElbkg(0);
    float maxElmc(0), maxElsig(0), maxElbkg(0); 
    float minElmc(0), minElsig(0), minElbkg(0);
    float max(0), min(0);
    fprintf( outTxt, "%6s", "PDF");
    for( int i=0; i<nPDF; i++ )
    {
        mcElPDF[i]  = getEvtPDF( f, "MC_El",    i+1 );
        sigElPDF[i] = getEvtPDF( f, "SigMC_El", i+1 );
        bkgElPDF[i] = getEvtPDF( f, "BkgMC_El", i+1 );
        fprintf( outTxt, "%6s", ("P_"+int2str(i+1)).c_str() );
    }
    fprintf( outTxt, "\n%6s", "MC:");
    for( int i=0; i<nPDF; i++ ){ float v=(mcElPDF[i][1]-mcElPDF[i][0])/mcElPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminElmc=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxElmc=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxElmc, max,  sqrt(sumw2El[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminElmc, min, -1*sqrt(sumw2El[1]+min*min));
    maxElmc=max; minElmc=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n%6s", "Sig:");
    for( int i=0; i<nPDF; i++ ){ float v=(sigElPDF[i][1]-sigElPDF[i][0])/sigElPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminElsig=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxElsig=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxElsig, max,  sqrt(sumw2sigEl[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminElsig, min, -1*sqrt(sumw2sigEl[1]+min*min));
    maxElsig=max; minElsig=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n%6s", "Bkg:");
    for( int i=0; i<nPDF; i++ ){ float v=(bkgElPDF[i][1]-bkgElPDF[i][0])/bkgElPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminElbkg=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxElbkg=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxElbkg, max,  sqrt(sumw2bkgEl[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminElbkg, min, -1*sqrt(sumw2bkgEl[1]+min*min));
    float bkgElMax=(bkgElPDF[imaxElsig][1]-bkgElPDF[imaxElsig][0])/bkgElPDF[imaxElsig][0]*100;
    float bkgElMin=(bkgElPDF[iminElsig][1]-bkgElPDF[iminElsig][0])/bkgElPDF[iminElsig][0]*100;
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxElsig, bkgElMax,    sqrt(sumw2bkgEl[0]+bkgElMax*bkgElMax));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminElsig, bkgElMin, -1*sqrt(sumw2bkgEl[1]+bkgElMin*bkgElMin));
    maxElbkg=max; minElbkg=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n\n\n");

    fprintf( outTxt, "Muon channel\n" );
    fprintf( outTxt, "%9s", " ");
    for( int i=0; i<nMuSyst; i++ )
    {
        mcMu[i]  = drawSyst( f, "MC_Mu",    systNameMu[i], output, xMuTitle, yTitle, logySyst);
        sigMu[i] = drawSyst( f, "SigMC_Mu", systNameMu[i], output, xMuTitle, yTitle, logySyst);
        bkgMu[i] = drawSyst( f, "BkgMC_Mu", systNameMu[i], output, xMuTitle, yTitle, logySyst);
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

    // Print PDF
    int iminMumc(0), iminMusig(0), iminMubkg(0);
    int imaxMumc(0), imaxMusig(0), imaxMubkg(0);
    float maxMumc(0), maxMusig(0), maxMubkg(0); 
    float minMumc(0), minMusig(0), minMubkg(0);
    max=0; min=0; 
    fprintf( outTxt, "%6s", "PDF");
    for( int i=0; i<nPDF; i++ )
    {
        mcMuPDF[i]  = getEvtPDF( f, "MC_Mu",    i+1 );
        sigMuPDF[i] = getEvtPDF( f, "SigMC_Mu", i+1 );
        bkgMuPDF[i] = getEvtPDF( f, "BkgMC_Mu", i+1 );
        fprintf( outTxt, "%6s", ("P_"+int2str(i+1)).c_str() );
    }
    fprintf( outTxt, "\n%6s", "MC:");
    for( int i=0; i<nPDF; i++ ){ float v=(mcMuPDF[i][1]-mcMuPDF[i][0])/mcMuPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminMumc=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxMumc=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxMumc, max,  sqrt(sumw2Mu[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminMumc, min, -1*sqrt(sumw2Mu[1]+min*min));
    maxMumc=max; minMumc=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n%6s", "Sig:");
    for( int i=0; i<nPDF; i++ ){ float v=(sigMuPDF[i][1]-sigMuPDF[i][0])/sigMuPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminMusig=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxMusig=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxMusig, max,  sqrt(sumw2sigMu[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminMusig, min, -1*sqrt(sumw2sigMu[1]+min*min));
    maxMusig=max; minMusig=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n%6s", "Bkg:");
    for( int i=0; i<nPDF; i++ ){ float v=(bkgMuPDF[i][1]-bkgMuPDF[i][0])/bkgMuPDF[i][0]*100; fprintf( outTxt, "%+6.2f", v ); if( v < min && fabs(v) < 10 ){ min=v; iminMubkg=i; } if( v > max && fabs(v) < 10 ){ max=v; imaxMubkg=i; } }  
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxMubkg, max,    sqrt(sumw2bkgMu[0]+max*max));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminMubkg, min, -1*sqrt(sumw2bkgMu[1]+min*min));
    float bkgMuMax=(bkgMuPDF[imaxMusig][1]-bkgMuPDF[imaxMusig][0])/bkgMuPDF[imaxMusig][0]*100;
    float bkgMuMin=(bkgMuPDF[iminMusig][1]-bkgMuPDF[iminMusig][0])/bkgMuPDF[iminMusig][0]*100;
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Up",   imaxMusig, bkgMuMax,    sqrt(sumw2bkgMu[0]+bkgMuMax*bkgMuMax));
    fprintf( outTxt, "\n%6s%6d%+6.2f, Total syst %+6.2f", "Down", iminMusig, bkgMuMin, -1*sqrt(sumw2bkgMu[1]+bkgMuMin*bkgMuMin));
    maxMubkg=max; minMubkg=min; 
    max=0; min=0; 
    fprintf( outTxt, "\n\n\n");
    fclose( outTxt );

    float minMlb=0, maxMlb=200;
    drawFittedStack( f, name, systNameEl, nElSyst, nPDFup, nPDFdown, 1, minMlb, maxMlb, output, xElTitle, yTitle, false, doTopMassUncRescale, false );
    drawFittedStack( f, name, systNameMu, nMuSyst, nPDFup, nPDFdown, 2, minMlb, maxMlb, output, xMuTitle, yTitle, false, doTopMassUncRescale, false );

    FILE* outTxt1;
    if( unBlind )
        outTxt1 = fopen((output+"/FinalAcpResults.txt").c_str(),"w");
    else
        outTxt1 = fopen((output+"/FinalBlindAcpResults.txt").c_str(),"w");
    getSubtractBkgResults( f, outTxt1, "O2", nElSyst, systNameEl, nPDFup, nPDFdown, 0, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O3", nElSyst, systNameEl, nPDFup, nPDFdown, 0, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O4", nElSyst, systNameEl, nPDFup, nPDFdown, 0, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O7", nElSyst, systNameEl, nPDFup, nPDFdown, 0, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O2", nMuSyst, systNameMu, nPDFup, nPDFdown, 1, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O3", nMuSyst, systNameMu, nPDFup, nPDFdown, 1, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O4", nMuSyst, systNameMu, nPDFup, nPDFdown, 1, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    getSubtractBkgResults( f, outTxt1, "O7", nMuSyst, systNameMu, nPDFup, nPDFdown, 1, minMlb, maxMlb, unBlind, doTopMassUncRescale, fCR );
    fclose( outTxt1 );

    FILE* outTxt2;
    if( unBlind )
        outTxt2 = fopen((output+"/FinalAcpResultsCombined.txt").c_str(),"w");
    else
        outTxt2 = fopen((output+"/FinalBlindAcpResultsCombined.txt").c_str(),"w");
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O2", "O_{2}", unBlind, 0, doTopMassUncRescale, fCR );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O3", "O_{3}", unBlind, 0, doTopMassUncRescale, fCR );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O4", "O_{4}", unBlind, 0, doTopMassUncRescale, fCR );
    getSubtractBkgResultsCombined( f, outTxt2, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O7", "O_{7}", unBlind, 0, doTopMassUncRescale, fCR );
    drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O2", "O_{2}", unBlind, 0, doTopMassUncRescale, fCR );
    drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O3", "O_{3}", unBlind, 0, doTopMassUncRescale, fCR );
    drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O4", "O_{4}", unBlind, 0, doTopMassUncRescale, fCR );
    drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O7", "O_{7}", unBlind, 0, doTopMassUncRescale, fCR );
    fclose( outTxt2 );

    if( !unBlind )
    {
        FILE* outTxt3;
        outTxt3 = fopen((output+"/FinalBlindAcpResultsCombined_Assume5.txt").c_str(),"w");
        getSubtractBkgResultsCombined( f, outTxt3, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O2", "O_{2}", false, 0.05, doTopMassUncRescale, fCR  );
        getSubtractBkgResultsCombined( f, outTxt3, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O3", "O_{3}", false, 0.05, doTopMassUncRescale, fCR  );
        getSubtractBkgResultsCombined( f, outTxt3, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O4", "O_{4}", false, 0.05, doTopMassUncRescale, fCR  );
        getSubtractBkgResultsCombined( f, outTxt3, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O7", "O_{7}", false, 0.05, doTopMassUncRescale, fCR  );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O2", "O_{2}", false, 0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O3", "O_{3}", false, 0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O4", "O_{4}", false, 0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O7", "O_{7}", false, 0.05, doTopMassUncRescale, fCR );
        fclose( outTxt3 );

        FILE* outTxt4;
        outTxt4 = fopen((output+"/FinalBlindAcpResultsCombined_Assume-5.txt").c_str(),"w");
        getSubtractBkgResultsCombined( f, outTxt4, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O2", "O_{2}", false, -0.05, doTopMassUncRescale, fCR );
        getSubtractBkgResultsCombined( f, outTxt4, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O3", "O_{3}", false, -0.05, doTopMassUncRescale, fCR );
        getSubtractBkgResultsCombined( f, outTxt4, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O4", "O_{4}", false, -0.05, doTopMassUncRescale, fCR );
        getSubtractBkgResultsCombined( f, outTxt4, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb,  output, "O7", "O_{7}", false, -0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O2", "O_{2}", false, -0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O3", "O_{3}", false, -0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O4", "O_{4}", false, -0.05, doTopMassUncRescale, fCR );
        drawSubtractBkgResultsCombined( f, systName, nSyst, nPDFup, nPDFdown, minMlb, maxMlb, output, "O7", "O_{7}", false, -0.05, doTopMassUncRescale, fCR );
        fclose( outTxt4 );
    }
}
