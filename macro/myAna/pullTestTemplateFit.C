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
#include "templateLFit.C"
using namespace std;

const int nCh=2;
int el=0;
int mu=1;
int nExp=1000;
string output="../results/15Dec_LepJet_MCDATA/FitResualts2";
string tmpPath="../results/15Dec_LepJet_MCDATA/TemplateSyst_EvtChi2_Top_Leptonic_Mbl.root";
string name="EvtChi2_Top_Leptonic_Mbl";
//string tmpPath="../results/15Dec_LepJet_MCDATA/TemplateSyst_EvtChi2_Top_Hadronic_Mass.root";
//string name="EvtChi2_Top_Hadronic_Mass";
//string tmpPath="../results/15Dec_LepJet_MCDATA/TemplateSyst_EvtChi2_Ht.root";
//string name="EvtChi2_Ht";

void pullTestTemplateFit()
{
    // EvtChi2_Top_Leptonic_Mbl
    float nSig[nCh]={30406.,36560.};
    float nBkg[nCh]={2063.,1829.};
    float eSig[nCh]={265.7,274.8};
    float eBkg[nCh]={205.3,201.4};
    string xTitle[nCh]={"M(ej) in e-ch", "M(#mub) in #mu-ch"};
    string yTitle="Events";

    //// EvtChi2_Top_Hadronic_Mass 
    //float nSig[nCh]={29743.,36398.};
    //float nBkg[nCh]={2727.,1992.};
    //float eSig[nCh]={433.0,394.1};
    //float eBkg[nCh]={400.2,346.5};
    //string xTitle[nCh]={"M_{top}(bjj) in e-ch", "M_{top}(bjj) in #mu-ch"};
    //string yTitle="Events";

    //// EvtChi2_Ht 
    //float nSig[nCh]={30891.,37905.6};
    //float nBkg[nCh]={1578.4,423.};
    //float eSig[nCh]={428.7,779.6};
    //float eBkg[nCh]={391.4,733.0};
    //string xTitle[nCh]={"H_{T} in e-ch", "H_{T} in #mu-ch"};
    //string yTitle="Events";

    TFile* f = new TFile(tmpPath.c_str());
    for( int ch=0; ch<nCh; ch++ )
    {
        pullTest(  f, name, ch+1, nExp, nSig[ch], nBkg[ch], eSig[ch], eBkg[ch], output);
        fitter_LH( f, name, ch+1, output, xTitle[ch], yTitle );
    }
}
