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
#include "sumTemplateInfo.C"
using namespace std;
void test()
{
    //const int nElSyst=6;
    //const int nMuSyst=7;
    string input="/Users/Alpha/myAna/TTBarCPV/TTBarCPVAnalysisRun1/macro/results/11Jan_LepJet_MCDATA/FitResults2/TemplateSyst_EvtChi2_Top_Leptonic_Mbl.root";
    //string systNameEl[nElSyst]={"Stat", "PU", "JER", "BTagSF", "TopPT", "elID"};
    //string systNameMu[nMuSyst]={"Stat", "PU", "JER", "BTagSF", "TopPT", "muID", "muISO"};
    TFile* f = new TFile(input.c_str());
    //drawFittedStack( f, "EvtChi2_Top_Leptonic_Mbl", systNameEl, nElSyst, 1 );
    sumTemplateInfo(f, "EvtChi2_Top_Leptonic_Mbl", "../results/11Jan_LepJet_MCDATA/FitResults2/Syst/", "M(e+j) [GeV]", "M(#mu+j) [GeV]", "Events", 1);

}
