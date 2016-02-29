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
    bool unblind=true;
    string input="/Users/Alpha/myAna/TTBarCPV/TTBarCPVAnalysisRun1/macro/results/11Jan_LepJet_MCDATA/FitResults4/TemplateSyst_EvtChi2_Top_Leptonic_Mbl.root";
    string inputCR="/Users/Alpha/myAna/TTBarCPV/TTBarCPVAnalysisRun1/macro/results/11Jan_LepJet_MCDATA/0bCR/Final_histograms_SemiLeptanic.root";
    TFile* f   = new TFile(input.c_str());
    TFile* fCR = new TFile(inputCR.c_str());
    sumTemplateInfo(f, "EvtChi2_Top_Leptonic_Mbl", "../results/11Jan_LepJet_MCDATA/FitResults4/Syst/Unblind_CR0b", "M_{eb} [GeV]", "M_{#mub} [GeV]", "Events", fCR, unblind );
    //sumTemplateInfo(f, "EvtChi2_Top_Leptonic_Mbl", "../results/11Jan_LepJet_MCDATA/FitResults4/Syst/CR0b_topMassReScale", "M_{eb} [GeV]", "M_{#mub} [GeV]", "Events", fCR );
    //sumTemplateInfo(f, "EvtChi2_Top_Leptonic_Mbl", "../results/11Jan_LepJet_MCDATA/FitResults4/Syst/CR0b", "M_{eb} [GeV]", "M_{#mub} [GeV]", "Events", fCR );
    //sumTemplateInfo(f, "EvtChi2_Top_Leptonic_Mbl", "../results/11Jan_LepJet_MCDATA/FitResults4/Syst/topMassReScale", "M_{eb} [GeV]", "M_{#mub} [GeV]", "Events");
}
