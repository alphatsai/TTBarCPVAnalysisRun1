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
#include "checkGenACP.C"
using namespace std;
void test2()
{
    string input="../results/20Jan_LepJet_GenACPCheck/TTJets_SemiLeptMGDecays.root";
    string output="../results/20Jan_LepJet_GenACPCheck/";
    TFile* f = new TFile(input.c_str());
    checkGenACP(f, "O2", 1, output);
    checkGenACP(f, "O3", 1, output);
    checkGenACP(f, "O4", 1, output);
    checkGenACP(f, "O7", 1, output);
    checkGenACP(f, "O2", 0, output);
    checkGenACP(f, "O3", 0, output);
    checkGenACP(f, "O4", 0, output);
    checkGenACP(f, "O7", 0, output);

}
