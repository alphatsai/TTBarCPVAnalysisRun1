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

#include "mkPlotForPredictedNonZeroACPScan.C"

using namespace std;

bool printTitle=false;
const int NOBS=4; 
const int NCH=3;

std::string fin="../results/30Aug_LepJet_AddO3O4/Final_PredictionTree.root";
std::string outpath="../results/30Aug_LepJet_AddO3O4";

void toyResults()
{
    TFile* f = new TFile(fin.c_str()); 
    for( int iobs=0; iobs<NOBS; iobs++ ){
        for( int ich=0; ich<NCH; ich++)
        {
            mkPlotForPrediction(f, outpath, "MC",      iobs, ich, printTitle);
            mkPlotForPrediction(f, outpath, "SudoExp", iobs, ich, printTitle);
            mkPlotForPrediction(f, outpath, "SudoSig", iobs, ich, printTitle);
            mkPull(f, (outpath+"/pull").c_str(), "SudoExp", iobs, ich);
            mkPull(f, (outpath+"/pull").c_str(), "SudoSig", iobs, ich);
        }
    }


}
