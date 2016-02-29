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
    //string input="../results/29Jan_LepJet_GenACP/ACP0/TTJets_SemiLeptMGDecays.root";
    //string output="../results/29Jan_LepJet_GenACP/ACP0";
    //TFile* f = new TFile(input.c_str());
    //checkGenACP(f, "O2", 1, output);
    //checkGenACP(f, "O3", 1, output);
    //checkGenACP(f, "O4", 1, output);
    //checkGenACP(f, "O7", 1, output);
    //checkGenACP(f, "O2", 0, output);
    //checkGenACP(f, "O3", 0, output);
    //checkGenACP(f, "O4", 0, output);
    //checkGenACP(f, "O7", 0, output);

    //string dirName="../results/29Jan_LepJet_GenACP";
    //string dirACP[11]={ "ACP-15", "ACP-10", "ACP-7.5", "ACP-5", "ACP-2.5", 
    //                    "ACP0",
    //                    "ACP+2.5", "ACP+5", "ACP+7.5", "ACP+10", "ACP+15"};
    // // Electron
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 1, 0.52 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 1, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 1, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 1, 0.72 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 0, 0.45 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 0, 0.26 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 0, 0.24 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 0, 0.70 );
    //
    // // Muon 
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 1, 0.52 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 1, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 1, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 1, 0.72 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 0, 0.45 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 0, 0.26 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 0, 0.24 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 0, 0.70 );

    //// Combined
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 1, 0.52 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 1, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 1, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 1, 0.72 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 0, 0.45 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 0, 0.26 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 0, 0.24 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 0, 0.70 );
    
    //for( int i=0; i<11; i++ )
    //{
    //    string input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
    //    string output = dirName+"/"+dirACP[i];
    //    TFile* f = new TFile(input.c_str());
    //    checkGenACP(f, "O2", 1, output);
    //    checkGenACP(f, "O3", 1, output);
    //    checkGenACP(f, "O4", 1, output);
    //    checkGenACP(f, "O7", 1, output);
    //    checkGenACP(f, "O2", 0, output);
    //    checkGenACP(f, "O3", 0, output);
    //    checkGenACP(f, "O4", 0, output);
    //    checkGenACP(f, "O7", 0, output);
    //    f->Close();
    //}

    string dirACP[3]={  "ACP-5", 
                        "ACP0",
                        "ACP+5"};
    string dirName="../results/25Feb_LepJet_GenACP_parton";
    for( int i=0; i<11; i++ )
    {
        string input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
        string output = dirName+"/"+dirACP[i];
        //string input  = dirName+"/TTJets_SemiLeptMGDecays.root";
        //string output = dirName;
        TFile* f = new TFile(input.c_str());
        checkGenACP(f, "O2", 1, output);
        checkGenACP(f, "O3", 1, output);
        checkGenACP(f, "O4", 1, output);
        checkGenACP(f, "O7", 1, output);
        checkGenACP(f, "O2", 0, output);
        checkGenACP(f, "O3", 0, output);
        checkGenACP(f, "O4", 0, output);
        checkGenACP(f, "O7", 0, output);
        checkParACP(f, "O2", 1, output);
        checkParACP(f, "O3", 1, output);
        checkParACP(f, "O4", 1, output);
        checkParACP(f, "O7", 1, output);
        checkParACP(f, "O2", 0, output);
        checkParACP(f, "O3", 0, output);
        checkParACP(f, "O4", 0, output);
        checkParACP(f, "O7", 0, output);
        f->Close();
    }
}
