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
    // Check syst unc.
    string dirName="../results/07Apr_GenAcp/ACP0";
    const int size=12;
    string dirSyst[size]={ "BTagSF", "JES", "JER", "elID", "muID", "muISO", "PU", "TopPT", "topMass", "topQ2", "topMEPS", "topPOWHEG"};

    FILE* outTxt;
    string hCat="EvtChi2";
    outTxt = fopen((dirName+"/DiluteSystUncertainties"+hCat+".txt").c_str(),"w");
    getChagnesSyst( outTxt, dirName, dirSyst, size, "O2", hCat );
    getChagnesSyst( outTxt, dirName, dirSyst, size, "O3", hCat );
    getChagnesSyst( outTxt, dirName, dirSyst, size, "O4", hCat );
    getChagnesSyst( outTxt, dirName, dirSyst, size, "O7", hCat );
    getChagnesSyst( outTxt, dirName, dirSyst, size, "Oa", hCat );
    getChagnesSyst( outTxt, dirName, dirSyst, size, "Ob", hCat );
    fclose( outTxt );

    const int size2 = size*2;
    string dirACP[size2]={ "BTagSFdown", "elIDdown", "JESdown", "JERdown", "muIDdown", "muISOdown", "PUdown", "TopPTdown", "topMassdown", "topQ2down", "topMEPSdown", "topPOWHEGdown",
                           "BTagSFup",   "elIDup",   "JESup",   "JERup",   "muIDup",   "muISOup",   "PUup",   "TopPTup",   "topMassup",   "topQ2up",   "topMEPSup",   "topPOWHEG" };
    for( int i=0; i<size2+1; i++ )
    {
        string input, output;
        if( i < size )
        {
            input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
            output = dirName+"/"+dirACP[i];
        }else{
            input  = dirName+"/nominal/TTJets_SemiLeptMGDecays.root";
            output = dirName+"/nominal";
        }
        TFile* f = new TFile(input.c_str());
        checkGenACP(f, "Oa", 1, output);
        checkGenACP(f, "Ob", 1, output);
        checkGenACP(f, "O2", 1, output);
        checkGenACP(f, "O3", 1, output);
        checkGenACP(f, "O4", 1, output);
        checkGenACP(f, "O7", 1, output);
        checkGenACP(f, "Oa", 0, output);
        checkGenACP(f, "Ob", 0, output);
        checkGenACP(f, "O2", 0, output);
        checkGenACP(f, "O3", 0, output);
        checkGenACP(f, "O4", 0, output);
        checkGenACP(f, "O7", 0, output);
        f->Close();
    }

    //// Make linearity test
    ////string input="../results/09Mar_LepJet_GenACPCheck/ACP0/TTJets_SemiLeptMGDecays.root";
    ////string output="../results/09Mar_LepJet_GenACPCheck/ACP0";
    ////TFile* f = new TFile(input.c_str());
    ////checkGenACP(f, "O2", 1, output);
    ////checkGenACP(f, "O3", 1, output);
    ////checkGenACP(f, "O4", 1, output);
    ////checkGenACP(f, "O7", 1, output);
    ////checkGenACP(f, "O2", 0, output);
    ////checkGenACP(f, "O3", 0, output);
    ////checkGenACP(f, "O4", 0, output);
    ////checkGenACP(f, "O7", 0, output);

    ////string dirName="../results/09Mar_LepJet_GenACPCheck";
    //string dirName="../results/16Mar_LepJet_ObOa/GenACPCheck";
    //string dirACP[11]={ "ACP-15", "ACP-10", "ACP-7.5", "ACP-5", "ACP-2.5", 
    //                    "ACP0",
    //                    "ACP+2.5", "ACP+5", "ACP+7.5", "ACP+10", "ACP+15"};
    // // Electron
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 1, 0.54 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 1, 0.36 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 1, 0.34 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 1, 0.73 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 0, 0.47 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 0, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 0, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 0, 0.72 );
    //
    // // Muon 
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 1, 0.54 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 1, 0.35 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 1, 0.34 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 1, 0.73 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 0, 0.47 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 0, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 0, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 0, 0.72 );

    //// Combined
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 1, 0.54 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 1, 0.36 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 1, 0.34 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 1, 0.73 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 0, 0.47 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 0, 0.31 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 0, 0.29 );
    // mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 0, 0.72 );

    //// Without Mlb cut
    //// // Electron
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 1, 0.52 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 1, 0.31 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 1, 0.29 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 1, 0.72 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 1, 0, 0.45 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 1, 0, 0.26 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 1, 0, 0.24 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 1, 0, 0.70 );
    ////
    //// // Muon 
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 1, 0.52 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 1, 0.31 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 1, 0.29 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 1, 0.72 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 2, 0, 0.45 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 2, 0, 0.26 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 2, 0, 0.24 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 2, 0, 0.70 );

    ////// Combined
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 1, 0.52 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 1, 0.31 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 1, 0.29 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 1, 0.72 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O2", "o_{2}", 0, 0, 0.45 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O3", "o_{3}", 0, 0, 0.26 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O4", "o_{4}", 0, 0, 0.24 );
    //// mkPlotSlopeACP( dirName, dirACP, 11, dirName, "O7", "o_{7}", 0, 0, 0.70 );
    //
    //for( int i=0; i<11; i++ )
    //{
    //    string input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
    //    string output = dirName+"/"+dirACP[i];
    //    TFile* f = new TFile(input.c_str());
    //    checkGenACP(f, "Oa", 1, output);
    //    checkGenACP(f, "Ob", 1, output);
    //    checkGenACP(f, "O2", 1, output);
    //    checkGenACP(f, "O3", 1, output);
    //    checkGenACP(f, "O4", 1, output);
    //    checkGenACP(f, "O7", 1, output);
    //    checkGenACP(f, "Oa", 0, output);
    //    checkGenACP(f, "Ob", 0, output);
    //    checkGenACP(f, "O2", 0, output);
    //    checkGenACP(f, "O3", 0, output);
    //    checkGenACP(f, "O4", 0, output);
    //    checkGenACP(f, "O7", 0, output);
    //    f->Close();
    //}
    //TFile* f = new TFile((dirName+"/ACP0/TTJets_SemiLeptMGDecays.root").c_str());
    //checkGenACP(f, "O2", 1, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O3", 1, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O4", 1, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O7", 1, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O2", 0, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O3", 0, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O4", 0, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O7", 0, dirName+"/BB",    "EvtBB", "GenBB");
    //checkGenACP(f, "O2", 1, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O3", 1, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O4", 1, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O7", 1, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O2", 0, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O3", 0, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O4", 0, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //checkGenACP(f, "O7", 0, dirName+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");

    //// For checking ttbar responce 
    //string dirACP[3]={  "ACP-5", 
    //                    "ACP0",
    //                    "ACP+5"};
    //string dirName="../results/10Mar_LepJet_GenACPCheck";
    //for( int i=0; i<11; i++ )
    //{
    //    string input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
    //    string output = dirName+"/"+dirACP[i];
    //    //string input  = dirName+"/TTJets_SemiLeptMGDecays.root";
    //    //string output = dirName;
    //    TFile* f = new TFile(input.c_str());
    //    checkGenACP(f, "O2", 1, output);
    //    checkGenACP(f, "O3", 1, output);
    //    checkGenACP(f, "O4", 1, output);
    //    checkGenACP(f, "O7", 1, output);
    //    checkGenACP(f, "O2", 0, output);
    //    checkGenACP(f, "O3", 0, output);
    //    checkGenACP(f, "O4", 0, output);
    //    checkGenACP(f, "O7", 0, output);
    //    checkParACP(f, "O2", 1, output);
    //    checkParACP(f, "O3", 1, output);
    //    checkParACP(f, "O4", 1, output);
    //    checkParACP(f, "O7", 1, output);
    //    checkParACP(f, "O2", 0, output);
    //    checkParACP(f, "O3", 0, output);
    //    checkParACP(f, "O4", 0, output);
    //    checkParACP(f, "O7", 0, output);
    //    f->Close();
    //}

    //// For checking W jets 
    ////string dirACP[1]={"ACP0"};
    //string dirName="../results/10Mar_LepJet_GenACPCheck";
    ////string dirName="../results/16Mar_LepJet_ObOa";
    ////for( int i=0; i<11; i++ )
    ////{
    //    //string input  = dirName+"/"+dirACP[i]+"/TTJets_SemiLeptMGDecays.root";
    //    //string output = dirName+"/"+dirACP[i];
    //    string input  = dirName+"/TTJets_SemiLeptMGDecays.root";
    //    string output = dirName;
    //    TFile* f = new TFile(input.c_str());
    //    checkGenACP(f, "O2", 1, output);
    //    checkGenACP(f, "O3", 1, output);
    //    checkGenACP(f, "O4", 1, output);
    //    checkGenACP(f, "O7", 1, output);
    //    checkGenACP(f, "O2", 0, output);
    //    checkGenACP(f, "O3", 0, output);
    //    checkGenACP(f, "O4", 0, output);
    //    checkGenACP(f, "O7", 0, output);
    //    checkGenACP(f, "O2", 1, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O3", 1, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O4", 1, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O7", 1, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O2", 0, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O3", 0, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O4", 0, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O7", 0, output+"/BB",    "EvtBB", "GenBB");
    //    checkGenACP(f, "O2", 1, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O3", 1, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O4", 1, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O7", 1, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O2", 0, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O3", 0, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O4", 0, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    checkGenACP(f, "O7", 0, output+"/J1Q1", "EvtJ1Q1", "GenJ1Q1");
    //    f->Close();
    ////}
}
