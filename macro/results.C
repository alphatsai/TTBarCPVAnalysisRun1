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
#include "plots.C"
#include "drawACP.C"
#include "stack.C"
using namespace std;
void results(TFile* f, string output="."){

    //drawObservable( f, "MC__Evt_O7Asym",      output, "O_{7}^{e+#mu}", 	"Combined", 0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym",      output, "O_{2}^{e+#mu}", 	"Combined", 0.05, 1);
    //drawObservable( f, "MC__Evt_O7Asym_Mu",   output, "O_{7}^{#mu}", 	"Muon", 	0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym_Mu",   output, "O_{2}^{#mu}", 	"Muon", 	0.05, 1);
    //drawObservable( f, "MC__Evt_O7Asym_El",   output, "O_{7}^{e}", 		"Electron", 0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym_El",   output, "O_{2}^{e}", 		"Electron", 0.05, 1);

    //drawObservableDist( f, output, "MC__Evt_O7", 	"O_{7}/M_{top}^{3}", "LepJets", 1);
    //drawObservableDist( f, output, "MC__Evt_O2", 	"O_{2}/M_{top}^{3}", "LepJets", 1);
    drawACP( f, 1, "", "BkgMC__Evt", "O2", output, "ACP", "O_{2}");
    drawACP( f, 1, "", "BkgMC__Evt", "O3", output, "ACP", "O_{3}");
    drawACP( f, 1, "", "BkgMC__Evt", "O4", output, "ACP", "O_{4}");
    drawACP( f, 1, "", "BkgMC__Evt", "O7", output, "ACP", "O_{7}");
    drawACP( f, 0, "", "TTJets_SemiLeptMGDecays__Evt", "O2", output, "ACP", "O_{2}");
    drawACP( f, 0, "", "TTJets_SemiLeptMGDecays__Evt", "O3", output, "ACP", "O_{3}");
    drawACP( f, 0, "", "TTJets_SemiLeptMGDecays__Evt", "O4", output, "ACP", "O_{4}");
    drawACP( f, 0, "", "TTJets_SemiLeptMGDecays__Evt", "O7", output, "ACP", "O_{7}");

    drawStack( f, "Evt_O2",    "O_{2}^{e+#mu}/M_{top}^{3}", "Events", output, 2, false );
    drawStack( f, "Evt_O2",    "O_{2}^{e+#mu}/M_{top}^{3}", "Events", output, 2, true );
    drawStack( f, "Evt_O3",    "O_{3}^{e+#mu}/M_{top}^{3}", "Events", output, 2, false );
    drawStack( f, "Evt_O3",    "O_{3}^{e+#mu}/M_{top}^{3}", "Events", output, 2, true );
    drawStack( f, "Evt_O4",    "O_{4}^{e+#mu}/M_{top}^{3}", "Events", output, 2, false );
    drawStack( f, "Evt_O4",    "O_{4}^{e+#mu}/M_{top}^{3}", "Events", output, 2, true );
    drawStack( f, "Evt_O7",    "O_{7}^{e+#mu}/M_{top}^{3}", "Events", output, 2, false );
    drawStack( f, "Evt_O7",    "O_{7}^{e+#mu}/M_{top}^{3}", "Events", output, 2, true );
    drawStack( f, "Evt_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}",   "Events", output, 2, false );
    drawStack( f, "Evt_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}",   "Events", output, 2, true );
    drawStack( f, "Evt_O3_Mu", "O_{3}^{#mu}/M_{top}^{3}",   "Events", output, 2, false );
    drawStack( f, "Evt_O3_Mu", "O_{3}^{#mu}/M_{top}^{3}",   "Events", output, 2, true );
    drawStack( f, "Evt_O4_Mu", "O_{4}^{#mu}/M_{top}^{3}",   "Events", output, 2, false );
    drawStack( f, "Evt_O4_Mu", "O_{4}^{#mu}/M_{top}^{3}",   "Events", output, 2, true );
    drawStack( f, "Evt_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}",   "Events", output, 2, false );
    drawStack( f, "Evt_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}",   "Events", output, 2, true );
    drawStack( f, "Evt_O2_El", "O_{2}^{e}/M_{top}^{3}",     "Events", output, 2, false );
    drawStack( f, "Evt_O2_El", "O_{2}^{e}/M_{top}^{3}",     "Events", output, 2, true );
    drawStack( f, "Evt_O3_El", "O_{3}^{e}/M_{top}^{3}",     "Events", output, 2, false );
    drawStack( f, "Evt_O3_El", "O_{3}^{e}/M_{top}^{3}",     "Events", output, 2, true );
    drawStack( f, "Evt_O4_El", "O_{4}^{e}/M_{top}^{3}",     "Events", output, 2, false );
    drawStack( f, "Evt_O4_El", "O_{4}^{e}/M_{top}^{3}",     "Events", output, 2, true );
    drawStack( f, "Evt_O7_El", "O_{7}^{e}/M_{top}^{3}",     "Events", output, 2, false );
    drawStack( f, "Evt_O7_El", "O_{7}^{e}/M_{top}^{3}",     "Events", output, 2, true );
   
    drawStack( f, "Evt_Top_Hadronic_Chi2", "#chi^{2}",           "Events", output, 10, true );
    drawStack( f, "Evt_Top_Hadronic_Chi2", "#chi^{2}",           "Events", output, 10, false );
    drawStack( f, "Evt_Top_Hadronic_Mass", "M_{bqq}^{top}[GeV]", "Events", output, 10, true );
    drawStack( f, "Evt_Top_Hadronic_Mass", "M_{bqq}^{top}[GeV]", "Events", output, 10, false );

    drawStack( f, "Evt_CutFlow_Mu", "", "Events", output, 1, true );
    drawStack( f, "Evt_CutFlow_El", "", "Events", output, 1, true );
    drawStack( f, "Evt_CutFlow_Mu", "", "Events", output, 1, false );
    drawStack( f, "Evt_CutFlow_El", "", "Events", output, 1, false );

    getCutFlowNum( f, "BkgMC__Evt_CutFlow_Mu", true, output);
    getCutFlowNum( f, "BkgMC__Evt_CutFlow_El", true, output);
    getCutFlowNum( f, "TTJets_SemiLeptMGDecays__Evt_CutFlow_Mu",    true, output);
    getCutFlowNum( f, "TTJets_SemiLeptMGDecays__Evt_CutFlow_El",    true, output);
    getCutFlowNum( f, "TTJets_NonSemiLeptMGDecays__Evt_CutFlow_Mu", true, output);
    getCutFlowNum( f, "TTJets_NonSemiLeptMGDecays__Evt_CutFlow_El", true, output);
    getCutFlowNum( f, "SingleT__Evt_CutFlow_Mu",                    true, output);
    getCutFlowNum( f, "SingleT__Evt_CutFlow_El",                    true, output);
    getCutFlowNum( f, "DiBoson__Evt_CutFlow_Mu",                    true, output);
    getCutFlowNum( f, "DiBoson__Evt_CutFlow_El",                    true, output);
    getCutFlowNum( f, "WJetsToLNu__Evt_CutFlow_Mu",                 true, output);
    getCutFlowNum( f, "WJetsToLNu__Evt_CutFlow_El",                 true, output);
    getCutFlowNum( f, "DYJetsToLL__Evt_CutFlow_Mu",                 true, output);
    getCutFlowNum( f, "DYJetsToLL__Evt_CutFlow_El",                 true, output);
}
