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
void getResultsLepJet(TFile* f, string output="."){

	drawObservable( f, "MC__Evt_O7Asym", output, "O_{7}", "Combined", 0.1, 0);
	drawObservable( f, "MC__Evt_O2Asym", output, "O_{2}", "Combined", 0.3, 1);
	drawObservable( f, "MC__Evt_O7Asym_Mu", output, "O_{7}", "Muon", 0.1, 0);
	drawObservable( f, "MC__Evt_O2Asym_Mu", output, "O_{2}", "Muon", 0.3, 1);
	drawObservable( f, "MC__Evt_O7Asym_El", output, "O_{7}", "Electron", 0.1, 0);
	drawObservable( f, "MC__Evt_O2Asym_El", output, "O_{2}", "Electron", 0.3, 1);
	drawObservable( f, "MC__Evt2b_O7Asym", output, "O_{7}", "Combined", 0.1, 0);
	drawObservable( f, "MC__Evt2b_O2Asym", output, "O_{2}", "Combined", 0.3, 1);
	drawObservable( f, "MC__Evt2b_O7Asym_Mu", output, "O_{7}", "Muon", 0.1, 0);
	drawObservable( f, "MC__Evt2b_O2Asym_Mu", output, "O_{2}", "Muon", 0.3, 1);
	drawObservable( f, "MC__Evt2b_O7Asym_El", output, "O_{7}", "Electron", 0.1, 0);
	drawObservable( f, "MC__Evt2b_O2Asym_El", output, "O_{2}", "Electron", 0.3, 1);

	drawObservableDist( f, output, "MC__Evt_O7", "O_{7}/M_{top}^{3}", "LepJets", 1);
	drawObservableDist( f, output, "MC__Evt_O2", "O_{2}/M_{top}^{3}", "LepJets", 1);
	drawObservableDist( f, output, "MC__Evt2b_O7", "O_{7}/M_{top}^{3}", "LepJets", 1);
	drawObservableDist( f, output, "MC__Evt2b_O2", "O_{2}/M_{top}^{3}", "LepJets", 1);

	drawACPLepJet( f, "MC__Evt", "", output, "ACP", 10, 1 );
	drawACPLepJet( f, "MC__Evt2b", "", output, "ACP", 10, 1 );

	drawStack( f, "Evt2b_O2", "O_{2}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O2", "O_{2}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt2b_O7", "O_{7}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O7", "O_{7}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O2", "O_{2}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O2", "O_{2}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O7", "O_{7}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O7", "O_{7}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt2b_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt2b_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O2_Mu", "O_{2}^{#mu}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O7_Mu", "O_{7}^{#mu}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt2b_O2_El", "O_{2}^{e}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O2_El", "O_{2}^{e}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt2b_O7_El", "O_{7}^{e}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt2b_O7_El", "O_{7}^{e}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O2_El", "O_{2}^{e}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O2_El", "O_{2}^{e}/M_{top}^{3}", "Events", output, 1, true );
	drawStack( f, "Evt_O7_El", "O_{7}^{e}/M_{top}^{3}", "Events", output, 1, false );
	drawStack( f, "Evt_O7_El", "O_{7}^{e}/M_{top}^{3}", "Events", output, 1, true );

	drawStack( f, "Evt_CutFlow_Mu", "", "Events", output, 1, true );
	drawStack( f, "Evt_CutFlow_El", "", "Events", output, 1, true );
	drawStack( f, "Evt_CutFlow_Mu", "", "Events", output, 1, false );
	drawStack( f, "Evt_CutFlow_El", "", "Events", output, 1, false );
	
}
