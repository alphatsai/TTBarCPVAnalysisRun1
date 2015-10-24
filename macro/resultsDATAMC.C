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
#include "draw2D.C"
using namespace std;
void resultsDATAMC(TFile* f, string output="."){

    //drawObservable( f, "MC__Evt_O7Asym",      output, "O_{7}^{e+#mu}", 	"Combined", 0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym",      output, "O_{2}^{e+#mu}", 	"Combined", 0.05, 1);
    //drawObservable( f, "MC__Evt_O7Asym_Mu",   output, "O_{7}^{#mu}", 	"Muon", 	0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym_Mu",   output, "O_{2}^{#mu}", 	"Muon", 	0.05, 1);
    //drawObservable( f, "MC__Evt_O7Asym_El",   output, "O_{7}^{e}", 		"Electron", 0.05, 0);
    //drawObservable( f, "MC__Evt_O2Asym_El",   output, "O_{2}^{e}", 		"Electron", 0.05, 1);

    //drawObservableDist( f, output, "MC__Evt_O7", 	"O_{7}/M_{top}^{3}", "LepJets", 1);
    //drawObservableDist( f, output, "MC__Evt_O2", 	"O_{2}/M_{top}^{3}", "LepJets", 1);
    drawACP( f, 1, "", "MC__Evt",    "O2", output, "ACP", "O_{2}");
    drawACP( f, 1, "", "MC__Evt",    "O3", output, "ACP", "O_{3}");
    drawACP( f, 1, "", "MC__Evt",    "O4", output, "ACP", "O_{4}");
    drawACP( f, 1, "", "MC__Evt",    "O7", output, "ACP", "O_{7}");
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
    // Linear
    drawStackWithData(_file0, "Evt_NSelJets_El",          "N(selected jets) in e-channel",                 "Events", output, 1 );  
    drawStackWithData(_file0, "Evt_NSelJets_Mu",          "N(selected jets) in #mu-channel",               "Events", output, 1 );  
    drawStackWithData(_file0, "Evt_NBJets_El",            "N(CSVM b-tagged jets) in e-channel",            "Events", output, 1 );  
    drawStackWithData(_file0, "Evt_NBJets_Mu",            "N(CSVM b-tagged jets) in #mu-channel",          "Events", output, 1 );  
    drawStackWithData(_file0, "Evt_HardJet_Pt_El",        "p_{T}(hard jet) in e-channel",                  "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardJet_Pt_Mu",        "p_{T}(hard jet) in #mu-channel",                "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardJet_E_El",         "Energy(hard jet) in e-channel",                 "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardJet_E_Mu",         "Energy(hard jet) in #mu-channel",               "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardJet_Eta_El",       "#eta(hard jet) in e-channel",                   "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_HardJet_Eta_Mu",       "#eta(hard jet) in #mu-channel",                 "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_Pt_El",   "p_{T}(hard non b-jet) in e-channel",            "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardNonBJet1_Pt_Mu",   "p_{T}(hard non b-jet) in #mu-channel",          "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_E_El",    "Energy(hard non b-jet) in e-channel",           "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardNonBJet1_E_Mu",    "Energy(hard non b-jet) in #mu-channel",         "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_Eta_El",  "#eta(hard non b-jet) in e-channel",             "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_HardNonBJet1_Eta_Mu",  "#eta(hard non b-jet) in #mu-channel",           "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_Pt_El",   "p_{T}(second non b-jet) in e-channel",          "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardNonBJet2_Pt_Mu",   "p_{T}(second non b-jet) in #mu-channel",        "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_E_El",    "Energy(second non b-jet) in e-channel",         "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_HardNonBJet2_E_Mu",    "Energy(second non b-jet) in #mu-channel",       "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_Eta_El",  "#eta(second non b-jet) in e-channel",           "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_HardNonBJet2_Eta_Mu",  "#eta(second non b-jet) in #mu-channel",         "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_Pt_El",    "p_{T}(hard non b-jet from t) in e-channel",     "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet from t) in #mu-channel",   "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_E_El",     "Energy(hard non b-jet from t) in e-channel",    "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_E_Mu",     "Energy(hard non b-jet from t) in #mu-channel",  "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_Eta_El",   "#eta(hard non b-jet from t) in e-channel",      "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_Eta_Mu",   "#eta(hard non b-jet from t) in #mu-channel",    "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_Pt_El",    "p_{T}(second non b-jet from t) in e-channel",   "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_Pt_Mu",    "p_{T}(second non b-jet from t) in #mu-channel", "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_E_El",     "Energy(second non b-jet from t) in e-channel",  "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_E_Mu",     "Energy(second non b-jet from t) in #mu-channel","Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_Eta_El",   "#eta(second non b-jet from t) in e-channel",    "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_Eta_Mu",   "#eta(second non b-jet from t) in #mu-channel",  "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_bJet_Pt_El",           "p_{T}(b-jet) in e-channel",                     "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_bJet_Pt_Mu",           "p_{T}(b-jet) in #mu-channel",                   "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_bJet_E_El",            "Energy(b-jet) in e-channel",                    "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_bJet_E_Mu",            "Energy(b-jet) in #mu-channel",                  "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_bJet_Eta_El",          "#eta(b-jet) in e-channel",                      "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_bJet_Eta_Mu",          "#eta(b-jet) in #mu-channel",                    "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_bbarJet_Pt_El",        "p_{T}(#bar{b}-jet) in e-channel",               "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_bbarJet_Pt_Mu",        "p_{T}(#bar{b}-jet) in #mu-channel",             "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_bbarJet_E_El",         "Energy(#bar{b}-jet) in e-channel",              "Events", output, 10 );  
    drawStackWithData(_file0, "Evt_bbarJet_E_Mu",         "Energy(#bar{b}-jet) in #mu-channel",            "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_bbarJet_Eta_El",       "#eta(#bar{b}-jet) in e-channel",                "Events", output, 2  );  
    drawStackWithData(_file0, "Evt_bbarJet_Eta_Mu",       "#eta(#bar{b}-jet) in #mu-channel",              "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_Top_Hadronic_Mass_El", "Inv. M_{top}(jjb) in e-channel",                "Events", output, 10 );
    drawStackWithData(_file0, "Evt_Top_Hadronic_Mass_Mu", "Inv. M_{top}(jjb) in #mu-channel",              "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_Top_Hadronic_Pt_El",   "Inv. p_{T, top}(jjb) in e-channel",             "Events", output, 10 );
    drawStackWithData(_file0, "Evt_Top_Hadronic_Pt_Mu",   "Inv. p_{T, top}(jjb) in #mu-channel",           "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_Top_Hadronic_Chi2_El", "#chi^{2}_{min}(jjb) in e-channel",              "Events", output, 5 );
    drawStackWithData(_file0, "Evt_Top_Hadronic_Chi2_Mu", "#chi^{2}_{min}(jjb) in #mu-channel",            "Events", output, 5 );
    drawStackWithData(_file0, "Evt_Top_Hadronic_Eta_El",  "#eta_{top}(jjb) in e-channel",                  "Events", output, 2  );
    drawStackWithData(_file0, "Evt_Top_Hadronic_Eta_Mu",  "#eta_{top}(jjb) in #mu-channel",                "Events", output, 2  );
    drawStackWithData(_file0, "Evt_isoLep_Pt_Mu",         "p_{T}(isolated #mu)",                           "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_isoLep_Pt_El",         "p_{T}(isolated e)",                             "Events", output, 10 ); 
    drawStackWithData(_file0, "Evt_isoLep_Eta_Mu",        "#eta(isolated #mu)",                            "Events", output, 2  ); 
    drawStackWithData(_file0, "Evt_isoLep_Eta_El",        "#eta(isolated e)",                              "Events", output, 2  );
    drawStackWithData(_file0, "Evt_Ht_El",                "H_{T}(selected jets) in e-channel",             "Events", output, 20 );  
    drawStackWithData(_file0, "Evt_Ht_Mu",                "H_{T}(selected jets) in #mu-channel",           "Events", output, 20 ); 
    // Log
    drawStackWithData(_file0, "Evt_NSelJets_El",          "N(selected jets) in e-channel",                 "Events", output, 1  ,1);  
    drawStackWithData(_file0, "Evt_NSelJets_Mu",          "N(selected jets) in #mu-channel",               "Events", output, 1  ,1);  
    drawStackWithData(_file0, "Evt_NBJets_El",            "N(CSVM b-tagged jets) in e-channel",            "Events", output, 1  ,1);  
    drawStackWithData(_file0, "Evt_NBJets_Mu",            "N(CSVM b-tagged jets) in #mu-channel",          "Events", output, 1  ,1);  
    drawStackWithData(_file0, "Evt_HardJet_Pt_El",        "p_{T}(hard jet) in e-channel",                  "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardJet_Pt_Mu",        "p_{T}(hard jet) in #mu-channel",                "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardJet_E_El",         "Energy(hard jet) in e-channel",                 "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardJet_E_Mu",         "Energy(hard jet) in #mu-channel",               "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardJet_Eta_El",       "#eta(hard jet) in e-channel",                   "Events", output, 2  ,1);  
    drawStackWithData(_file0, "Evt_HardJet_Eta_Mu",       "#eta(hard jet) in #mu-channel",                 "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_Pt_El",    "p_{T}(hard non b-jet) in e-channel",           "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet) in #mu-channel",         "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_E_El",     "Energy(hard non b-jet) in e-channel",          "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet1_E_Mu",     "Energy(hard non b-jet) in #mu-channel",        "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet1_Eta_El",   "#eta(hard b-jet) in e-channel",                "Events", output, 2  ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet1_Eta_Mu",   "#eta(hard b-jet) in #mu-channel",              "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_Pt_El",   "p_{T}(second non b-jet) in e-channel",          "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet2_Pt_Mu",   "p_{T}(second non b-jet) in #mu-channel",        "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_E_El",    "Energy(second non b-jet) in e-channel",         "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet2_E_Mu",    "Energy(second non b-jet) in #mu-channel",       "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_HardNonBJet2_Eta_El",  "#eta(second non b-jet) in e-channel",           "Events", output, 2  ,1);  
    drawStackWithData(_file0, "Evt_HardNonBJet2_Eta_Mu",  "#eta(second non b-jet) in #mu-channel",         "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_Pt_El",    "p_{T}(hard non b-jet from t) in e-channel",     "Events", output, 10,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet from t) in #mu-channel",   "Events", output, 10,1 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_E_El",     "Energy(hard non b-jet from t) in e-channel",    "Events", output, 10,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_E_Mu",     "Energy(hard non b-jet from t) in #mu-channel",  "Events", output, 10,1 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet1_Eta_El",   "#eta(hard non b-jet from t) in e-channel",      "Events", output, 2 ,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet1_Eta_Mu",   "#eta(hard non b-jet from t) in #mu-channel",    "Events", output, 2 ,1 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_Pt_El",    "p_{T}(second non b-jet from t) in e-channel",   "Events", output, 10,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_Pt_Mu",    "p_{T}(second non b-jet from t) in #mu-channel", "Events", output, 10,1 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_E_El",     "Energy(second non b-jet from t) in e-channel",  "Events", output, 10,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_E_Mu",     "Energy(second non b-jet from t) in #mu-channel","Events", output, 10,1 ); 
    drawStackWithData(_file0, "Evt_TopNonBJet2_Eta_El",   "#eta(second non b-jet from t) in e-channel",    "Events", output, 2 ,1 );  
    drawStackWithData(_file0, "Evt_TopNonBJet2_Eta_Mu",   "#eta(second non b-jet from t) in #mu-channel",  "Events", output, 2 ,1 ); 
    drawStackWithData(_file0, "Evt_bJet_Pt_El",           "p_{T}(b-jet) in e-channel",                     "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_bJet_Pt_Mu",           "p_{T}(b-jet) in #mu-channel",                   "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_bJet_E_El",            "Energy(b-jet) in e-channel",                    "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_bJet_E_Mu",            "Energy(b-jet) in #mu-channel",                  "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_bJet_Eta_El",          "#eta(b-jet) in e-channel",                      "Events", output, 2  ,1);  
    drawStackWithData(_file0, "Evt_bJet_Eta_Mu",          "#eta(b-jet) in #mu-channel",                    "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_bbarJet_Pt_El",        "p_{T}(#bar{b}-jet) in e-channel",               "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_bbarJet_Pt_Mu",        "p_{T}(#bar{b}-jet) in #mu-channel",             "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_bbarJet_E_El",         "Energy(#bar{b}-jet) in e-channel",              "Events", output, 10 ,1);  
    drawStackWithData(_file0, "Evt_bbarJet_E_Mu",         "Energy(#bar{b}-jet) in #mu-channel",            "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_bbarJet_Eta_El",       "#eta(#bar{b}-jet) in e-channel",                "Events", output, 2  ,1);  
    drawStackWithData(_file0, "Evt_bbarJet_Eta_Mu",       "#eta(#bar{b}-jet) in #mu-channel",              "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_Top_Hadronic_Chi2_El", "#chi^{2}_{min}(jjb) in e-channel",              "Events", output, 5  ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Chi2_Mu", "#chi^{2}_{min}(jjb) in #mu-channel",            "Events", output, 5  ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Mass_El", "Inv. M_{top}(jjb) in e-channel",                "Events", output, 10 ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Mass_Mu", "Inv. M_{top}(jjb) in #mu-channel",              "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_Top_Hadronic_Eta_El",  "#eta_{top}(jjb) in e-channel",                  "Events", output, 2  ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Eta_Mu",  "#eta_{top}(jjb) in #mu-channel",                "Events", output, 2  ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Pt_El",   "p_{T, top}(jjb) in e-channel",                  "Events", output, 10 ,1);
    drawStackWithData(_file0, "Evt_Top_Hadronic_Pt_Mu",   "p_{T, top}(jjb) in #mu-channel",                "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_isoLep_Pt_Mu",         "p_{T}(isolated #mu)",                           "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_isoLep_Pt_El",         "p_{T}(isolated e)",                             "Events", output, 10 ,1); 
    drawStackWithData(_file0, "Evt_isoLep_Eta_Mu",        "#eta(isolated #mu)",                            "Events", output, 2  ,1); 
    drawStackWithData(_file0, "Evt_isoLep_Eta_El",        "#eta(isolated e)",                              "Events", output, 2  ,1);
    drawStackWithData(_file0, "Evt_Ht_El",                "H_{T}(selected jets) in e-channel",             "Events", output, 20 ,1 );  
    drawStackWithData(_file0, "Evt_Ht_Mu",                "H_{T}(selected jets) in #mu-channel",           "Events", output, 20 ,1 ); 
    // Linear chi2
    drawStackWithData(_file0, "EvtChi2_NSelJets_El",          "N(selected jets) in e-channel",                 "Events", output, 1 );  
    drawStackWithData(_file0, "EvtChi2_NSelJets_Mu",          "N(selected jets) in #mu-channel",               "Events", output, 1 );  
    drawStackWithData(_file0, "EvtChi2_NBJets_El",            "N(CSVM b-tagged jets) in e-channel",            "Events", output, 1 );  
    drawStackWithData(_file0, "EvtChi2_NBJets_Mu",            "N(CSVM b-tagged jets) in #mu-channel",          "Events", output, 1 );  
    drawStackWithData(_file0, "EvtChi2_HardJet_Pt_El",        "p_{T}(hard jet) in e-channel",                  "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardJet_Pt_Mu",        "p_{T}(hard jet) in #mu-channel",                "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardJet_E_El",         "Energy(hard jet) in e-channel",                 "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardJet_E_Mu",         "Energy(hard jet) in #mu-channel",               "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardJet_Eta_El",       "#eta(hard jet) in e-channel",                   "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_HardJet_Eta_Mu",       "#eta(hard jet) in #mu-channel",                 "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Pt_El",   "p_{T}(hard non b-jet) in e-channel",            "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Pt_Mu",   "p_{T}(hard non b-jet) in #mu-channel",          "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_E_El",    "Energy(hard non b-jet) in e-channel",           "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_E_Mu",    "Energy(hard non b-jet) in #mu-channel",         "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Eta_El",  "#eta(hard non b-jet) in e-channel",             "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Eta_Mu",  "#eta(hard non b-jet) in #mu-channel",           "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Pt_El",   "p_{T}(second non b-jet) in e-channel",          "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Pt_Mu",   "p_{T}(second non b-jet) in #mu-channel",        "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_E_El",    "Energy(second non b-jet) in e-channel",         "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_E_Mu",    "Energy(second non b-jet) in #mu-channel",       "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Eta_El",  "#eta(second non b-jet) in e-channel",           "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Eta_Mu",  "#eta(second non b-jet) in #mu-channel",         "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Pt_El",    "p_{T}(hard non b-jet from t) in e-channel",     "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet from t) in #mu-channel",   "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_E_El",     "Energy(hard non b-jet from t) in e-channel",    "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_E_Mu",     "Energy(hard non b-jet from t) in #mu-channel",  "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Eta_El",   "#eta(hard non b-jet from t) in e-channel",      "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Eta_Mu",   "#eta(hard non b-jet from t) in #mu-channel",    "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Pt_El",    "p_{T}(second non b-jet from t) in e-channel",   "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Pt_Mu",    "p_{T}(second non b-jet from t) in #mu-channel", "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_E_El",     "Energy(second non b-jet from t) in e-channel",  "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_E_Mu",     "Energy(second non b-jet from t) in #mu-channel","Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Eta_El",   "#eta(second non b-jet from t) in e-channel",    "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Eta_Mu",   "#eta(second non b-jet from t) in #mu-channel",  "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_bJet_Pt_El",           "p_{T}(b-jet) in e-channel",                     "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_bJet_Pt_Mu",           "p_{T}(b-jet) in #mu-channel",                   "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_bJet_E_El",            "Energy(b-jet) in e-channel",                    "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_bJet_E_Mu",            "Energy(b-jet) in #mu-channel",                  "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_bJet_Eta_El",          "#eta(b-jet) in e-channel",                      "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_bJet_Eta_Mu",          "#eta(b-jet) in #mu-channel",                    "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_Pt_El",        "p_{T}(#bar{b}-jet) in e-channel",               "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_bbarJet_Pt_Mu",        "p_{T}(#bar{b}-jet) in #mu-channel",             "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_E_El",         "Energy(#bar{b}-jet) in e-channel",              "Events", output, 10 );  
    drawStackWithData(_file0, "EvtChi2_bbarJet_E_Mu",         "Energy(#bar{b}-jet) in #mu-channel",            "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_Eta_El",       "#eta(#bar{b}-jet) in e-channel",                "Events", output, 2  );  
    drawStackWithData(_file0, "EvtChi2_bbarJet_Eta_Mu",       "#eta(#bar{b}-jet) in #mu-channel",              "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Mass_El", "Inv. M_{top}(jjb) in e-channel",                "Events", output, 10 );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Mass_Mu", "Inv. M_{top}(jjb) in #mu-channel",              "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Chi2_El", "#chi^{2}_{min}(jjb) in e-channel",              "Events", output, 5 );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Chi2_Mu", "#chi^{2}_{min}(jjb) in #mu-channel",            "Events", output, 5 );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Eta_El",  "#eta_{top}(jjb) in e-channel",                  "Events", output, 2  );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Eta_Mu",  "#eta_{top}(jjb) in #mu-channel",                "Events", output, 2  );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Pt_El",   "p_{T, top}(jjb) in e-channel",             "Events", output, 10 );
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Pt_Mu",   "p_{T, top}(jjb) in #mu-channel",           "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Pt_Mu",         "p_{T}(isolated #mu)",                           "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Pt_El",         "p_{T}(isolated e)",                             "Events", output, 10 ); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Eta_Mu",        "#eta(isolated #mu)",                            "Events", output, 2  ); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Eta_El",        "#eta(isolated e)",                              "Events", output, 2  );
    drawStackWithData(_file0, "EvtChi2_Ht_El",                "H_{T}(selected jets) in e-channel",             "Events", output, 20 );  
    drawStackWithData(_file0, "EvtChi2_Ht_Mu",                "H_{T}(selected jets) in #mu-channel",           "Events", output, 20 ); 
    // Log
    drawStackWithData(_file0, "EvtChi2_NSelJets_El",          "N(selected jets) in e-channel",                 "Events", output, 1  ,1);  
    drawStackWithData(_file0, "EvtChi2_NSelJets_Mu",          "N(selected jets) in #mu-channel",               "Events", output, 1  ,1);  
    drawStackWithData(_file0, "EvtChi2_NBJets_El",            "N(CSVM b-tagged jets) in e-channel",            "Events", output, 1  ,1);  
    drawStackWithData(_file0, "EvtChi2_NBJets_Mu",            "N(CSVM b-tagged jets) in #mu-channel",          "Events", output, 1  ,1);  
    drawStackWithData(_file0, "EvtChi2_HardJet_Pt_El",        "p_{T}(hard jet) in e-channel",                  "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardJet_Pt_Mu",        "p_{T}(hard jet) in #mu-channel",                "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardJet_E_El",         "Energy(hard jet) in e-channel",                 "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardJet_E_Mu",         "Energy(hard jet) in #mu-channel",               "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardJet_Eta_El",       "#eta(hard jet) in e-channel",                   "Events", output, 2  ,1);  
    drawStackWithData(_file0, "EvtChi2_HardJet_Eta_Mu",       "#eta(hard jet) in #mu-channel",                 "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Pt_El",    "p_{T}(hard non b-jet) in e-channel",           "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet) in #mu-channel",         "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_E_El",     "Energy(hard non b-jet) in e-channel",          "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_E_Mu",     "Energy(hard non b-jet) in #mu-channel",        "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Eta_El",   "#eta(hard b-jet) in e-channel",                "Events", output, 2  ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet1_Eta_Mu",   "#eta(hard b-jet) in #mu-channel",              "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Pt_El",   "p_{T}(second non b-jet) in e-channel",          "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Pt_Mu",   "p_{T}(second non b-jet) in #mu-channel",        "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_E_El",    "Energy(second non b-jet) in e-channel",         "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_E_Mu",    "Energy(second non b-jet) in #mu-channel",       "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Eta_El",  "#eta(second non b-jet) in e-channel",           "Events", output, 2  ,1);  
    drawStackWithData(_file0, "EvtChi2_HardNonBJet2_Eta_Mu",  "#eta(second non b-jet) in #mu-channel",         "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Pt_El",    "p_{T}(hard non b-jet from t) in e-channel",     "Events", output, 10,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Pt_Mu",    "p_{T}(hard non b-jet from t) in #mu-channel",   "Events", output, 10,1 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_E_El",     "Energy(hard non b-jet from t) in e-channel",    "Events", output, 10,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_E_Mu",     "Energy(hard non b-jet from t) in #mu-channel",  "Events", output, 10,1 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Eta_El",   "#eta(hard non b-jet from t) in e-channel",      "Events", output, 2 ,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet1_Eta_Mu",   "#eta(hard non b-jet from t) in #mu-channel",    "Events", output, 2 ,1 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Pt_El",    "p_{T}(second non b-jet from t) in e-channel",   "Events", output, 10,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Pt_Mu",    "p_{T}(second non b-jet from t) in #mu-channel", "Events", output, 10,1 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_E_El",     "Energy(second non b-jet from t) in e-channel",  "Events", output, 10,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_E_Mu",     "Energy(second non b-jet from t) in #mu-channel","Events", output, 10,1 ); 
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Eta_El",   "#eta(second non b-jet from t) in e-channel",    "Events", output, 2 ,1 );  
    drawStackWithData(_file0, "EvtChi2_TopNonBJet2_Eta_Mu",   "#eta(second non b-jet from t) in #mu-channel",  "Events", output, 2 ,1 ); 
    drawStackWithData(_file0, "EvtChi2_bJet_Pt_El",           "p_{T}(b-jet) in e-channel",                     "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_bJet_Pt_Mu",           "p_{T}(b-jet) in #mu-channel",                   "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_bJet_E_El",            "Energy(b-jet) in e-channel",                    "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_bJet_E_Mu",            "Energy(b-jet) in #mu-channel",                  "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_bJet_Eta_El",          "#eta(b-jet) in e-channel",                      "Events", output, 2  ,1);  
    drawStackWithData(_file0, "EvtChi2_bJet_Eta_Mu",          "#eta(b-jet) in #mu-channel",                    "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_Pt_El",        "p_{T}(#bar{b}-jet) in e-channel",               "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_bbarJet_Pt_Mu",        "p_{T}(#bar{b}-jet) in #mu-channel",             "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_E_El",         "Energy(#bar{b}-jet) in e-channel",              "Events", output, 10 ,1);  
    drawStackWithData(_file0, "EvtChi2_bbarJet_E_Mu",         "Energy(#bar{b}-jet) in #mu-channel",            "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_bbarJet_Eta_El",       "#eta(#bar{b}-jet) in e-channel",                "Events", output, 2  ,1);  
    drawStackWithData(_file0, "EvtChi2_bbarJet_Eta_Mu",       "#eta(#bar{b}-jet) in #mu-channel",              "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Chi2_El", "#chi^{2}_{min}(jjb) in e-channel",              "Events", output, 5  ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Chi2_Mu", "#chi^{2}_{min}(jjb) in #mu-channel",            "Events", output, 5  ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Mass_El", "Inv. M_{top}(jjb) in e-channel",                "Events", output, 10 ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Mass_Mu", "Inv. M_{top}(jjb) in #mu-channel",              "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Eta_El",  "#eta_{top}(jjb) in e-channel",                  "Events", output, 2  ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Eta_Mu",  "#eta_{top}(jjb) in #mu-channel",                "Events", output, 2  ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Pt_El",   "p_{T, top}(jjb) in e-channel",                  "Events", output, 10 ,1);
    drawStackWithData(_file0, "EvtChi2_Top_Hadronic_Pt_Mu",   "p_{T, top}(jjb) in #mu-channel",                "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Pt_Mu",         "p_{T}(isolated #mu)",                           "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Pt_El",         "p_{T}(isolated e)",                             "Events", output, 10 ,1); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Eta_Mu",        "#eta(isolated #mu)",                            "Events", output, 2  ,1); 
    drawStackWithData(_file0, "EvtChi2_isoLep_Eta_El",        "#eta(isolated e)",                              "Events", output, 2  ,1);
    drawStackWithData(_file0, "EvtChi2_Ht_El",                "H_{T}(selected jets) in e-channel",             "Events", output, 20 ,1 );  
    drawStackWithData(_file0, "EvtChi2_Ht_Mu",                "H_{T}(selected jets) in #mu-channel",           "Events", output, 20 ,1 );

    drawStackWithData( f, "Evt_CutFlow_Mu", "", "Events", output, 1, true );
    drawStackWithData( f, "Evt_CutFlow_El", "", "Events", output, 1, true );
    drawStackWithData( f, "Evt_CutFlow_Mu", "", "Events", output, 1, false );
    drawStackWithData( f, "Evt_CutFlow_El", "", "Events", output, 1, false );

    draw2D(_file0, "MC__TH2_Chi2_vs_TopHadronicMass_El",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 5500   );
    draw2D(_file0, "MC__TH2_Chi2_vs_TopHadronicMass_Mu",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 5500   );
    draw2D(_file0, "DATA_Electron__TH2_Chi2_vs_TopHadronicMass_El", "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 5500, 1);
    draw2D(_file0, "DATA_Muon__TH2_Chi2_vs_TopHadronicMass_Mu",     "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 5500, 1);
    draw2D(_file0, "MC__TH2_Chi2_vs_Ht_El",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20, 1900   );
    draw2D(_file0, "MC__TH2_Chi2_vs_Ht_Mu",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20, 1900   );
    draw2D(_file0, "DATA_Electron__TH2_Chi2_vs_Ht_El",              "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20, 1900, 1);
    draw2D(_file0, "DATA_Muon__TH2_Chi2_vs_Ht_Mu",                  "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20, 1900, 1);
    draw2D(_file0, "MC__TH2_TopHadronicMass_vs_Ht_El",              "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 800 );
    draw2D(_file0, "MC__TH2_TopHadronicMass_vs_Ht_Mu",              "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 800 );
    draw2D(_file0, "DATA_Electron__TH2_TopHadronicMass_vs_Ht_El",   "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 800, 1);
    draw2D(_file0, "DATA_Muon__TH2_TopHadronicMass_vs_Ht_Mu",       "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 800, 1);

    draw2D(_file0, "MC__TH2Chi2_Chi2_vs_TopHadronicMass_El",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 2000   );
    draw2D(_file0, "MC__TH2Chi2_Chi2_vs_TopHadronicMass_Mu",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 2000   );
    draw2D(_file0, "DATA_Electron__TH2Chi2_Chi2_vs_TopHadronicMass_El", "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 2000, 1);
    draw2D(_file0, "DATA_Muon__TH2Chi2_Chi2_vs_TopHadronicMass_Mu",     "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 2000, 1);
    draw2D(_file0, "MC__TH2Chi2_Chi2_vs_Ht_El",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20,  900  );
    draw2D(_file0, "MC__TH2Chi2_Chi2_vs_Ht_Mu",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20,  900  );
    draw2D(_file0, "DATA_Electron__TH2Chi2_Chi2_vs_Ht_El",              "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20, 900, 1);
    draw2D(_file0, "DATA_Muon__TH2Chi2_Chi2_vs_Ht_Mu",                  "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20, 900, 1);
    draw2D(_file0, "MC__TH2Chi2_TopHadronicMass_vs_Ht_El",              "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 300 );
    draw2D(_file0, "MC__TH2Chi2_TopHadronicMass_vs_Ht_Mu",              "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 300 );
    draw2D(_file0, "DATA_Electron__TH2Chi2_TopHadronicMass_vs_Ht_El",   "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 300, 1);
    draw2D(_file0, "DATA_Muon__TH2Chi2_TopHadronicMass_vs_Ht_Mu",       "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 300,  1);

    //draw2D(_file0, "MC__TH2_Chi2_vs_TopHadronicMass_El",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 2000   );
    //draw2D(_file0, "MC__TH2_Chi2_vs_TopHadronicMass_Mu",            "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 2000   );
    //draw2D(_file0, "DATA_Electron__TH2_Chi2_vs_TopHadronicMass_El", "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in e-channel",                         5, 10, 2000, 1);
    //draw2D(_file0, "DATA_Muon__TH2_Chi2_vs_TopHadronicMass_Mu",     "", output, "#chi^{2}_{min}", "Inv. M_{top}(jjb) in #mu-channel",                       5, 10, 2000, 1);
    //draw2D(_file0, "MC__TH2_Chi2_vs_Ht_El",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20,  900  );
    //draw2D(_file0, "MC__TH2_Chi2_vs_Ht_Mu",                         "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20,  900  );
    //draw2D(_file0, "DATA_Electron__TH2_Chi2_vs_Ht_El",              "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in e-channel",                      5, 20, 900, 1);
    //draw2D(_file0, "DATA_Muon__TH2_Chi2_vs_Ht_Mu",                  "", output, "#chi^{2}_{min}", "H_{T}(selected jets) in #mu-channel",                    5, 20, 900, 1);
    //draw2D(_file0, "MC__TH2_TopHadronicMass_vs_Ht_El",              "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 300 );
    //draw2D(_file0, "MC__TH2_TopHadronicMass_vs_Ht_Mu",              "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 300 );
    //draw2D(_file0, "DATA_Electron__TH2_TopHadronicMass_vs_Ht_El",   "", output, "Inv. M_{top}(jjb) in e-channel",   "H_{T}(selected jets) in e-channel",   10, 20, 300, 1);
    //draw2D(_file0, "DATA_Muon__TH2_TopHadronicMass_vs_Ht_Mu",       "", output, "Inv. M_{top}(jjb) in #mu-channel", "H_{T}(selected jets) in #mu-channel", 10, 20, 300,  1);

    getCutFlowNum( f, "DATA_Muon__Evt_CutFlow_Mu", true, output);
    getCutFlowNum( f, "DATA_Electron__Evt_CutFlow_El", true, output);
    getCutFlowNum( f, "MC__Evt_CutFlow_Mu", true, output);
    getCutFlowNum( f, "MC__Evt_CutFlow_El", true, output);
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
