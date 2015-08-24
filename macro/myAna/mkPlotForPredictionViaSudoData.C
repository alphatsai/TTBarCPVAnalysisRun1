void mkPlotForPredictionViaSudoData( TFile* fin, std::string outpath=".", int iobs=0, int ich=0, bool printTitle=true, const int NPOINTS=12 )
{
    std::string title1, title2; 

    const int NOBS=2; // O2=0, O7=1
    const int OBS_O2=0;
    const int OBS_O7=1;
    int iOBS;
         if( iobs == OBS_O2 ){ iOBS=OBS_O2; title1="O2"; }
    else if( iobs == OBS_O7 ){ iOBS=OBS_O7; title1="O7"; }  
    else{
        printf(">> [ERROR] Please give right observable index:\n");
        printf("           O2:%d, O7:%d\n", OBS_O2, OBS_O7);
        return;
    }

    const int NCH=3;  // Combined=2, Electron=0, Muon=1
    const int CH_Electron=0;
    const int CH_Muon=1;
    const int CH_Combined=2;
    int iCH;
         if( ich == CH_Electron ){ iCH=CH_Electron; title2="Electron_channel"; }
    else if( ich == CH_Muon     ){ iCH=CH_Muon;     title2="Muon_channel";     }  
    else if( ich == CH_Combined ){ iCH=CH_Combined; title2="Combined_channel"; }  
    else{
        printf(">> [ERROR] Please give right channel index:\n");
        printf("           Electron:%d, Muon:%d, Combined:%d\n", CH_Electron, CH_Muon, CH_Combined);
        return;
    }

    const int UP=0;
    const int LOW=1;

    TTree* tin =(TTree*)fin->Get("tree");
    if( NPOINTS != tin->GetEntries() ){ printf(">> [ERROR] NPOINTS should be %d, insted of %d\n", tin->GetEntries(), NPOINTS); return; }

    float predictNonZeroACPMean;
    float predictNonZeroACPUncs[NCH];
    float acpMean[NOBS][NCH]; 
    float acpUncs[NOBS][NCH];
    tin->SetBranchAddress("predictNonZeroACPMean", &predictNonZeroACPMean  );
    tin->SetBranchAddress("predictNonZeroACPUncs", &predictNonZeroACPUncs[0]);
    tin->SetBranchAddress("acpMean", &acpMean[0][0]);
    tin->SetBranchAddress("acpUncs", &acpUncs[0][0]);

    float nonZeroACP[NPOINTS];
    float nonZeroACPUnc[NPOINTS][2];
    float acp[NPOINTS]; 
    float unc_1s[NPOINTS][2];
    float unc_2s[NPOINTS][2];

    for(int idx=0; idx<NPOINTS; idx++){
        tin->GetEntry(idx);

        nonZeroACP[idx]    = predictNonZeroACPMean*100;
        nonZeroACPUnc[idx][UP]  = nonZeroACP[idx] + predictNonZeroACPUncs[iCH]*100;
        nonZeroACPUnc[idx][LOW] = nonZeroACP[idx] - predictNonZeroACPUncs[iCH]*100;

        acp[idx] = acpMean[iOBS][iCH]*100;

        unc_1s[idx][UP]  = acp[idx] + acpUncs[iOBS][iCH]*100;
        unc_2s[idx][UP]  = acp[idx] + acpUncs[iOBS][iCH]*100*2;
        unc_1s[idx][LOW] = acp[idx] - acpUncs[iOBS][iCH]*100;
        unc_2s[idx][LOW] = acp[idx] - acpUncs[iOBS][iCH]*100*2;
        //printf("%d\n", idx);
        //printf("+ 1s %.2f, 2s %.2f\n", unc_1s[idx][UP], unc_2s[idx][UP]);
        //printf("- 1s %.2f, 2s %.2f\n", unc_1s[idx][LOW], unc_2s[idx][LOW]);
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH2F *frame = new TH2F( ("frame"+title1+title2).c_str(), "frame", NPOINTS, nonZeroACP[0], nonZeroACP[NPOINTS-1], NPOINTS, nonZeroACP[0], nonZeroACP[NPOINTS-1]);

    frame->SetStats(kFALSE);
    frame->SetXTitle("Ideal non-Zero ACP [%]");
    frame->SetYTitle("Predicted non-Zero ACP [%]");
    if( printTitle ) frame->SetTitle((title1+"_"+title2).c_str());
    frame->Draw();

    c1->SetGridx();
    c1->SetGridy();

    TPolyLine *pl_2s = new TPolyLine(NPOINTS*2);
    TPolyLine *pl_1s = new TPolyLine(NPOINTS*2);
    TPolyLine *o_1s  = new TPolyLine(NPOINTS*2);

    TGraph *pl_med = new TGraph(NPOINTS);
    TGraph *pl_obs = new TGraph(NPOINTS);

    for(int i=0; i<NPOINTS; i++) {
        pl_2s->SetNextPoint( nonZeroACP[i], unc_2s[i][UP] );
        pl_1s->SetNextPoint( nonZeroACP[i], unc_1s[i][UP] );
        o_1s ->SetNextPoint( nonZeroACP[i], nonZeroACPUnc[i][UP]);

        pl_med->SetPoint( i, nonZeroACP[i], nonZeroACP[i]);
        pl_obs->SetPoint( i, nonZeroACP[i], acp[i]       );

    }
    for(int i=NPOINTS-1; i>=0; i--) {
        pl_2s->SetNextPoint(  nonZeroACP[i], unc_2s[i][LOW] );
        pl_1s->SetNextPoint(  nonZeroACP[i], unc_1s[i][LOW] );
        o_1s ->SetNextPoint(  nonZeroACP[i], nonZeroACPUnc[i][LOW]);
    }

    pl_2s->SetFillColor(46);
    pl_2s->Draw("f");
    pl_1s->SetFillColor(kGreen);
    pl_1s->Draw("f");
    o_1s->SetFillStyle(3244);
    o_1s->SetFillColor(13);
    o_1s->Draw("f");
    pl_med->SetLineStyle(7);
    pl_med->SetLineWidth(2);
    pl_med->Draw("l");
    pl_obs->Draw("*lsame");

    c1->SaveAs((outpath+"/PredictedACPScan_"+title1+"_"+title2+".pdf").c_str());
}
