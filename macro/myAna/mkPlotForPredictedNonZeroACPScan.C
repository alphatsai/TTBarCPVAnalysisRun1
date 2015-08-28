// TFile should be from the results of mkTreeForPredictedNonZeroACPScan.C or mkTreeForPredictedNonZeroACPScanViaSudoExp.C
void mkPlotForPrediction( TFile* fin, std::string outpath=".", std::string cat="MC", int iobs=0, int ich=0, bool printTitle=true, const int NPOINTS=13 )
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

    float assumedACPMean;
    float assumedACPUncs[NOBS][NCH];
    float acpMean[NOBS][NCH]; 
    float acpUncs[NOBS][NCH];
    tin->SetBranchAddress("assumedACPMean",         &assumedACPMean      );
    tin->SetBranchAddress("assumedACPUncs",         &assumedACPUncs[0][0]);
    tin->SetBranchAddress((cat+"_acpMean").c_str(), &acpMean[0][0]);
    tin->SetBranchAddress((cat+"_acpUncs").c_str(), &acpUncs[0][0]);

    float assumedACP[NPOINTS];
    float assumedACPUnc[NPOINTS][2];
    float acp[NPOINTS]; 
    float unc_1s[NPOINTS][2];
    float unc_2s[NPOINTS][2];

    for(int idx=0; idx<NPOINTS; idx++){
        tin->GetEntry(idx);

        assumedACP[idx]    = assumedACPMean*100;
        assumedACPUnc[idx][UP]  = assumedACP[idx] + assumedACPUncs[iOBS][iCH]*100;
        assumedACPUnc[idx][LOW] = assumedACP[idx] - assumedACPUncs[iOBS][iCH]*100;

        acp[idx] = acpMean[iOBS][iCH]*100;

        unc_1s[idx][UP]  = acp[idx] + acpUncs[iOBS][iCH]*100;
        unc_2s[idx][UP]  = acp[idx] + acpUncs[iOBS][iCH]*100*2;
        unc_1s[idx][LOW] = acp[idx] - acpUncs[iOBS][iCH]*100;
        unc_2s[idx][LOW] = acp[idx] - acpUncs[iOBS][iCH]*100*2;
        printf("%d\n", idx);
        printf("+ 1s %.2f, 2s %.2f\n", unc_1s[idx][UP], unc_2s[idx][UP]);
        printf("- 1s %.2f, 2s %.2f\n", unc_1s[idx][LOW], unc_2s[idx][LOW]);
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH2F *frame = new TH2F( ("frame"+cat+title1+title2).c_str(), "frame", NPOINTS, assumedACP[0], assumedACP[NPOINTS-1], NPOINTS, assumedACP[0], assumedACP[NPOINTS-1]);

    frame->SetStats(kFALSE);
    frame->SetXTitle("Ideal non-Zero ACP [%]");
    frame->SetYTitle("Predicted non-Zero ACP [%]");
    if( printTitle ) frame->SetTitle((cat+" "+title1+"_"+title2).c_str());
    frame->Draw();

    c1->SetGridx();
    c1->SetGridy();

    TPolyLine *pl_2s = new TPolyLine(NPOINTS*2);
    TPolyLine *pl_1s = new TPolyLine(NPOINTS*2);
    TPolyLine *o_1s  = new TPolyLine(NPOINTS*2);

    TGraph *pl_med = new TGraph(NPOINTS);
    TGraph *pl_obs = new TGraph(NPOINTS);

    for(int i=0; i<NPOINTS; i++) {
        pl_2s->SetNextPoint( assumedACP[i], unc_2s[i][UP] );
        pl_1s->SetNextPoint( assumedACP[i], unc_1s[i][UP] );
        o_1s ->SetNextPoint( assumedACP[i], assumedACPUnc[i][UP]);

        pl_med->SetPoint( i, assumedACP[i], assumedACP[i]);
        pl_obs->SetPoint( i, assumedACP[i], acp[i]       );

    }
    for(int i=NPOINTS-1; i>=0; i--) {
        pl_2s->SetNextPoint(  assumedACP[i], unc_2s[i][LOW] );
        pl_1s->SetNextPoint(  assumedACP[i], unc_1s[i][LOW] );
        o_1s ->SetNextPoint(  assumedACP[i], assumedACPUnc[i][LOW]);
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

    c1->SaveAs((outpath+"/PredictedACPScan_"+cat+"_"+title1+"_"+title2+".pdf").c_str());
}

void mkPull( TFile* fin, std::string outpath=".", std::string cat="SudoExp", int iobs=0, int ich=0, bool printTitle=true )
{
    std::string title1, title2; 

    const int NENTRY=13;
    float  assumedACP[NENTRY]  ={  -30,   -25,   -20,   -15,   -10,   -5,   0,    5,    10,    15,    20,    25,    30 };
    string assumedACP_s[NENTRY]={ "m30", "m25", "m20", "m15", "m10", "m5", "0", "p5", "p10", "p15", "p20", "p25", "p30"};
 
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
         if( ich == CH_Electron ){ iCH=CH_Electron; title2="_El"; }
    else if( ich == CH_Muon     ){ iCH=CH_Muon;     title2="_Mu"; }  
    else if( ich == CH_Combined ){ iCH=CH_Combined; title2=""; }  
    else{
        printf(">> [ERROR] Please give right channel index:\n");
        printf("           Electron:%d, Muon:%d, Combined:%d\n", CH_Electron, CH_Muon, CH_Combined);
        return;
    }
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(101);
    gStyle->SetFitFormat("3.3g");
    TCanvas *c1 = new TCanvas("c1pull","c1",800,600);
    TCanvas *c2 = new TCanvas("c1all","c1", 1100,800);
    c2->Divide(5,3);
    TH1D* h[NENTRY];
    for( int i=0; i<NENTRY; i++)
    {
        c1->Clear();
        c1->cd();

        h[i] = (TH1D*)fin->Get((cat+"_"+assumedACP_s[i]+"ACPmean_"+title1+title2).c_str());
        float max = assumedACP[i]+10;
        float min = assumedACP[i]-10;
        float maxy = h[i]->GetMaximum()*1.5;
        h[i]->GetXaxis()->SetRangeUser(min, max);
        h[i]->GetXaxis()->SetNdivisions(509);
        h[i]->GetXaxis()->SetLabelFont(42);
        h[i]->GetXaxis()->SetLabelSize(0.05);
        h[i]->GetXaxis()->SetTitleSize(0.05);
        h[i]->GetXaxis()->SetTitleOffset(0.9);
        
        h[i]->GetYaxis()->SetRangeUser( 0., maxy);
        h[i]->GetYaxis()->SetLabelFont(42);
        h[i]->GetYaxis()->SetTitleSize(0.05);
        h[i]->GetYaxis()->SetTitleOffset(0.94);
        h[i]->SetXTitle(("ACP("+title1+")"+title2+" [%]").c_str());
        h[i]->SetYTitle("# of Experiment");
        h[i]->Draw("HISTE"); 

        TF1* gaus = new TF1("gaus", "gaus", min, max);
        gaus->SetLineColor(2);
        h[i]->Fit( gaus, "WR"); cout<<endl;
        gaus->Draw("SAME");
        
        TLine* l = new TLine(assumedACP[i], 0, assumedACP[i], maxy);
        l->SetLineWidth(2);
        l->SetLineStyle(7);
        l->Draw("SAME");
    
        c1->SaveAs((outpath+"/Pull_"+cat+"_"+assumedACP_s[i]+"ACP_"+title1+title2+".pdf").c_str());

        c2->cd(i+1);
        h[i]->Draw("HISTE"); 
        gaus->Draw("SAME");
        l->Draw("SAME");
    }
    c2->SaveAs((outpath+"/PullAll_"+cat+"_"+title1+title2+".pdf").c_str());
}
