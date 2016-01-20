#include "../caculate.C"
char* getNumber( TH1D* h_asym )
{
    char text[128];
    float pN = h_asym->GetBinContent(2);
    float nN = h_asym->GetBinContent(1);
    float pE = h_asym->GetBinError(2);
    float nE = h_asym->GetBinError(1);
    float acp  = caculateACP(h_asym)*100;
    float acpE = caculateACPerrorWrt(h_asym)*100;
    sprintf( text, "%9s %9.2f %9.2f %9.0f %9.0f %9.0f %9.0f", " ", acp, acpE, pN, nN, pE, nE);
    std::cout<<text<<std::endl;
    return text; 
}
char* getChanges( TH1D* h_change )
{
    char text[128];
    float ppN = h_change->GetBinContent(1);
    float nnN = h_change->GetBinContent(2);
    float npN = h_change->GetBinContent(3);
    float pnN = h_change->GetBinContent(4);
    float ppE = h_change->GetBinError(1);
    float nnE = h_change->GetBinError(2);
    float npE = h_change->GetBinError(3);
    float pnE = h_change->GetBinError(4);
    sprintf( text, "%9.1f(%4.2f) %9.1f(%4.2f) %9.1f(%4.2f) %9.1f(%4.2f)", ppN, ppE, nnN, nnE, npN, npE, pnN, pnE);
    std::cout<<text<<std::endl;
    return text; 
}
void checkGenACP( TFile* f, std::string Obs="O2", bool chi2Cut=1, std::string output=".")
{
    FILE* outTxt;
    std::string dirName = "CheckEventsOfLepJets/";
    TH1D *h_reco_asym, *h_reco_asym_el, *h_reco_asym_mu;
    TH1D *h_gen_asym,  *h_gen_asym_el,  *h_gen_asym_mu;
    TH1D *h_change, *h_change_el, *h_change_mu;
    if( chi2Cut )
    {
        outTxt = fopen((output+"/CheckGenObsChi2_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+"EvtChi2_"+Obs+"Asym").c_str());  
        h_reco_asym_el = (TH1D*)f->Get((dirName+"EvtChi2_"+Obs+"Asym_El").c_str()); 
        h_reco_asym_mu = (TH1D*)f->Get((dirName+"EvtChi2_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+"GenChi2_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+"GenChi2_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+"GenChi2_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+"EvtChi2_Change"+Obs).c_str());          
        h_change_mu    = (TH1D*)f->Get((dirName+"EvtChi2_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+"EvtChi2_Change"+Obs+"_El").c_str());
    }else{
        outTxt = fopen((output+"/CheckGenObs_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+"Evt_"+Obs+"Asym").c_str());
        h_reco_asym_el = (TH1D*)f->Get((dirName+"Evt_"+Obs+"Asym_El").c_str());
        h_reco_asym_mu = (TH1D*)f->Get((dirName+"Evt_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+"Gen_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+"Gen_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+"Gen_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+"Evt_Change"+Obs).c_str());
        h_change_mu    = (TH1D*)f->Get((dirName+"Evt_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+"Evt_Change"+Obs+"_El").c_str());
    }

    // Get ACP number
    char text[128];
    sprintf( text, "Electron channel"); 
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "%9s %6s(%%) %9s %9s %9s %9s %9s", "RECO", "Acp", "err(Acp)", "Positive", "Negative", "error(P)", "error(n)"); 
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_reco_asym_el) );
    fprintf( outTxt, "\n" );
    sprintf( text, "%9s ", "GEN");  
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym_el));
    fprintf( outTxt, "\n" );

    sprintf( text, "Muon channel" );
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "%9s %6s(%%) %9s %9s %9s %9s %9s\n", "RECO", "Acp", "err(Acp)", "Positive", "Negative", "error(P)", "error(n)");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_reco_asym_mu));
    fprintf( outTxt, "\n" );
    sprintf( text, "%9s ", "GEN");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym_mu));
    fprintf( outTxt, "\n" );

    sprintf( text, "Combined channel" );
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "%9s %6s(%%) %9s %9s %9s %9s %9s\n", "RECO", "Acp", "err(Acp)", "Positive", "Negative", "error(P)", "error(n)");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_reco_asym));
    fprintf( outTxt, "\n" );
    sprintf( text, "%9s ", "GEN");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym));
    fprintf( outTxt, "\n" );
    fprintf( outTxt, "\n" );
    fprintf( outTxt, "\n" );

    // Get changes
    h_change->Scale(100/(h_change->Integral()));
    h_change_el->Scale(100/(h_change_el->Integral()));
    h_change_mu->Scale(100/(h_change_mu->Integral()));

    sprintf( text, "Event changes from GEN -> RECO"); 
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "%9s(%4s) %9s(%4s) %9s(%4s) %9s(%4s)", "p->p", "err.", "n->n", "err.", "n->p", "err.", "p->n", "err.");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "Electron channel"); 
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getChanges(h_change_el) );
    fprintf( outTxt, "\n" );
    sprintf( text, "Muon channel"); 
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getChanges(h_change_mu) );
    fprintf( outTxt, "\n" );
    sprintf( text, "Combined channel"); 
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getChanges(h_change) );
    fprintf( outTxt, "\n" );
    fclose( outTxt );
}
