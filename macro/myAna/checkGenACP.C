#include "TGraphErrors.h"
#include "../caculate.C"
#include "../CMS_lumi.C"
#include "../help.C"
const float nomWrt=0.08714819442;
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
char* getAcpChangeRate( TH1D* h_asym1, TH1D* h_asym2, std::string title="RECO/GEN" )
{
    char text[128];
    float acp1  = caculateACP(h_asym1)*100;
    float acpE1 = caculateACPerrorWrt(h_asym1)*100;
    float acp2  = caculateACP(h_asym2)*100;
    float acpE2 = caculateACPerrorWrt(h_asym2)*100;
    float rate=acp1/acp2;
    float rateE=rate*sqrt(acpE1*acpE1/acp1/acp1+acpE2*acpE2/acp2/acp2);
    //sprintf( text, "%12s %.2f (%.2f)", title.c_str(), 1-rate, rateE );
    sprintf( text, "%12s %.2f (%.2f)", title.c_str(), rate, rateE );
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
    float pD  = (ppN-pnN)/(ppN+pnN);
    float nD  = (nnN-npN)/(nnN+npN);
    float Df  = (pD+nD)/2;
    sprintf( text, "%9.1f(%4.2f) %9.1f(%4.2f) %9.1f(%4.2f) %9.1f(%4.2f) %9.2f", ppN, ppE, nnN, nnE, npN, npE, pnN, pnE, Df);
    std::cout<<text<<std::endl;
    return text; 
}
float getSumChangeValue( TH1D* h_change, int idx=0 )
{
    float value=0;
    if( idx == 0 ) value = h_change->GetBinContent(1) + h_change->GetBinContent(2);
    else           value = h_change->GetBinContent(3) + h_change->GetBinContent(4);
    return value;
}
float getDilution( float sameN, float exN )
{
    return (sameN-exN)/(sameN+exN);
}
void getChagnesSyst( FILE* outTxt, std::string dirPath, std::string* systName, int systNi, std::string Obs="O2", std::string nreco="EvtChi2" )
{
    const int systN = systNi;
    std::string tune[2]={"up","down"};
    std::string chName[3]={"","_El","_Mu"};
    std::string chName_[3]={"Combined","Electron","Muon"};
    TH1D *h[3];
    TH1D *h_syst[3][systN][2];
    float v[3][3];
    float v_syst[3][systN][2][3];

    std::string hdirName = "CheckEventsOfLepJets/";
    std::string rootname="TTJets_SemiLeptMGDecays.root"; 
    std::string nominalDir="nominal"; 

    TFile* f = new TFile((dirPath+"/"+nominalDir+"/"+rootname).c_str());
    for( int c=0; c<3; c++ )
    { 
        h[c] = (TH1D*)((TH1D*)f->Get((hdirName+nreco+"_Change"+Obs+chName[c]).c_str()))->Clone();
        h[c]->Scale(100/(h[c]->Integral())); 
        v[c][0] = getSumChangeValue( h[c], 0 );
        v[c][1] = getSumChangeValue( h[c], 1 );
        v[c][2] = getDilution( v[c][0], v[c][1] );
    } 

    int iTopMass=-1;
    for( int i=0; i<systN; i++ )
    {
        if( systName[i].find("topMass")!=std::string::npos ) iTopMass=i; 
        for( int j=0; j<2; j++)
        {
            TFile* fs = new TFile((dirPath+"/"+systName[i]+tune[j]+"/"+rootname).c_str());
            for( int c=0; c<3; c++ )
            { 
                h_syst[c][i][j] = (TH1D*)((TH1D*)fs->Get((hdirName+nreco+"_Change"+Obs+chName[c]).c_str()))->Clone(); 
                h_syst[c][i][j]->Scale(100/(h_syst[c][i][j]->Integral()));
                float rateS = getSumChangeValue( h_syst[c][i][j], 0 );
                float rateE = getSumChangeValue( h_syst[c][i][j], 1 );
                if( iTopMass == i ) // rescale top mass uncertaity to 1 GeV
                {
                    rateS = (rateS+5*v[c][0])/6;
                    rateE = (rateE+5*v[c][1])/6;
                }
                v_syst[c][i][j][2] = getDilution( rateS, rateE );
                v_syst[c][i][j][2] = v_syst[c][i][j][2] - v[c][2];
                v_syst[c][i][j][0] = rateS - v[c][0];
                v_syst[c][i][j][1] = rateE - v[c][1];
            }          
        }
    }
   
    float sumw2Dsyst[3][2]; 
    // Print texts
    fprintf( outTxt, "### Observable %s\n", Obs.c_str() );
    for( int c=0; c<3; c++ )
    { 
        sumw2Dsyst[c][0] = 0;
        sumw2Dsyst[c][1] = 0;
        for( int i=0; i<systN; i++ )
        {
            if( v_syst[c][i][0][2]*v_syst[c][i][1][2] < 0 )
            {
                if( v_syst[c][i][0][2] > v_syst[c][i][1][2] ){
                    sumw2Dsyst[c][0] += v_syst[c][i][0][2]*v_syst[c][i][0][2];
                    sumw2Dsyst[c][1] += v_syst[c][i][1][2]*v_syst[c][i][1][2];
                }else{
                    sumw2Dsyst[c][0] += v_syst[c][i][1][2]*v_syst[c][i][1][2];
                    sumw2Dsyst[c][1] += v_syst[c][i][0][2]*v_syst[c][i][0][2];
                }
            }
            else
            {
                if( v_syst[c][i][0][2] > 0 ){
                    if( v_syst[c][i][0][2] > v_syst[c][i][1][2] ){
                        sumw2Dsyst[c][0] += v_syst[c][i][0][2]*v_syst[c][i][0][2];
                    }else{
                        sumw2Dsyst[c][0] += v_syst[c][i][1][2]*v_syst[c][i][1][2];
                    }
                }else{
                    if( v_syst[c][i][0][2] < v_syst[c][i][1][2] ){
                        sumw2Dsyst[c][1] += v_syst[c][i][0][2]*v_syst[c][i][0][2];
                    }else{
                        sumw2Dsyst[c][1] += v_syst[c][i][1][2]*v_syst[c][i][1][2];
                    }
                }

            }
        }

        fprintf( outTxt, "## %s channel\n", chName_[c].c_str() );
        fprintf( outTxt, "Nominal R(no)/R(ex)/D : %.2f/%.2f/%.3f\n", v[c][0], v[c][1], v[c][2] );
        fprintf( outTxt, "SystUnc ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7s ", systName[i].c_str());
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dR(no)+ ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.2f ", v_syst[c][i][0][0] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dR(no)- ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.2f ", v_syst[c][i][1][0] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dR(ex)+ ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.2f ", v_syst[c][i][0][1] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dR(ex)- ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.2f ", v_syst[c][i][1][1] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dD+     ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.3f ", v_syst[c][i][0][2] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "dD-     ");
        for( int i=0; i<systN; i++ ) fprintf( outTxt, "%7.3f ", v_syst[c][i][1][2] );
        fprintf( outTxt, "\n" );
        fprintf( outTxt, "Total D +/- : %.3f/%.3f\n", sqrt(sumw2Dsyst[c][0]), -1*sqrt(sumw2Dsyst[c][1]));
        fprintf( outTxt, "\n" );
    }
    fprintf( outTxt, "\n" );
}

void checkGenACP( TFile* f, std::string Obs="O2", bool chi2Cut=1, std::string output=".", std::string nreco="Evt", std::string ngen="Gen" )
{
    cout<<"[INFO] checkGenACP "<<Obs<<endl;
    FILE* outTxt;
    std::string dirName = "CheckEventsOfLepJets/";
    TH1D *h_reco_asym, *h_reco_asym_el, *h_reco_asym_mu;
    TH1D *h_gen_asym,  *h_gen_asym_el,  *h_gen_asym_mu;
    TH1D *h_change, *h_change_el, *h_change_mu;
    if( chi2Cut )
    {
        outTxt = fopen((output+"/CheckGenObsChi2_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+nreco+"Chi2_"+Obs+"Asym").c_str());  
        h_reco_asym_el = (TH1D*)f->Get((dirName+nreco+"Chi2_"+Obs+"Asym_El").c_str()); 
        h_reco_asym_mu = (TH1D*)f->Get((dirName+nreco+"Chi2_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+ngen+"Chi2_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+ngen+"Chi2_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+ngen+"Chi2_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+nreco+"Chi2_Change"+Obs).c_str());          
        h_change_mu    = (TH1D*)f->Get((dirName+nreco+"Chi2_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+nreco+"Chi2_Change"+Obs+"_El").c_str());
    }else{
        outTxt = fopen((output+"/CheckGenObs_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+nreco+"_"+Obs+"Asym").c_str());
        h_reco_asym_el = (TH1D*)f->Get((dirName+nreco+"_"+Obs+"Asym_El").c_str());
        h_reco_asym_mu = (TH1D*)f->Get((dirName+nreco+"_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+ngen+"_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+ngen+"_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+ngen+"_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+nreco+"_Change"+Obs).c_str());
        h_change_mu    = (TH1D*)f->Get((dirName+nreco+"_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+nreco+"_Change"+Obs+"_El").c_str());
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
    fprintf( outTxt, getAcpChangeRate(h_reco_asym_el,h_gen_asym_el));
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
    fprintf( outTxt, getAcpChangeRate(h_reco_asym_mu,h_gen_asym_mu));
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
    fprintf( outTxt, getAcpChangeRate(h_reco_asym,h_gen_asym));
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
    sprintf( text, "%9s(%4s) %9s(%4s) %9s(%4s) %9s(%4s) %9s", "p->p", "err.", "n->n", "err.", "n->p", "err.", "p->n", "err.", "D-factor");
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

void checkParACP( TFile* f, std::string Obs="O2", bool chi2Cut=1, std::string output=".")
{
    cout<<"[INFO] checkParACP "<<Obs<<endl;
    FILE* outTxt;
    std::string dirName = "CheckEventsOfLepJets/";
    TH1D *h_reco_asym, *h_reco_asym_el, *h_reco_asym_mu;
    TH1D *h_gen_asym,  *h_gen_asym_el,  *h_gen_asym_mu;
    TH1D *h_change, *h_change_el, *h_change_mu;
    if( chi2Cut )
    {
        outTxt = fopen((output+"/CheckParObsChi2_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+"EvtBBChi2_"+Obs+"Asym").c_str());  
        h_reco_asym_el = (TH1D*)f->Get((dirName+"EvtBBChi2_"+Obs+"Asym_El").c_str()); 
        h_reco_asym_mu = (TH1D*)f->Get((dirName+"EvtBBChi2_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+"ParChi2_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+"ParChi2_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+"ParChi2_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+"ParChi2_Change"+Obs).c_str());          
        h_change_mu    = (TH1D*)f->Get((dirName+"ParChi2_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+"ParChi2_Change"+Obs+"_El").c_str());
    }else{
        outTxt = fopen((output+"/CheckParObs_"+Obs+".txt").c_str(),"w");
        h_reco_asym    = (TH1D*)f->Get((dirName+"EvtBB_"+Obs+"Asym").c_str());
        h_reco_asym_el = (TH1D*)f->Get((dirName+"EvtBB_"+Obs+"Asym_El").c_str());
        h_reco_asym_mu = (TH1D*)f->Get((dirName+"EvtBB_"+Obs+"Asym_Mu").c_str());
        h_gen_asym     = (TH1D*)f->Get((dirName+"Par_"+Obs+"Asym").c_str());
        h_gen_asym_el  = (TH1D*)f->Get((dirName+"Par_"+Obs+"Asym_El").c_str());
        h_gen_asym_mu  = (TH1D*)f->Get((dirName+"Par_"+Obs+"Asym_Mu").c_str());
        h_change       = (TH1D*)f->Get((dirName+"Par_Change"+Obs).c_str());
        h_change_mu    = (TH1D*)f->Get((dirName+"Par_Change"+Obs+"_Mu").c_str());
        h_change_el    = (TH1D*)f->Get((dirName+"Par_Change"+Obs+"_El").c_str());
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
    sprintf( text, "%9s ", "Parton");  
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym_el));
    fprintf( outTxt, "\n" );
    fprintf( outTxt, getAcpChangeRate(h_reco_asym_el,h_gen_asym_el, "RECO/PAR"));
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
    sprintf( text, "%9s ", "Parton");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym_mu));
    fprintf( outTxt, "\n" );
    fprintf( outTxt, getAcpChangeRate(h_reco_asym_mu,h_gen_asym_mu, "RECO/PAR"));
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
    sprintf( text, "%9s ", "Parton");
    fprintf( outTxt, text );
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    fprintf( outTxt, getNumber(h_gen_asym));
    fprintf( outTxt, "\n" );
    fprintf( outTxt, getAcpChangeRate(h_reco_asym,h_gen_asym, "RECO/PAR"));
    fprintf( outTxt, "\n" );
    fprintf( outTxt, "\n" );
    fprintf( outTxt, "\n" );

    // Get changes
    h_change->Scale(100/(h_change->Integral()));
    h_change_el->Scale(100/(h_change_el->Integral()));
    h_change_mu->Scale(100/(h_change_mu->Integral()));

    sprintf( text, "Event changes from Parton -> RECO"); 
    fprintf( outTxt, "\n" );
    std::cout<<text<<std::endl;
    sprintf( text, "%9s(%4s) %9s(%4s) %9s(%4s) %9s(%4s) %9s", "p->p", "err.", "n->n", "err.", "n->p", "err.", "p->n", "err.", "D-factor");
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

void mkPlotSlopeACP( std::string inputDir, std::string* inputACPDir, int nGenACP, std::string outpath=".", std::string Obs="O2", std::string nameObs="O_{2}", int ich=1, bool chi2Cut=1, float Dfactor=1 )
{    
    std::string ch;
    std::string channel;
    if( ich == 1 ){
        ch="_El";
        channel="e+jets";
    }else if( ich == 2){
        ch="_Mu";
        channel="#mu+jets";
    }else{
        ch="";
        channel="l+jets";
    }
    TH1D* h_asym[2][nGenACP];
    float acp[3][nGenACP];
    float acp1E[3][nGenACP];
    float acp2E[3][nGenACP];
    for( int i=0; i<nGenACP; i++ )
    {
        TFile* f = new TFile((inputDir+"/"+inputACPDir[i]+"/TTJets_SemiLeptMGDecays.root").c_str());
        if( chi2Cut ){
            h_asym[0][i] = (TH1D*)((TH1D*)f->Get(("CheckEventsOfLepJets/GenChi2_"+Obs+"Asym"+ch).c_str()))->Clone();
            h_asym[1][i] = (TH1D*)((TH1D*)f->Get(("CheckEventsOfLepJets/EvtChi2_"+Obs+"Asym"+ch).c_str()))->Clone();
        }else{
            h_asym[0][i] = (TH1D*)((TH1D*)f->Get(("CheckEventsOfLepJets/Gen_"+Obs+"Asym"+ch).c_str()))->Clone();
            h_asym[1][i] = (TH1D*)((TH1D*)f->Get(("CheckEventsOfLepJets/Evt_"+Obs+"Asym"+ch).c_str()))->Clone();
        }
        h_asym[0][i]->Scale(nomWrt);
        h_asym[1][i]->Scale(nomWrt);
        acp[0][i]  = caculateACP(h_asym[0][i])*100;
        acp[1][i]  = caculateACP(h_asym[1][i])*100;
        acp1E[0][i] = caculateACPerrorWrt(h_asym[0][i])*100;
        acp1E[1][i] = caculateACPerrorWrt(h_asym[1][i])*100;
        //acp1E[0][i] = caculateACPerror(h_asym[0][i])*100;
        //acp1E[1][i] = caculateACPerror(h_asym[1][i])*100;
        //if( i == (nGenACP-1)/2 ){
        //    acp[2][i]   = acp[1][i];
        //    acp1E[2][i] = acp1E[1][i];
        //}else{
            acp[2][i]   = acp[1][i]/Dfactor;
            acp1E[2][i] = acp1E[1][i]/Dfactor;
        //}
        acp2E[0][i] = acp1E[0][i]*2;
        acp2E[1][i] = acp1E[1][i]*2;
        acp2E[2][i] = acp1E[2][i]*2;
        f->Close();
    } 

    TH2D* h = new TH2D( "base", "", 400, -20, 20, 400, -20, 20 );
    TLine* line = new TLine( -20, -20, 20, 20 );
    line->SetLineStyle(3); 

    TGraphErrors* hg0 = new TGraphErrors( nGenACP, acp[0], acp[1], acp1E[0], acp1E[1] );
    TGraphErrors* hgD = new TGraphErrors( nGenACP, acp[0], acp[2], acp1E[0], acp1E[2] );
    TGraphErrors* hg0_2s = new TGraphErrors( nGenACP, acp[0], acp[1], acp2E[0], acp2E[1] );
    TGraphErrors* hgD_2s = new TGraphErrors( nGenACP, acp[0], acp[2], acp2E[0], acp2E[2] );

    //TCanvas *c1 = new TCanvas(("c1"+ch+Obs).c_str(),"", 234,169,800,600); 
    TCanvas *c1 = new TCanvas(("c1"+ch+Obs).c_str(),"", 234,169,W,H); 
    c1->Clear();
    c1->Range(-28.67314,-28.6836,22.8479,24.24942);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);

   c1->SetLeftMargin(0.1683417);
   c1->SetRightMargin(0.05527638);
   c1->SetTopMargin(0.08027923);
   c1->SetBottomMargin(0.1640489);
   //c1->SetLeftMargin(0.1683417);
   //c1->SetRightMargin(0.07286432);
   //c1->SetTopMargin(0.08027923);
   //c1->SetBottomMargin(0.1640489);
    //c1->SetLeftMargin(0.1582915);
    //c1->SetRightMargin(0.07286432);
    //c1->SetTopMargin(0.05043478);
    //c1->SetBottomMargin(0.1634783);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
    c1->SetGridx();
    c1->SetGridy();

    gStyle->SetOptStat(0);

    hgD->SetLineColor(4);
    hgD->SetLineWidth(2);
    hg0->SetLineColor(2);
    hg0->SetLineWidth(2);
    //hg0_2s->SetLineWidth(2);
    //hgD_2s->SetLineWidth(2);

    h->GetXaxis()->SetTitle(("A_{CP}("+nameObs+") in Gen. [%]").c_str());
    h->GetYaxis()->SetTitle(("A'_{CP}("+nameObs+") in Reco. [%]").c_str());
    h->GetXaxis()->SetNdivisions(507);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.07);
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetNdivisions(507);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.07);
    h->GetYaxis()->SetTitleOffset(0.87);
    h->GetYaxis()->SetTitleFont(42);
    h->Draw();

    line->Draw("SAME");
    hgD_2s->Draw("LESAME");
    hgD->Draw("LESAME");
    hg0_2s->Draw("LESAME");
    hg0->Draw("LESAME");

    extraText="Simulation Preliminary";
    writeExtraText=true;
    CMS_lumi( c1, 2, 0, true, false );

    //TLegend *leg = new TLegend(0.1821608,0.6845754,0.5238693,0.9376083,NULL,"brNDC");
    TLegend *leg = new TLegend(0.1896985,0.6579407,0.531407,0.9109948,NULL,"brNDC");
    //leg->SetHeader(("Simulation in "+channel).c_str());
    leg->SetHeader((channel+" channel").c_str());
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05217391);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    char text[200];
    //sprintf(text, "Apply dilution factor D=%.2f", Dfactor );
    sprintf(text, "After dilution factor correction", Dfactor );
    leg->AddEntry( hgD,    text,                 "lpe");
    leg->AddEntry( hg0,    "No dilution factor", "lpe");
    leg->AddEntry( hg0_2s, "#pm2#sigma stat. error",        "lpe");
    leg->Draw();
    if( chi2Cut )
        c1->SaveAs((outpath+"/GenAcpVsRecoAcpChi2_"+Obs+ch+".pdf").c_str());
    else
        c1->SaveAs((outpath+"/GenAcpVsRecoAcp_"+Obs+ch+".pdf").c_str());
    delete h;
}
