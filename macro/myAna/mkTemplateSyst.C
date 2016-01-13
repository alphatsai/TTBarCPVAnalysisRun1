#include <string>
#include "../help.C"
#include "templateLFit.C"
#define NPAR 2

std::string fileName="Final_histograms_SemiLeptanic.root";
std::string systDir="Syst";
std::string fName="../results/15Dec_LepJet_MCDATA";
std::string hName="EvtChi2_Top_Leptonic_Mbl";
//std::string hName="EvtChi2_Top_Hadronic_Mass";
//std::string hName="EvtChi2_Ht";
std::string output=fName;
bool doFit=true;
int rebin=10;
//int rebin=20;
int CoCH=0;
int ElCH=1;
int MuCH=2;
int sig=0;
int bkg=1;
void getHistStatUnc( TH1D* h, TH1D* h_u, TH1D* h_d )
{
    int bins=h->GetNbinsX();
    for( int i=1; i<=bins; i++)
    {
        h_u->SetBinContent(i, h->GetBinContent(i)+h->GetBinError(i));
        h_d->SetBinContent(i, h->GetBinContent(i)-h->GetBinError(i));
    }
}
void getHistStatUncNomalized( TH1D* h, TH1D* h_u, TH1D* h_d )
{
    int bins=h->GetNbinsX();
    double sumEvts=0;
    double sumw2=0;
    for( int i=1; i<=bins; i++)
    {
        sumEvts += h->GetBinContent(i);
        sumw2 += h->GetBinError(i)*h->GetBinError(i);
    }
    h_u->Scale((sumEvts+sqrt(sumw2))/h_u->Integral());
    h_d->Scale((sumEvts-sqrt(sumw2))/h_d->Integral());
}
void mkTemplateSyst()
{
    const int nCh=3, nSyst=2+7, nMC=3, nTune=2;
    std::string mcName[nMC]={"TTJets__", "BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__", "MC__"}; // sig=0, bkg=1
    std::string chName[nCh]={"", "_El", "_Mu"};
    std::string tuneName[nTune]={"up", "down"};
    std::string systName[nSyst]={"", "Stat", "PU", "JER", "BTagSF", "TopPT", "elID", "muID", "muISO"};
    
    TFile* fin[nSyst][nTune];
    fin[0][0] = new TFile((fName+"/"+fileName).c_str()); //mean hist
    for( int s=2; s<nSyst; s++ ){
        for( int t=0; t<nTune; t++ ){
            std::string fname=fName+"/"+systDir+"/"+systName[s]+tuneName[t]+"/"+fileName;
            fin[s][t] = new TFile(fname.c_str());
        } 
    }

    TFile* fout  = new TFile((output+"/TemplateSyst_"+hName+".root").c_str(), "RECREATE");

    // Data
    std::cout<<"[INFO] Copying data..."<<std::endl;
    TH1D *h_data[nCh];
    h_data[CoCH] = (TH1D*)((TH1D*)fin[0][0]->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA");
    h_data[ElCH] = (TH1D*)((TH1D*)fin[0][0]->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA_El");
    h_data[MuCH] = (TH1D*)((TH1D*)fin[0][0]->Get(("DATA_Muon__"+hName+"_Mu").c_str()))->Clone("DATA_Mu");
    for( int ch=0; ch<nCh; ch++ )
    { 
        h_data[ch]->Rebin(rebin);  
        fix(h_data[ch]); 
    } 
    h_data[CoCH]->Add(h_data[MuCH]);

    // MC
    std::cout<<"[INFO] Copying MC..."<<std::endl;
    TH1D *h_mc[nMC][nCh][nSyst][nTune];
    for( int mc=0; mc<nMC; mc++ )
    {
        std::string name0;
        if( mc==sig )      name0="SigMC"; 
        else if( mc==bkg ) name0="BkgMC"; 
        else name0="MC"; 
        std::cout<<"       "<<name0<<"..."<<std::endl;
        for( int ch=0; ch<nCh; ch++ )
        {
            std::cout<<"        Ch"<<chName[ch]<<std::endl<<"          ";
            std::string hname = mcName[mc]+hName+chName[ch];
            std::string name1 = name0+chName[ch];
            h_mc[mc][ch][0][0] = (TH1D*)((TH1D*)fin[0][0]->Get(hname.c_str()))->Clone(name1.c_str());
            h_mc[mc][ch][0][0]->Rebin(rebin);
            fix(h_mc[mc][ch][0][0]);
            for( int s=1; s<nSyst; s++ )
            {
                std::cout<<" "<<systName[s];
                std::string name2 = name1+"_"+systName[s];
                if( s==1 )
                {
                    std::string nameu = name2+tuneName[0];
                    std::string named = name2+tuneName[1];
                    h_mc[mc][ch][s][0] = (TH1D*)h_mc[mc][ch][0][0]->Clone(nameu.c_str());
                    h_mc[mc][ch][s][1] = (TH1D*)h_mc[mc][ch][0][0]->Clone(named.c_str());
                    getHistStatUncNomalized( h_mc[mc][ch][0][0], h_mc[mc][ch][s][0], h_mc[mc][ch][s][1] );
                    //getHistStatUnc( h_mc[mc][ch][0][0], h_mc[mc][ch][s][0], h_mc[mc][ch][s][1] );
                }
                else
                {
                    for( int t=0; t<nTune; t++ )
                    {
                        std::string name3 = name2+tuneName[t];
                        h_mc[mc][ch][s][t] = (TH1D*)((TH1D*)fin[s][t]->Get(hname.c_str()))->Clone(name3.c_str());
                        h_mc[mc][ch][s][t]->Rebin(rebin);
                        fix(h_mc[mc][ch][s][t]);
                    } 
                }
            }
            std::cout<<std::endl;
        }
    }

    // Fitting
    if( doFit )
    {
        std::string nameSig0 = "SigFitted";
        std::string nameBkg0 = "BkgFitted";
        std::cout<<"[INFO] Fitting MC..."<<std::endl;
        TH1D *h_fittedMC[nMC][nCh][nSyst][nTune];

        for( int ch=0; ch<nCh; ch++ )
        {
            std::cout<<"[INFO] Ch"<<chName[ch]<<" ====================================================================================="<<std::endl;
            std::string nameSig1 = nameSig0+chName[ch];
            std::string nameBkg1 = nameBkg0+chName[ch];
            h_fittedMC[sig][ch][0][0] = (TH1D*)h_mc[sig][ch][0][0]->Clone(nameSig1.c_str());
            h_fittedMC[bkg][ch][0][0] = (TH1D*)h_mc[bkg][ch][0][0]->Clone(nameBkg1.c_str());
            fitter_LH( h_fittedMC[sig][ch][0][0], h_fittedMC[bkg][ch][0][0], h_data[ch] );

            for( int s=1; s<nSyst; s++ )
            {
                std::cout<<"[INFO] Ch"<<chName[ch]<<"::"<<systName[s]<<" ========================================================================================== "<<std::endl;
                for( int t=0; t<nTune; t++ )
                {
                    std::string nameSig2 = nameSig1+"_"+systName[s]+tuneName[t];
                    std::string nameBkg2 = nameBkg1+"_"+systName[s]+tuneName[t];
                    h_fittedMC[sig][ch][s][t] = (TH1D*)h_mc[sig][ch][s][t]->Clone(nameSig2.c_str());
                    h_fittedMC[bkg][ch][s][t] = (TH1D*)h_mc[bkg][ch][s][t]->Clone(nameBkg2.c_str());
                    fitter_LH( h_fittedMC[sig][ch][s][t], h_fittedMC[bkg][ch][s][t], h_data[ch] );
                } 
            }
            std::cout<<std::endl;
        }
    }
    fout->Write(); 
}
