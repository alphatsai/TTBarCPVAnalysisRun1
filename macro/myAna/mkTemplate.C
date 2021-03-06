#include <string>
#include "../help.C"
#include "templateLFit.C"
#define NPAR 2

std::string fileName="Final_histograms_SemiLeptanic.root";
std::string systDir="Syst";
int CoCH=0;
int ElCH=1;
int MuCH=2;

void mkTemplate( std::string fName, std::string hName, std::string output=".", int rebin=1, bool doFit=true )
{
    const int nCh=3, nSyst=2+7, nMC=2, nTune=2;
    std::string mcName[nMC]={"TTJets__", "BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"}; // sig=0, bkg=1
    std::string chName[nCh]={"", "_El", "_Mu"};
    std::string tuneName[nTune]={"up", "down"};
    std::string systName[nSyst]={"", "Stat", "PU", "JER", "BTagSF", "TopPT", "elID", "muID", "muISO"};
    
    TFile* fin[nSyst][nTune] = new TFile();
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
        if( mc==0 ) name0="SigMC"; 
        else name0="BkgMC"; 
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
                    getHistStatUnc( h_mc[mc][ch][0][0], h_mc[mc][ch][s][0], h_mc[mc][ch][s][1] );
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
        const int sig=0, bkg=1;
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

void mkTemplate( TFile* f, std::string hName, std::string output=".", int rebin=1 )
{
    TFile* fout = new TFile((output+"/Template_"+hName+".root").c_str(), "RECREATE");
    TH1D *h_data, *h_sigMC, *h_bkgMC;
    TH1D *h_data_El, *h_sigMC_El, *h_bkgMC_El;
    TH1D *h_data_Mu, *h_sigMC_Mu, *h_bkgMC_Mu;

    h_data  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA");
    h_sigMC = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_El").c_str()))->Clone("SigMC");
    h_bkgMC = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_El").c_str()))->Clone("BkgMC");

    h_data_El  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hName+"_El").c_str()))->Clone("DATA_El");
    h_sigMC_El = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_El").c_str()))->Clone("SigMC_El");
    h_bkgMC_El = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_El").c_str()))->Clone("BkgMC_El");

    h_data_Mu  = (TH1D*)((TH1D*)f->Get(("DATA_Muon__"+hName+"_Mu").c_str()))->Clone("DATA_Mu");
    h_sigMC_Mu = (TH1D*)((TH1D*)f->Get(("TTJets__"+hName+"_Mu").c_str()))->Clone("SigMC_Mu");
    h_bkgMC_Mu = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hName+"_Mu").c_str()))->Clone("BkgMC_Mu");

    h_data->Rebin(rebin);  h_data_El->Rebin(rebin);  h_data_Mu->Rebin(rebin);
    h_sigMC->Rebin(rebin); h_sigMC_El->Rebin(rebin); h_sigMC_Mu->Rebin(rebin);
    h_bkgMC->Rebin(rebin); h_bkgMC_El->Rebin(rebin); h_bkgMC_Mu->Rebin(rebin);
    
    fix(h_data);  fix(h_data_El);  fix(h_data_Mu);
    fix(h_sigMC); fix(h_sigMC_El); fix(h_sigMC_Mu);
    fix(h_bkgMC); fix(h_bkgMC_El); fix(h_bkgMC_Mu);

    h_data->Add(h_data_Mu);
    h_sigMC->Add(h_sigMC_Mu);
    h_bkgMC->Add(h_bkgMC_Mu);

    getHistStatUnc( h_data,     "DATA_Stat");
    getHistStatUnc( h_data_El,  "DATA_El_Stat");
    getHistStatUnc( h_data_Mu,  "DATA_Mu_Stat");
    getHistStatUnc( h_sigMC,    "SigMC_Stat");
    getHistStatUnc( h_sigMC_El, "SigMC_El_Stat");
    getHistStatUnc( h_sigMC_Mu, "SigMC_Mu_Stat");
    getHistStatUnc( h_bkgMC,    "BkgMC_Stat");
    getHistStatUnc( h_bkgMC_El, "BkgMC_El_Stat");
    getHistStatUnc( h_bkgMC_Mu, "BkgMC_Mu_Stat");

    fout->Write(); 
}
void getHistStatUnc( TH1D* h, TH1D* h_u, TH1D* h_d )
{
    int bins=h->GetNbinsX();
    for( int i=1; i<=bins; i++)
    {
        h_u->SetBinContent(i, h->GetBinContent(i)+h->GetBinError(i));
        h_d->SetBinContent(i, h->GetBinContent(i)-h->GetBinError(i));
    }
}
void getHistStatUnc( TH1D* h, std::string name )
{
    TH1D* h_u = (TH1D*)h->Clone((name+"up").c_str());
    TH1D* h_d = (TH1D*)h->Clone((name+"down").c_str());
    int bins=h->GetNbinsX();
    for( int i=1; i<=bins; i++)
    {
        h_u->SetBinContent(i, h->GetBinContent(i)+h->GetBinError(i));
        h_d->SetBinContent(i, h->GetBinContent(i)-h->GetBinError(i));
    }
}
