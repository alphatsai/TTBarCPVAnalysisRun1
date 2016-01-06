void mkTemplate( TFile* f, std::string hNameIn, std::string hNameOut="", std::string output="." )
{
    TFile* fout = new TFile((output+"/Template_"+hNameOut+".root").c_str(), "RECREATE");
    TH1D *h_data, *h_sigMC, *h_bkgMC;
    //TH1D *h_data_El, *h_sigMC_El, *h_bkgMC_El;
    //TH1D *h_data_Mu, *h_sigMC_Mu, *h_bkgMC_Mu;

    h_data  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hNameIn+"_El").c_str()))->Clone("DATA");
    h_sigMC = (TH1D*)((TH1D*)f->Get(("TTJets__"+hNameIn+"_El").c_str()))->Clone("SigMC");
    h_bkgMC = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hNameIn+"_El").c_str()))->Clone("BkgMC");

    //h_data_El  = (TH1D*)((TH1D*)f->Get(("DATA_Electron__"+hNameIn+"_El").c_str()))->Clone("DATA_El");
    //h_sigMC_El = (TH1D*)((TH1D*)f->Get(("TTJets__"+hNameIn+"_El").c_str()))->Clone("SigMC_El");
    //h_bkgMC_El = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hNameIn+"_El").c_str()))->Clone("BkgMC_El");

    //h_data_Mu  = (TH1D*)((TH1D*)f->Get(("DATA_Muon__"+hNameIn+"_Mu").c_str()))->Clone("DATA_Mu");
    //h_sigMC_Mu = (TH1D*)((TH1D*)f->Get(("TTJets__"+hNameIn+"_Mu").c_str()))->Clone("SigMC_Mu");
    //h_bkgMC_Mu = (TH1D*)((TH1D*)f->Get(("BkgMC_TTJetsNonSemiLeptMGDecaysExcluded__"+hNameIn+"_Mu").c_str()))->Clone("BkgMC_Mu");;

    //h_data->Add(h_data_Mu);
    //h_sigMC->Add(h_sigMC_Mu);
    //h_bkgMC->Add(h_bkgMC_Mu);

    fout->Write(); 
}
