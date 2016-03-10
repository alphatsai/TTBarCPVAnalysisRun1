void mkCombinedHist( TFile* fin )
{
    TH1D* h_sigEl, *h_sigMu, *h_dataEl, *h_dataMu, *h_bkgEl, *h_bkgMu;
    h_sigEl  = (TH1D*)fin->Get("SigMC_El");
    h_sigMu  = (TH1D*)fin->Get("SigMC_Mu");
    h_bkgEl  = (TH1D*)fin->Get("BkgMC_El");
    h_bkgMu  = (TH1D*)fin->Get("BkgMC_Mu");
    h_dataEl = (TH1D*)fin->Get("DATA_El");
    h_dataMu = (TH1D*)fin->Get("DATA_Mu");

    TFile* fout = new TFile("HistCombined.root", "RECREATE");

    // NOTE: all the hist shall have same binning
    int bins   = h_sigEl->GetNbinsX();
    float xmin = h_sigEl->GetXaxis()->GetBinLowEdge(h_sigEl->GetXaxis()->GetFirst());
    float xmax = h_sigEl->GetXaxis()->GetBinUpEdge(h_sigEl->GetXaxis()->GetLast());
    TH1D* h_sig  = new TH1D("SigMC", "", 2*bins, xmin, 2*xmax);
    TH1D* h_bkg  = new TH1D("BkgMC", "", 2*bins, xmin, 2*xmax);
    TH1D* h_data = new TH1D("DATA",  "", 2*bins, xmin, 2*xmax);

    //h_sig ->Add(h_sigEl);
    //h_bkg ->Add(h_bkgEl);
    //h_data->Add(h_dataEl);

    for( int b=1; b<=bins; b++ )
    {
        int bb = b+bins;
        h_sig ->SetBinContent(  b, h_sigEl ->GetBinContent(b));
        h_sig ->SetBinContent( bb, h_sigMu ->GetBinContent(b));
        h_bkg ->SetBinContent(  b, h_bkgEl ->GetBinContent(b));
        h_bkg ->SetBinContent( bb, h_bkgMu ->GetBinContent(b));
        h_data->SetBinContent(  b, h_dataEl->GetBinContent(b));
        h_data->SetBinContent( bb, h_dataMu->GetBinContent(b));
    }
    fout->Write();
}
void mkCuttedHist( TFile* fin, float cmin=0, float cmax=0 )
{   
    if( cmin == cmax ) return;
    TFile* fout = new TFile("HistCutted.root", "RECREATE");
    TH1D* h_sigEl, *h_sigMu, *h_dataEl, *h_dataMu, *h_bkgEl, *h_bkgMu;
    h_sigEl  = (TH1D*)((TH1D*)fin->Get("SigMC_El"))->Clone();
    h_sigMu  = (TH1D*)((TH1D*)fin->Get("SigMC_Mu"))->Clone();
    h_bkgEl  = (TH1D*)((TH1D*)fin->Get("BkgMC_El"))->Clone();
    h_bkgMu  = (TH1D*)((TH1D*)fin->Get("BkgMC_Mu"))->Clone();
    h_dataEl = (TH1D*)((TH1D*)fin->Get("DATA_El"))->Clone();
    h_dataMu = (TH1D*)((TH1D*)fin->Get("DATA_Mu"))->Clone();  

    int bins   = h_sigEl->GetNbinsX();
    float xmin = h_sigEl->GetXaxis()->GetBinLowEdge(h_sigEl->GetXaxis()->GetFirst());
    float xmax = h_sigEl->GetXaxis()->GetBinUpEdge(h_sigEl->GetXaxis()->GetLast());
    for( int b=1; b<=bins; b++ )
    {
        if( h_sigEl->GetXaxis()->GetBinLowEdge(b) >= cmax || h_sigEl->GetXaxis()->GetBinLowEdge(b) < cmin )
        {
            h_sigEl ->SetBinContent( b, 0 );  
            h_sigMu ->SetBinContent( b, 0 );  
            h_bkgEl ->SetBinContent( b, 0 );  
            h_bkgMu ->SetBinContent( b, 0 );  
            h_dataEl->SetBinContent( b, 0 );  
            h_dataMu->SetBinContent( b, 0 );  
            h_sigEl ->SetBinError( b, 0 );  
            h_sigMu ->SetBinError( b, 0 );  
            h_bkgEl ->SetBinError( b, 0 );  
            h_bkgMu ->SetBinError( b, 0 );  
            h_dataEl->SetBinError( b, 0 );  
            h_dataMu->SetBinError( b, 0 );  
        }
    }
    fout->Write();
}
