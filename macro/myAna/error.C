void error()
{
    // Electron
    float sig0 = 31354.3;
    float bkg0 = 3322.4;
    float sig  = 28934;
    float bkg  = 2346;
    float cEff = 0.73;
    float sigE0 = 270.1;
    float bkgE0 = 211.9;

    ////Muon 
    //float sig0 = 36346.1;
    //float bkg0 = 3636.9;
    //float sig  = 33745;
    //float bkg  = 2764;
    //float cEff = 0.79;
    //float sigE0 = 324.8;
    //float bkgE0 = 269.7;

    float sum0 = sig0+bkg0;
    float sum  = sig+bkg;
    float sigE = sigE0*sig/sig0;
    float bkgE = bkgE0*bkg/bkg0;
    float err0 = sqrt(sigE0*sigE0+bkgE0*bkgE0-2*cEff*sigE0*bkgE0);
    float err  = sqrt(sigE*sigE+bkgE*bkgE-2*cEff*sigE*bkgE);
    float pure0 = sig0/sum0;
    float pure  = sig/sum;
    float pErr0 = sqrt(bkg0*bkg0*bkgE0*bkgE0+sig0*sig0*bkgE0*bkgE0+2*sig0*bkg0*cEff*sigE0*bkgE0)/sum0/sum0;
    float pErr = sqrt(bkg*bkg*bkgE*bkgE+sig*sig*bkgE*bkgE+2*sig*bkg*cEff*sigE*bkgE)/sum/sum;
    
    printf("1. Sig: %.0f (%.0f), bkg: %.0f (%.0f), sum: %.0f (%.0f %.0f), purity: %.1f (%.3f)\n", sig0, sigE0, bkg0, bkgE0, sum0, err0, sqrt(sum0), pure0*100, pErr0*100);     
    printf("2. Sig: %.0f (%.0f), bkg: %.0f (%.0f), sum: %.0f (%.0f %.0f), purity: %.1f (%.3f)\n", sig, sigE, bkg, bkgE, sum, err, sqrt(sum), pure*100, pErr*100);     

}
