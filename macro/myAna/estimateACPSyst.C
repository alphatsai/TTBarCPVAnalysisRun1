float estimateACP( float nData, float nBkgP, float nBkgN, float AcpData=0 )
{
    float nDataO[2]={ nData*(1.+AcpData/100.)/2., nData*(1.-AcpData/100.)/2. };
    float nSigO[2] = { nDataO[0]-nBkgP, nDataO[1]-nBkgN };
    return (nSigO[0]-nSigO[1])/(nSigO[0]+nSigO[1])*100;
}
void estimateACPSyst( std::string output="./SystUncertaintiesEst.txt", float AcpData=0 )
{
    std::string sysNameEl[8]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "elID" }; 
    std::string sysNameMu[9]={"TopMatch", "TopScale", "TopPT", "PU", "JER", "JES", "BTagSF", "muID", "muISO"}; 
    float sysEl[2][8]={ {2.66, 41.57, 23.78,-1.20,5.36,-6.61,-5.50,-0.15},
                        {8.27,-30.61,-24.57, 0.30,4.12,16.02, 5.73, 0.15}}; 
    float sysMu[2][9]={ {10.50, 45.00, 29.60, 0.55,-2.59,-5.96,-5.37,-2.14,-0.27},
                        { 9.17,-31.81,-30.96,-0.93, 4.45,22.91, 5.64, 2.16, 0.21}}; 
  
    std::string oName[4]={"O2","O3","O4","O7"};
    float nBkgElP[4]={889.319,948.25,889.616,931.472}; 
    float nBkgElN[4]={938.568,879.637,938.271,896.415}; 
    float nBkgMuP[4]={1034.71,1076.58,1041.12,1068.25}; 
    float nBkgMuN[4]={1137.22,1095.35,1130.8,1103.68};
    float wrtEl=2063./1828.; 
    float wrtMu=1829./2171.;

    FILE* outTxt;
    outTxt = fopen(output.c_str(),"w");
    fprintf( outTxt, "Electron channel\n" );
    fprintf( outTxt, "%9s", " " );
    for( int s=0; s<8; s++){
        fprintf( outTxt, "%9s", sysNameEl[s].c_str() );
    }
    fprintf( outTxt, "\n");
    // Electron
    for( int o=0; o<4; o++){ 
        fprintf( outTxt, "%9s\n", oName[o].c_str() );
        float acp = estimateACP( 32469., nBkgElP[o]*wrtEl, nBkgElN[o]*wrtEl, AcpData );
        for( int t=0; t<2; t++){
            if( t==0 )
                fprintf( outTxt, "%9s", "+1sigma" );
            else
                fprintf( outTxt, "%9s", "-1sigma" );
            for( int s=0; s<8; s++){
                float acpSys = estimateACP( 32469., nBkgElP[o]*wrtEl*(1.+sysEl[t][s]/100.), nBkgElN[o]*wrtEl*(1.+sysEl[t][s]/100.), AcpData );
                fprintf( outTxt, "%+9.2f", acpSys-acp );
            }
            fprintf( outTxt, "\n");
        }
    }
    fprintf( outTxt, "\n");
    // Muon 
    fprintf( outTxt, "Muon channel\n" );
    fprintf( outTxt, "%9s", " " );
    for( int s=0; s<9; s++){
        fprintf( outTxt, "%9s", sysNameMu[s].c_str() );
    }
    fprintf( outTxt, "\n");
    for( int o=0; o<4; o++){ 
        fprintf( outTxt, "%9s\n", oName[o].c_str() );
        float acp = estimateACP( 38389., nBkgMuP[o]*wrtMu, nBkgMuN[o]*wrtMu, AcpData );
        for( int t=0; t<2; t++){
            if( t==0 )
                fprintf( outTxt, "%9s", "+1sigma" );
            else
                fprintf( outTxt, "%9s", "-1sigma" );
            for( int s=0; s<9; s++){
                float acpSys = estimateACP( 38389., nBkgMuP[o]*wrtMu*(1.+sysMu[t][s]/100.), nBkgMuN[o]*wrtMu*(1.+sysMu[t][s]/100.), AcpData );
                fprintf( outTxt, "%+9.2f", acpSys-acp );
            }
            fprintf( outTxt, "\n");
        }
    }
}
