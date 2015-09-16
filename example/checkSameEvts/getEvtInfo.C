void getEvtInfo( TFile* f, std::string dataname; std::string treeName="SemiLeptanic/tree" )
{
        int RunNo;
        int LumiNo;
        int isEleEvt;
        int isMuonEvt;
        long long int EvtNo;

        TTree* t =(TTree*)f->Get(treeName.c_str());
        t->SetBranchAddress("EvtInfo.RunNo",     &RunNo     );
        t->SetBranchAddress("EvtInfo.LumiNo",    &LumiNo    );
        t->SetBranchAddress("EvtInfo.EvtNo",     &EvtNo     );
        t->SetBranchAddress("EvtInfo.isEleEvt",  &isEleEvt  );
        t->SetBranchAddress("EvtInfo.isMuonEvt", &isMuonEvt );

        cout<<">> Writing to "<<dataname<<"_MuChannel.txt..."<<endl;
        cout<<">> Writing to "<<dataname<<"_ElChannel.txt..."<<endl;

        int numElEvts=0;
        int numMuEvts=0;
        FILE * outEl = fopen((dataname+"_ElChannel.txt").c_str(),"w");
        FILE * outMu = fopen((dataname+"_MuChannel.txt").c_str(),"w");
        fprintf( outEl, "# %6s %7s %12s\n", "RunNo", "LumiNo", "EvtNo");
        fprintf( outMu, "# %6s %7s %12s\n", "RunNo", "LumiNo", "EvtNo");
        for( int evt=0; evt<t->GetEntries(); evt++)
        {
            t->GetEntry(evt);
            if( isEleEvt && !isMuonEvt )
            {
                fprintf( outEl, "%8d %7d %12ld\n", RunNo, LumiNo, EvtNo);
                numElEvts++;
            }
            else if( !isEleEvt && isMuonEvt )
            {
                fprintf( outMu, "%8d %7d %12ld\n", RunNo, LumiNo, EvtNo);
                numMuEvts++;
            }
            else{ cout<<">> [ERROR] There is ambiguous channel[El, Mu] ["<<isEleEvt<<", "<<isMuonEvt<<"]"<<endl; }
        }
        fprintf( outEl, "# Total events %d\n", numElEvts);
        fprintf( outMu, "# Total events %d\n", numMuEvts);
        fclose(outEl);
        fclose(outMu);
        delete t;

        cout<<">> Done!"<<endl;
}
