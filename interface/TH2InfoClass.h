#ifndef TH2INFOCLASS_H
#define TH2INFOCLASS_H

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/TH2Info.h"

using namespace std;

template<typename TH2> 
class TH2InfoClass
{
    public:
        TH2InfoClass( bool deBug=false );
        void addNewTH2( std::string name, std::string title, std::string xtitle, std::string ytitle, std::string xunit, std::string yunit, int binX, double minX, double maxX, int binY, double minY, double maxY);
        void defaultTH2Info(); 
        void ClearTH2Info();
        void CreateTH2();
        void CreateTH2( TFile* f, std::string dir_name=""); 
        void CreateTH2( edm::Service<TFileService> f); 
        void SetTitles(); 
        void Sumw2();
        TH2* GetTH2(std::string name);
        TH2Info GetInfo(std::string name);
        TH2Info GetInfo(int index);

        bool debug;
        int size;

    private:
        map<std::string, TH2*> mapTH2;
        map<std::string, int> indexTH2;
        vector<TH2Info> Info; 
};

    template<typename TH2>
TH2InfoClass<TH2>::TH2InfoClass( bool deBug )
{
    debug=deBug;    
    defaultTH2Info();
    size=Info.size();
}
    template<typename TH2>
void TH2InfoClass<TH2>::ClearTH2Info()
{
    int presize=Info.size();
    Info.clear();
    if( debug ) printf(">> [DEBUG] Clear vector<TH2> Info, size:%d -> %d\n", presize, Info.size());
}
    template<typename TH2>
void TH2InfoClass<TH2>::addNewTH2(std::string name, std::string title, std::string xtitle, std::string ytitle, std::string xunit, std::string yunit, int binX, double minX, double maxX, int binY, double minY, double maxY)
{
    Info.push_back( TH2Info( name, title, xtitle, ytitle, xunit, yunit, binX, minX, maxX, binY, minY, maxY ));
    size=Info.size();
    if( debug ) printf(">> [DEBUG] Add new TH2 %-24s X:bin/min/max:%d/%.2f/%.2f, Y:bin/min/max:%d/%.2f/%.2f, new size is %d\n", name.c_str(), binX, minX, maxX, binY, minY, maxY, size);
}
    template<typename TH2>
void TH2InfoClass<TH2>::defaultTH2Info()
{
  //Info.push_back( TH2Info( Name,           Title,                    xTitle,               yTitle,  xUnit, yUnit,Bin,  Min, Max ));
}

//* Create Histogram
    template<typename TH2> 
void TH2InfoClass<TH2>::CreateTH2()
{
    if( debug ) printf(">> [DEBUG] %d TH2s will be created...\n", size);
    for(int i=0; i<size; i++){
        if( debug ) printf(">> [DEBUG] #%3d %-24s created. X:Bin/Min/Max:%d/%.2f/%.2f, Y:Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].BinX, Info[i].MinX, Info[i].MaxX, Info[i].BinY, Info[i].MinY, Info[i].MaxY );
        indexTH2[Info[i].Name] = i;
        mapTH2[Info[i].Name] = new TH2(Info[i].Name.c_str(),"", Info[i].BinX, Info[i].MinX, Info[i].MaxX, Info[i].BinY, Info[i].MinY, Info[i].MaxY );
    }
}
    template<typename TH2> 
void TH2InfoClass<TH2>::CreateTH2( edm::Service<TFileService> f )
{
    if( debug ) printf(">> [DEBUG] %d TH2s will be created...\n", size);
    for(int i=0; i<size; i++){
        if( debug ) printf(">> [DEBUG] #%3d %-24s created. X:Bin/Min/Max:%d/%.2f/%.2f, Y:Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].BinX, Info[i].MinX, Info[i].MaxX, Info[i].BinY, Info[i].MinY, Info[i].MaxY );
        indexTH2[Info[i].Name] = i;
        mapTH2[Info[i].Name] = f->make<TH2>(Info[i].Name.c_str(),"", Info[i].BinX, Info[i].MinX, Info[i].MaxX, Info[i].BinY, Info[i].MinY, Info[i].MaxY );
    }
}
    template<typename TH2> 
void TH2InfoClass<TH2>::CreateTH2( TFile* f, std::string dir_name)
{
    printf(">> %d TH2s will be got...\n", size);
    for(int i=0; i<size; i++){ 
        if( debug ) printf(">> [DEBUG] #%3d %-24s created. X:Bin/Min/Max:%d/%.2f/%.2f, Y:Bin/Min/Max:%d/%.2f/%.2f\n",i, Info[i].Name.c_str(), Info[i].BinX, Info[i].MinX, Info[i].MaxX, Info[i].BinY, Info[i].MinY, Info[i].MaxY);
        indexTH2[Info[i].Name] = i;
        mapTH2[Info[i].Name] =(TH2*)f->Get( (dir_name+Info[i].Name).c_str() );
    }
}

//* Set some option for Histogram
template<typename TH2> 
void TH2InfoClass<TH2>::SetTitles(){
    for(int i=0; i<size; i++){ 
        //mapTH2[Info[i].Name]->SetTile(Info[i]._title.c_str());
        string xt, yt;
        if( Info[i].xUnit.size() == 0 ) xt = Info[i].xTitle;
        else xt = Info[i].xTitle+" ["+Info[i].xUnit+"]";    
        if( Info[i].yUnit.size() == 0 ) yt = Info[i].yTitle;
        else yt = Info[i].yTitle+" ["+Info[i].yUnit+"]";    
        mapTH2[Info[i].Name]->SetXtitle( xt.c_str() );
        mapTH2[Info[i].Name]->SetYtitle( yt.c_str() );
    }
}

template<typename TH2> 
void TH2InfoClass<TH2>::Sumw2(){
    for(int i=0; i<size; i++){ 
        mapTH2.find(Info[i].Name)->second->Sumw2();
    }
}

//* Get Histogram
template<typename TH2> 
TH2* TH2InfoClass<TH2>::GetTH2(std::string name){
    if( mapTH2.find(name) == mapTH2.end() ){
        printf(">> [ERROR] %s is not found in TH2InfoClass::GetTH2(std::string)\n", name.c_str());
        return NULL;
    }else{
        return mapTH2.find(name)->second;
    }
}

//* Get Info variables
template<typename TH2> 
TH2Info TH2InfoClass<TH2>::GetInfo(std::string name){
    TH2Info info;
    if( mapTH2.find(name) == mapTH2.end() ){
        printf(">> [ERROR] %s is not found in TH2InfoClass::GetInfo(std::string)\n", name.c_str());
    }else{    
        info=Info[indexTH2.find(name)->second];
    }
    return info;
}
template<typename TH2> 
TH2Info TH2InfoClass<TH2>::GetInfo(int index){
    TH2Info info;
    if( index >= size ){
        printf(">> [ERROR] %d is over the size(%d) in TH2InfoClass::GetInfo(int)\n", index, size);
    }else{
        info=Info[index];
    }
    return info;
}

#endif
//
