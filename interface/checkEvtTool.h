#ifndef CHECKEVTTOOL_H
#define CHECKEVTTOOL_H

#include "iostream"
#include "fstream"
#include "string"
#include <map>

class checkEvtTool
{
    public:
        checkEvtTool( bool debug=false );
        checkEvtTool( std::string myJson, bool debug=false ); 
        ~checkEvtTool(){}

        void addJson( std::string myJson );
        void makeJsonMap();
        void listMyJsons();
        void listGoodLumis();

        bool isGoodEvt( int runNumber,int LumiSec );

    private:
        bool hasJson;
        bool deBug;
        std::vector<std::string> myJsons;
        std::multimap< int, std::pair<int, int> > JsonMap;    // <runNumber, pair<startLumiSec, endLumiSec> >, multimap for same runNumber

        void checkChars( char *nameInput, std::string &nameOutput, int &status );
};

// Constructors 
checkEvtTool::checkEvtTool( bool debug )
{
    deBug = debug; 
    hasJson=false; 
}

checkEvtTool::checkEvtTool( std::string myJson, bool debug )
{
    deBug = debug;
    addJson(myJson); 
}

// Public functions
void checkEvtTool::addJson( std::string myJson )
{
    hasJson=true; 
    myJsons.push_back(myJson);
}

void checkEvtTool::makeJsonMap()
{
    if( !hasJson ){
        printf(">> [ERROR] No any JSON input.\n"); 
        printf("           Please use checkEvt::addJson( std::string &myJson ) to add your JSON file.\n");
        exit(0);
    }
    if( deBug ) listMyJsons(); 

    for( std::vector<std::string>::iterator ijson = myJsons.begin(); ijson != myJsons.end(); ++ijson )
    {
        ifstream JSON(ijson->c_str());
        if(!JSON) {
            printf("[ERROR] Can not found JSON file, %s. Please check if the JSON file exists.\n", ijson->c_str());
            exit(0);
        }   

        char name[128];
        int runNumber = 0;
        int startLumiSec = 0;
        int endLumiSec = 0;
        while(!JSON.eof())
        {
            JSON >> name;
            std::string nameOutput;
            int status(0);
            checkChars( name, nameOutput, status );

            if(status==1){
                runNumber = atoi(nameOutput.c_str());
            }else if(status==2){
                startLumiSec = atoi(nameOutput.c_str());
            }else if(status==3){
                endLumiSec = atoi(nameOutput.c_str());
                JsonMap.insert( std::make_pair( runNumber, std::make_pair(startLumiSec,endLumiSec) ));
            }
        }
        JSON.close();
    }
    if( deBug ) listGoodLumis(); 
}

void checkEvtTool::listMyJsons()
{
    printf(">> [DEBUG] Input Json files:\n");
    for( std::vector<std::string>::iterator ijson = myJsons.begin(); ijson != myJsons.end(); ++ijson )
    {
        printf("           %s\n", ijson->c_str()); 
    }
}

void checkEvtTool::listGoodLumis()
{
    printf(">> [DEBUG] Total Good Run:[LumiStart, lumiEnd]:\n");
    for( std::map< int, std::pair<int, int> >::iterator imap = JsonMap.begin(); imap != JsonMap.end(); ++imap )
    {
        printf("%17d:[%4d, %4d]\n", imap->first, imap->second.first, imap->second.second ); 
    }
}

bool checkEvtTool::isGoodEvt( int runNumber,int LumiSec )
{
    bool isGoodEvt_ = false;
    std::multimap< int, std::pair<int, int> >::iterator JsonMapItr_;
    JsonMapItr_ = JsonMap.find(runNumber);

    if( JsonMapItr_ == JsonMap.end() ){
        isGoodEvt_ = false;
    }else{
        while ( JsonMapItr_!= JsonMap.end()) // loop for same runNumber with different of lumi sections
        {
            if( JsonMapItr_->first != runNumber ) break;
            if( LumiSec >= JsonMapItr_->second.first && 
                    LumiSec <= JsonMapItr_->second.second )
            {
                isGoodEvt_ = true;
            }
            if( isGoodEvt_ ) break;
            JsonMapItr_++;
        }
    }
    return isGoodEvt_;
}

// Private functions
void checkEvtTool::checkChars( char *nameInput, std::string &nameOutput, int &status )
{
    for( unsigned int size_=0; size_<strlen(nameInput); size_++ )
    {
        if( nameInput[0]=='{' || nameInput[0]=='"' ){
            status = 1; // for run number
        }else if( nameInput[0]=='[' ){
            status = 2; // for start evt number
        }else{
            status = 3; // for end evt number
        }
        if( !( nameInput[size_]=='{' || nameInput[size_]=='"' || nameInput[size_]=='[' ||
               nameInput[size_]=='}' || nameInput[size_]==',' || nameInput[size_]==']' ||
               nameInput[size_]==':' )){
            nameOutput += nameInput[size_];
        }
    }
}

#endif
