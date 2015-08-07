#ifndef TH2INFO_H
#define TH2INFO_H

#include <map>
#include <string>
#include <TFile.h>

using namespace std;

class TH2Info{
    public:
        TH2Info();
        TH2Info(
                std::string name, 
                std::string title, 
                std::string xtitle, 
                std::string ytitle, 
                std::string xunit, 
                std::string yunit, 
                int     binX, 
                double  minX, 
                double  maxX,  
                int     binY, 
                double  minY, 
                double  maxY  
               ) : 
            Name  (name), 
            Title (title), 
            xTitle(xtitle), 
            yTitle(ytitle), 
            xUnit (xunit), 
            yUnit (yunit), 
            BinX   (binX), 
            MinX   (minX), 
            MaxX   (maxX), 
            BinY   (binY), 
            MinY   (minY), 
            MaxY   (maxY){} 

        std::string Name;
        std::string Title;
        std::string xTitle;
        std::string yTitle;
        std::string xUnit;
        std::string yUnit;
        int     BinX;
        double  MinX;
        double  MaxX;
        int     BinY;
        double  MinY;
        double  MaxY;
};
#endif 

