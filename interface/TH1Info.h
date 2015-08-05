#ifndef TH1INFO_H
#define TH1INFO_H

#include <map>
#include <string>
#include <TFile.h>

using namespace std;

class TH1Info{
    public:
        TH1Info();
        TH1Info(
                std::string name, 
                std::string title, 
                std::string xtitle, 
                std::string ytitle, 
                std::string xunit, 
                std::string yunit, 
                int     bin, 
                double  min, 
                double  max  
               ) : 
            Name  (name), 
            Title (title), 
            xTitle(xtitle), 
            yTitle(ytitle), 
            xUnit (xunit), 
            yUnit (yunit), 
            Bin   (bin), 
            Min   (min), 
            Max   (max){} 

        std::string Name;
        std::string Title;
        std::string xTitle;
        std::string yTitle;
        std::string xUnit;
        std::string yUnit;
        int     Bin;
        double  Min;
        double  Max;
};
#endif 

