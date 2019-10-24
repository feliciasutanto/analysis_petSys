#ifndef _PSDEQUALIZER_H_
#define _PSDEQUALIZER_H_

#include <TH2.h>

using namespace std;

class PSDequalizer
{
private:
    
public:
    
    TH2D *hSimpEq0;
    TH2D *hSimpEq1;
    
    void    ReadDataEqualizer();
    double  GetEqualizedErg(double erg, int x, int y, int sipmNum);

};

#endif
