#ifndef PSDanalyze_h
#define PSDanalyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "PSDhisto.h"
#include "PSDequalizer.h"

using namespace std;

class PSDanalyze
{
    private :
    
    double tot,tail,timestamp,erg;
    int xpos,ypos,SipmNum,asicNum,digiNum;
    int chNum;
    
    public :
    
    double tNow;
    PSDanalyze(string filename);
    virtual ~PSDanalyze();
    virtual void Loop(string fileName, PSDhisto histos, PSDequalizer myEqual);
    
};

#endif
