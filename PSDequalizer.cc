#include "PSDconstants.h"
#include "PSDequalizer.h"

#include <fstream>
#include <TH2.h>
#include <iostream>

using namespace std;

void PSDequalizer::ReadDataEqualizer()
{
    
    hSimpEq0 = new TH2D("hSimpEq0","",numRodX,0,numRodX,numRodY,0,numRodY);
    hSimpEq1 = new TH2D("hSimpEq1","",numRodX,0,numRodX,numRodY,0,numRodY);
    
    //////////////////////////////////
    //Flat field equalizer
    //////////////////////////////////
    
    //open data file
    ifstream rodEqCalFileIn(flatEqualFile.c_str());
    
    //if the flat field equalizer data exists
    if ( rodEqCalFileIn.good() ){
        
        double nextval=0.0;
        
        //need to ask Steven, why the x and y are flipped..?
        for (int y=1;y<numRodX+1;y++) {
            for (int x=1;x<numRodY+1;x++) {
                rodEqCalFileIn >> nextval;
                hSimpEq0->SetBinContent(x,y,nextval);
            }
        }
        for (int y=1;y<numRodX+1;y++) {
            for (int x=1;x<numRodY+1;x++) {
                rodEqCalFileIn >> nextval;
                hSimpEq1->SetBinContent(x,y,nextval);
            }
        }
        
    }
    //if the flat field equalizer data does not exists
    else {
        for (int x=1;x<numRodX+1;x++) {
            for (int y=1;y<numRodY+1;y++) {
                hSimpEq0->SetBinContent(x,y,1.0);
                hSimpEq1->SetBinContent(x,y,1.0);
            }
        }
    }
    
    //close file
    rodEqCalFileIn.close();
    
    //////////////////////////////////
    //Cs137 center run equalizer
    //////////////////////////////////
    
    //open data file
    ifstream calfilein(csEqualFile.c_str());
    
    //if the cs137 center run equalizer exists
    //replace the equalizer with the Cs center data
    if (calfilein.good()) {
        
        double FFP=0.0;
        
        for (int x=1;x<numRodX+1;x++) {
            for (int y=1;y<numRodY+1;y++) {
                calfilein >> FFP;
                hSimpEq0->SetBinContent(x,y,FFP);
            }
        }
        
        for (int x=1;x<numRodX+1;x++) {
            for (int y=1;y<numRodY+1;y++) {
                calfilein >> FFP;
                hSimpEq1->SetBinContent(x,y,FFP);
            }
        }
    }
    
    //close file
    calfilein.close();
    
}

double PSDequalizer::GetEqualizedErg(double erg, int x, int y, int sipmNum)
{
    double ergEq = 0.0;
    
    if( sipmNum == 0 ) ergEq = erg * hSimpEq0->GetBinContent(x+1,y+1);
    else               ergEq = erg * hSimpEq1->GetBinContent(x+1,y+1);
    
    return ergEq;
}

