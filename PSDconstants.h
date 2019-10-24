#ifndef PSDconstants_H_
#define PSDconstants_H_

#include <string>

//num of rods
const int numRodX = 8;
const int numRodY = 8;

//conversion from ADC to MeVee (MeVee = convFac * ADC + consCon)
//32 rods, noTef, noAmp , 40cm = (2.30628e-05,6.66246e-02) (fac,con)
const double convFac = 1.0;
const double convCon = 0.0;

//Time to group events in mutiple channel as one channel
const double timeOneEv = 0.000001; //1 microsec

//File names of equalizer data
//To get these files, do not use any of the equalizer files that are available in the folder
const std::string flatEqualFile = ""; //if none, then all equalizer is set to 1.0
const std::string csEqualFile   = "";

#endif
