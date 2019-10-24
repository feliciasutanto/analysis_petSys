#ifndef _PSDHISTO_H_
#define _PSDHISTO_H_

#include <TCutG.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include "PSDconstants.h"

using namespace std;

class PSDhisto
{
private:
    
    double calTot;
    double rawTot,rawTail;
    double psp,pspTot;
    double sumTot,sumTail,sumCalTot;
    int mulNeu,mulGam;
    double sumTot0,sumTot1;
    double minTime0,minTime1;
    
public:
    
    //Fiducial cuts are applied to the equalized spectrum in adc
    TCutG *gammaCut;
    TCutG *neutronCut;
    TCutG *tot_gammaCut;
    TCutG *tot_neutronCut;
    
    //Raw spectrum in ADC for each channel
    TH1F *energyHisto0[numRodX][numRodY]; //SiPM 0
    TH1F *energyHisto1[numRodX][numRodY]; //SiPM 1
    
    //Raw spectrum in ADC for each rod
    TH2F *scatterHisto[numRodX][numRodY];
    TH1F *energyHisto [numRodX][numRodY];
    TH1F *neutronHisto[numRodX][numRodY];
    TH1F *gammaHisto  [numRodX][numRodY];
    TH1F *qRatioHisto [numRodX][numRodY];
    
    //Calibrated spectrum in MeVee for each rod
    TH2F *cal_scatterHisto[numRodX][numRodY];
    TH1F *cal_energyHisto [numRodX][numRodY];
    TH1F *cal_neutronHisto[numRodX][numRodY];
    TH1F *cal_gammaHisto  [numRodX][numRodY];
    
    //Raw spectrum in ADC for all rods
    TH2F *tot_scatterHisto;
    TH1F *tot_energyHisto;
    TH1F *tot_neutronHisto;
    TH1F *tot_gammaHisto;
    TH1F *tot_qRatioHisto;
    TH2F *tot_qRatioErgHisto;
    TH2F *tot_erg0erg1Histo;
    
    //Calibrated spectrum in MeVee for all rods
    TH2F *tot_cal_scatterHisto;
    TH1F *tot_cal_energyHisto;
    TH1F *tot_cal_neutronHisto;
    TH1F *tot_cal_gammaHisto;
    
    //Multiplicity vs Energy spectrum
    TH2F *mul_erg_gammaHisto;
    TH2F *mul_erg_neutronHisto;
    
    //Time between two SiPM
    TH1F * timeBetSipm;
    
    //Active digitizer
    TH1F *digiActive;
    TH2F *posActive0;
    TH2F *posActive1;
    
    //Functions
    void createHistos();
    void fillHistos(vector< vector<double> > &adcTail0,
                    vector< vector<double> > &adcTot0,
                    vector< vector<double> > &adcTail1,
                    vector< vector<double> > &adcTot1,
                    vector< vector<double> > &time0,
                    vector< vector<double> > &time1);
    void createAndWriteRootFile( string finalOutPath );
    void MakeCuts();
    void fillDigi(double digiNum){ digiActive->Fill(digiNum); };
    void fillPos0(double x, double y){ posActive0->Fill(x,y); };
    void fillPos1(double x, double y){ posActive1->Fill(x,y); };
    
    
};

#endif





