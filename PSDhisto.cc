#include "PSDhisto.h"
#include "PSDconstants.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCutG.h>

using namespace std;

void PSDhisto::createHistos()
{
    
    for(int i=0; i<numRodX; i++){
        for(int j=0; j<numRodY; j++){
            
            stringstream SSname;
            TString histoName;
            
            ///////////////////////////////////////
            //Raw spectrum in ADC for each channel
            ///////////////////////////////////////
            
            SSname << "rod_" << i <<"_"<< j <<"_SiPM0_energyHisto";
            histoName = SSname.str();
            energyHisto0[i][j] = new TH1F(histoName,";Light output (ADC);Counts",1e4,0,1e2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_SiPM1_energyHisto";
            histoName = SSname.str();
            energyHisto1[i][j] = new TH1F(histoName,";Light output (ADC);Counts",1e4,0,1e2);
            
            ///////////////////////////////////////
            //Raw spectrum in ADC for each rod
            ///////////////////////////////////////
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_scatterHisto";
            histoName = SSname.str();
            scatterHisto[i][j] = new TH2F(histoName,";Light output (ADC);PSP",1e4,0,1e2,1e3,-2,2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_energyHisto";
            histoName = SSname.str();
            energyHisto[i][j] = new TH1F(histoName,";Light output (ADC);Counts",1e4,0,1e2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_gammaHisto";
            histoName = SSname.str();
            gammaHisto[i][j] = new TH1F(histoName,";Light output (ADC);Counts",1e4,0,1e2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_neutronHisto";
            histoName = SSname.str();
            neutronHisto[i][j] = new TH1F(histoName,";Light output (ADC);Counts",1e4,0,1e2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_qRatioHisto";
            histoName = SSname.str();
            qRatioHisto[i][j] = new TH1F(histoName,";Charge ratio;Counts",1e3,-2,2);
            
            ///////////////////////////////////////
            //Calibrated spectrum in MeVee for each rod
            ///////////////////////////////////////
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_cal_scatterHisto";
            histoName = SSname.str();
            cal_scatterHisto[i][j] = new TH2F(histoName,";Light output (MeVee);PSP",1e4,0,50,1e3,-2,2);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_cal_energyHisto";
            histoName = SSname.str();
            cal_energyHisto[i][j] = new TH1F(histoName,";Light output (MeVee);Counts",1e4,0,50);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_cal_gammaHisto";
            histoName = SSname.str();
            cal_gammaHisto[i][j] = new TH1F(histoName,";Light output (MeVee);Counts",1e4,0,50);
            
            SSname.str("");
            SSname << "rod_" << i <<"_"<< j <<"_cal_neutronHisto";
            histoName = SSname.str();
            cal_neutronHisto[i][j] = new TH1F(histoName,";Light output (MeVee);Counts",1e4,0,50);
        }
    }
    
    ///////////////////////////////////////
    //Raw spectrum in ADC for all rod
    ///////////////////////////////////////
    
    tot_scatterHisto= new TH2F("tot_scatterHisto",";Light output (ADC);PSP"   ,1e4,0,1e2,1e3,-2,2);
    tot_energyHisto = new TH1F("tot_energyHisto" ,";Light output (ADC);Counts",1e4,0,1e2);
    tot_gammaHisto  = new TH1F("tot_gammaHisto"  ,";Light output (ADC);Counts",1e4,0,1e2);
    tot_neutronHisto= new TH1F("tot_neutronHisto",";Light output (ADC);Counts",1e4,0,1e2);
    tot_qRatioHisto = new TH1F("tot_qRatioHisto",";Charge ratio;Counts",1e3,-2,2);
    tot_qRatioErgHisto= new TH2F("tot_qRatioErgHisto",";Light output (ADC);Charge ratio",1e4,0,1e2,1e3,-2,2);
    tot_erg0erg1Histo = new TH2F("tot_erg0erg1Histo",";Light output 0 (ADC);Light output 1 (ADC)",1e4,0,1e2,1e4,0,1e2);
    
    ///////////////////////////////////////
    //Calibrated spectrum in MeVee for all rod
    ///////////////////////////////////////
    
    tot_cal_scatterHisto= new TH2F("tot_cal_scatterHisto",";Light output (MeVee);PSP"   ,1e4,0,50,1e3,-2,2);
    tot_cal_energyHisto = new TH1F("tot_cal_energyHisto" ,";Light output (MeVee);Counts",1e4,0,50);
    tot_cal_gammaHisto  = new TH1F("tot_cal_gammaHisto"  ,";Light output (MeVee);Counts",1e4,0,50);
    tot_cal_neutronHisto= new TH1F("tot_cal_neutronHisto",";Light output (MeVee);Counts",1e4,0,50);
    
    ///////////////////////////////////////
    //Multiplicity per gamma/neutron event
    ///////////////////////////////////////
    
    mul_erg_gammaHisto  = new TH2F("mul_erg_gammaHisto"  ,";Light output (MeVee);Rod multiplicity",
                                   1e4,0,50,numRodX*numRodY+1,-0.5,numRodX*numRodY+0.5);
    mul_erg_neutronHisto= new TH2F("mul_erg_neutronHisto",";Light output (MeVee);Rod multiplicity",
                                   1e4,0,50,numRodX*numRodY+1,-0.5,numRodX*numRodY+0.5);
    
    ///////////////////////////////////////
    //Time between two SiPM
    ///////////////////////////////////////
    
    timeBetSipm = new TH1F("timeBetSipm",";Time (ps); Counts",1e4,-5e3,5e3);
    
    //digitizer channel active
    digiActive = new TH1F("digiActive",";ch no;Counts",1e4,-10,1e3);
    posActive0 = new TH2F("posActive0",";x;y",10,-1.5,8.5,10,-1.5,8.5);
    posActive1 = new TH2F("posActive1",";x;y",10,-1.5,8.5,10,-1.5,8.5);
}

void PSDhisto::MakeCuts()
{
    
    gammaCut  = new TCutG("gammaCut",0);
    neutronCut= new TCutG("neutronCut",0);
    tot_gammaCut  = new TCutG("tot_gammaCut",0);
    tot_neutronCut= new TCutG("tot_neutronCut",0);
    
}

void PSDhisto::fillHistos(vector< vector<double> > &adcTail0,
                          vector< vector<double> > &adcTot0,
                          vector< vector<double> > &adcTail1,
                          vector< vector<double> > &adcTot1,
                          vector< vector<double> > &time0,
                          vector< vector<double> > &time1)
{
    //Initialization
    sumTot = 0.0; sumTail = 0.0;
    mulNeu = 0; mulGam = 0;
    sumTot0 = 0.0; sumTot1 = 0.0;
    vector <int> ineu,jneu;
    vector <double> rawNeu;
    minTime0 = 1.e30; minTime1 = 1.e30;
    double tot0,tot1,tail0,tail1;
    
    //Fill up histogram for each rod
    for(int i=0; i<numRodX; i++){
        for(int j=0; j<numRodY; j++){
            
            //get the minimum time, but not zero time
            //also do we want to limit the energy as well to just include the photopeak of 0.511
            if( time0[i][j]<minTime0 && time0[i][j]!=0.0)// && adcTot0[i][j]>170. && adcTot0[i][j]<175. )
            { minTime0=time0[i][j]; }
            if( time1[i][j]<minTime1 && time1[i][j]!=0.0)// && adcTot1[i][j]>170. && adcTot1[i][j]<175. )
            { minTime1=time1[i][j]; }
            
            //get Raw total and tail
            rawTot = sqrt( adcTot0[i][j]  * adcTot1[i][j]  );
            rawTail= sqrt( adcTail0[i][j] * adcTail1[i][j] );
            //sometimes, the tail and total are swapped...
            /*if(adcTot0[i][j]>adcTail0[i][j]){tot0=adcTot0[i][j];tail0=adcTail0[i][j];}
             else                            {tail0=adcTot0[i][j];tot0=adcTail0[i][j];}
             if(adcTot1[i][j]>adcTail1[i][j]){tot1=adcTot1[i][j];tail1=adcTail1[i][j];}
             else                            {tail1=adcTot1[i][j];tot1=adcTail1[i][j];}
             rawTot = sqrt( tot0  * tot1  );
             rawTail= sqrt( tail0 * tail1 );*/
            
            /*//Fill up histo for each channel (ADC)
             if( adcTot0[i][j] > 0.0 ) energyHisto0[i][j]->Fill( tot0 ); //adcTot0[i][j] );
             if( adcTot1[i][j] > 0.0 ) energyHisto1[i][j]->Fill( tot1 ); //adcTot1[i][j] );*/
            if( adcTot0[i][j] > 0.0 ) energyHisto0[i][j]->Fill( adcTot0[i][j] );
            if( adcTot1[i][j] > 0.0 ) energyHisto1[i][j]->Fill( adcTot1[i][j] );
            
            //Fill up histo for each rod
            if( rawTot > 0.0 ){
                
                psp    = (rawTot-rawTail)/rawTot;
                calTot = (convFac*rawTot)+convCon;
                
                //cout << "psp,rawStart,rawTot: " << psp << "," << rawTail << "," << rawTot << endl;
                //cout << "x,y: " << i<< "," << j << endl;
                
                //Fill up histo for each rod (ADC)
                energyHisto[i][j]  ->Fill(rawTot     );
                scatterHisto[i][j] ->Fill(rawTot, psp);
                if(gammaCut  ->IsInside(rawTot,psp)){ gammaHisto[i][j]  ->Fill(rawTot); mulGam++; }
                if(neutronCut->IsInside(rawTot,psp)){ neutronHisto[i][j]->Fill(rawTot); mulNeu++; }
                qRatioHisto[i][j]->Fill( adcTot0[i][j] / (adcTot0[i][j]+adcTot1[i][j]) );
                
                //Fill up histo for each rod (MeVee)
                cal_energyHisto[i][j] ->Fill(calTot);
                cal_scatterHisto[i][j]->Fill(calTot, psp);
                if(gammaCut  ->IsInside(rawTot,psp)){ cal_gammaHisto[i][j]  ->Fill(calTot); }
                if(neutronCut->IsInside(rawTot,psp)){ cal_neutronHisto[i][j]->Fill(calTot); }
                
                sumTot  += rawTot;
                sumTail += rawTail;
                
                if(rawTot>sumTot){sumTot=rawTot;}
                if(rawTail>sumTail){sumTail=rawTail;}
                
                //sumTot0 += adcTot0[i][j];
                //sumTot1 += adcTot1[i][j];
                
            }
        }
    }
    
    pspTot = (sumTot-sumTail)/sumTot;
    sumCalTot = (sumTot*convFac)+convCon;
    
    if( sumTot0>0.0 || sumTot1>0.0 ){
        tot_qRatioHisto->Fill(sumTot0/(sumTot0+sumTot1));
        tot_qRatioErgHisto->Fill(sumTot,sumTot0/(sumTot0+sumTot1));
        tot_erg0erg1Histo ->Fill(sumTot0,sumTot1);
    }
    
    if( sumTot > 0.0 ){
        
        //Fill up histogram for sum of rods (ADC)
        tot_energyHisto ->Fill(sumTot);
        tot_scatterHisto->Fill(sumTot, pspTot);
        if(tot_gammaCut  ->IsInside(sumTot,pspTot)) tot_gammaHisto  ->Fill(sumTot);
        if(tot_neutronCut->IsInside(sumTot,pspTot)) tot_neutronHisto->Fill(sumTot);
        
        //Fill up histogram for sum of rods (MeVee)
        tot_cal_energyHisto ->Fill(sumCalTot);
        tot_cal_scatterHisto->Fill(sumCalTot, pspTot);
        if(tot_gammaCut  ->IsInside(sumTot,pspTot)){
            tot_cal_gammaHisto  ->Fill(sumCalTot);
            mul_erg_gammaHisto  ->Fill( sumCalTot , mulGam );
        }
        if(tot_neutronCut->IsInside(sumTot,pspTot)){
            tot_cal_neutronHisto->Fill(sumCalTot);
            mul_erg_neutronHisto->Fill( sumCalTot , mulNeu );
        }
        
    }
    
    //fill up teh time difference between the two SiPM
    //if( minTime0<1.e29 && minTime1<1.e29){ timeBetSipm->Fill( (minTime0-minTime1)*1.e12 ); } //ps
    
}

void PSDhisto::createAndWriteRootFile( string finalOutPath )
{
    //write root file
    TFile *run = new TFile(finalOutPath.c_str(),"recreate");
    run->SetCompressionSettings(3);
    
    gammaCut  ->Write();
    neutronCut->Write();
    
    tot_gammaCut  ->Write();
    tot_neutronCut->Write();
    
    timeBetSipm ->Write();
    digiActive  ->Write();
    posActive0 -> Write();
    posActive1 -> Write();
    
    tot_scatterHisto->Write();
    tot_energyHisto ->Write();
    tot_gammaHisto  ->Write();
    tot_neutronHisto->Write();
    tot_qRatioHisto ->Write();
    tot_qRatioErgHisto->Write();
    tot_erg0erg1Histo->Write();
    
    tot_cal_scatterHisto->Write();
    tot_cal_energyHisto ->Write();
    tot_cal_gammaHisto  ->Write();
    tot_cal_neutronHisto->Write();
    
    mul_erg_gammaHisto  ->Write();
    mul_erg_neutronHisto->Write();
    
    for(int i=0; i<numRodX; i++){
        for(int j=0; j<numRodY; j++){
            
            scatterHisto[i][j]->Write();
            energyHisto[i][j] ->Write();
            energyHisto0[i][j]->Write();
            energyHisto1[i][j]->Write();
            gammaHisto[i][j]  ->Write();
            neutronHisto[i][j]->Write();
            qRatioHisto[i][j] ->Write();
            
            cal_scatterHisto[i][j]->Write();
            cal_energyHisto[i][j] ->Write();
            cal_gammaHisto[i][j]  ->Write();
            cal_neutronHisto[i][j]->Write();
        }
    }
    
    run->Close();
}


