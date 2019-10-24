#define PSDanalyze_cxx

#include "PSDanalyze.h"
#include "PSDconstants.h"
#include "PSDhisto.h"
#include "PSDequalizer.h"

#include <TCutG.h>
#include <TTreeIndex.h>
#include <TFile.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void PSDanalyze::Loop(string fileName, PSDhisto histos, PSDequalizer myEqual)
{
    //read the root data file
    vector<vector<double>> myData;
    vector<double> tempVec;
    TFile *fin =  TFile::Open(fileName.c_str());
    TTree *tree = (TTree*)fin->Get("data");
    long long time; float energy; UInt_t channelID;
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("channelID", &channelID);
    Long64_t nentries = tree->GetEntriesFast();
    
    for ( Long64_t jentry = 0 ; jentry < nentries ; jentry++ ) {
        tree->GetEntry(jentry);
        tempVec.resize(0);
        tempVec.push_back(double(time)/1.e12);    //time in s
        tempVec.push_back(double(energy));        //integral
        tempVec.push_back(double(channelID));     //digi num
        myData.push_back(tempVec);
    }
    
    //sort the vector based on time
    std::sort (myData.begin(), myData.end());
    
    //Resize the time for new data file
    tNow = 0.0;
    
    //map
    int q_start_channels[32] = {18,20,24,30,33,39,42,47,19,22,28,32,35,40,46,49,0,5,15,2,9,59,51,55,3,7,4,27,58,57,53,63};
    int q_total_channels[32] = {16,21,25,29,34,37,43,45,17,23,26,31,36,41,44,48,1,12,8,13,38,61,52,56,10,14,6,11,62,50,54,60};
    //int q_start_channels[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
    //int q_total_channels[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
    
    //Set vectors that accumulate tail and total charge for each event
    vector< vector<double> > adcTail0 ( numRodX,vector<double>(numRodY,0.0) );
    vector< vector<double> > adcTot0  ( numRodX,vector<double>(numRodY,0.0) );
    vector< vector<double> > adcTail1 ( numRodX,vector<double>(numRodY,0.0) );
    vector< vector<double> > adcTot1  ( numRodX,vector<double>(numRodY,0.0) );
    vector< vector<double> > time0    ( numRodX,vector<double>(numRodY,0.0) );
    vector< vector<double> > time1    ( numRodX,vector<double>(numRodY,0.0) );
    
    //Start the loop
    for (Long64_t jentry=0 ; jentry<nentries ; jentry++)
    {
        //Check if code is still running
        if( jentry % 10 == 0 )cout<<"\r"<<100*(float)(jentry)/(float)nentries<<"% "<<flush;
        
        //get info
        tot = 0.0;tail = 0.0;
        timestamp = myData[jentry][0];
        erg       = myData[jentry][1];
        digiNum   = int(myData[jentry][2]);
        histos.fillDigi(myData[jentry][2]);
        
        //digi contains q tail and q total
        //diginum goes from 0 to 64 and some of them is tot, and some is start
        if     (                digiNum<64  ){SipmNum = 0; asicNum=0;                        }
        else if(digiNum>=64  && digiNum<64*2){SipmNum = 1; asicNum=0; digiNum = digiNum-64;  }
        else if(digiNum>=64*2&& digiNum<64*3){SipmNum = 0; asicNum=1; digiNum = digiNum-64*2;}
        else                                 {SipmNum = 1; asicNum=1; digiNum = digiNum-64*3;}
        /*
        if     (                digiNum<64  ){SipmNum = 0; asicNum=0;                        }
        else if(digiNum>=64  && digiNum<64*2){SipmNum = 0; asicNum=1; digiNum = digiNum-64;  }
        else if(digiNum>=64*2&& digiNum<64*3){SipmNum = 1; asicNum=0; digiNum = digiNum-64*2;}
        else                                 {SipmNum = 1; asicNum=1; digiNum = digiNum-64*3;}*/
        
        
        //find th channel number
        for(int k=0;k<32;k++){
            
            //this is a start integral
            if( digiNum == q_start_channels[k]){
                
                if(asicNum==0){chNum = k;   }
                else          {chNum = 32 + (31-k);}
                tail = erg; //tail is actualy start
                break;
            }
            //this is a total integral
            else if( digiNum == q_total_channels[k]){
                
                if(asicNum==0){chNum = k;   }
                else          {chNum = 32 + (31-k);}
                tot = erg;
                break;
            }
        }
        
        //determine x and y
        xpos = int(chNum) / 8 ; //multipication factor
        ypos = int(chNum) % 8 ; //remiander
        if(SipmNum==1){xpos = 7-xpos;} //mirrored one y axis because SiPM B is mirrored of SiPM A
        
        if(SipmNum==0){histos.fillPos0(double(xpos),double(ypos));}
        else{histos.fillPos1(double(xpos),double(ypos));}
        
        //if(tot>0.0 ){
        
        //pulses are still one event
        if( timestamp-tNow < timeOneEv ){
            
            if(SipmNum==0){
                adcTail0[xpos][ypos] += tail;
                adcTot0 [xpos][ypos] += tot;
                if(time0[xpos][ypos]==0.0){time0[xpos][ypos]=timestamp;} //take the first time channel detects signal
            }
            if(SipmNum==1){
                adcTail1[xpos][ypos] += tail;
                adcTot1 [xpos][ypos] += tot;
                if(time1[xpos][ypos]==0.0){time1[xpos][ypos]=timestamp;} //take the first time channel detects signal
            }
            
            cout << "(time,tot,tail):(" << timestamp <<","<< tot<<","<<tail << ") (digiNum,chNum):("<<digiNum<<","<< chNum << ") (SiPMnum,asicNum,xpos,ypos):("<< SipmNum <<","<<asicNum << ","<< xpos<<","<<ypos <<")"<< endl;
            
        }
        
        //this is a new event
        else{
            
            cout << "new event =========================" << endl;
            
            cout << "(time,tot,tail):(" << timestamp <<","<< tot<<","<<tail << ") (digiNum,chNum):("<<digiNum<<","<< chNum << ") (SiPMnum,asicNum,xpos,ypos):("<< SipmNum <<","<<asicNum << ","<< xpos<<","<<ypos <<")"<< endl;
            
            //process the old event
            histos.fillHistos(adcTail0,adcTot0,adcTail1,adcTot1,time0,time1);
            
            //set for new event
            tNow = timestamp;
            for (int i = 0; i < numRodX; i++){
                for (int j = 0; j < numRodY; j++){
                    adcTail0[i][j]= 0.0; adcTot0[i][j]= 0.0;
                    adcTail1[i][j]= 0.0; adcTot1[i][j]= 0.0;
                    time0[i][j]   = 0.0; time1[i][j]  = 0.0;
                }
            }
            
            //include this new event
            if(SipmNum==0){
                adcTail0[xpos][ypos] += tail;
                adcTot0 [xpos][ypos] += tot;
                if(time0[xpos][ypos]==0.0){time0[xpos][ypos]=timestamp;} //take the first time channel detects signal
            }
            if(SipmNum==1){
                adcTail1[xpos][ypos] += tail;
                adcTot1 [xpos][ypos] += tot;
                if(time1[xpos][ypos]==0.0){time1[xpos][ypos]=timestamp;} //take the first time channel detects signal
            }
            
        }//new event
        //}//good pulse
    }//loop of pulse
}

PSDanalyze::PSDanalyze(string filename)
{
    cout << "file name is: " << filename << endl;
}

PSDanalyze::~PSDanalyze()
{
    //empty
}





