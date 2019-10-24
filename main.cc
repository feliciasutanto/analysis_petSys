// F. Sutanto (Oct 2019)
// Analysis for PET sys daq format
// Class structure is largely borrowed from J.Nattress analysis for CAEN CoMPASS

// Recommended usage:
// g++ main.cc PSDhisto.cc PSDanalyze.cc PSDequalizer.cc `root-config --libs --cflags` -o main -O3 -ftree-vectorize
// ./main "path/to/file"

#include "PSDanalyze.h"
#include "PSDhisto.h"
#include "PSDequalizer.h"
#include <ctime>

int main(int argc, char** argv)
{
    //Get the common name of file
    cout << "What's the file's name?" << endl;
    string commonName; cin >> commonName;
    
    //Get the number of files we're processing
    //cout << "How many files are you going to process?" << endl;
    //int numOfFilesToProcess; cin >> numOfFilesToProcess;
    
    //start clock to estimate execution time
    int start_s=clock();
    
    //Create histograms to collect data
    PSDhisto histos;
    histos.createHistos();
    histos.MakeCuts();
    
    //Read equalizer data
    PSDequalizer myEqual;
    myEqual.ReadDataEqualizer();
    
    //Start to perform analysis
    //for( int fileNum=0; fileNum < numOfFilesToProcess; fileNum++ )
    //{
    //Get the name of teh file to proecess
    stringstream ss;
    //ss << argv[1] << "/" << commonName << fileNum << ".root";
    ss << argv[1] << "/" << commonName;
    string fileName = ss.str();
    
    //Start analyzing the data for each pulse
    PSDanalyze data(fileName);
    data.Loop(fileName, histos, myEqual);
    
    //print out the total measurement time
    cout << endl << "Time total: " << data.tNow <<" s"<< endl;
    //}
    
    //Write resuts to a root file
    cout << "Writing output root file" << endl;
    stringstream sss;
    sss << argv[1] << "/o.root";
    histos.createAndWriteRootFile( sss.str() );
    
    //Stop clock and show total execution time
    int stop_s=clock();
    cout << "Total execution time: " <<
    (stop_s-start_s)/double(CLOCKS_PER_SEC) <<" s"<< endl;
    
    return 0;
}
