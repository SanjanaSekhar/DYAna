#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>


using namespace std;

  float BinHisto(string filename, string histname){
    cout << "in class" << endl;
    TFile *F = TFile::Open(filename.c_str());
    TH1F *hist = (TH1F *)F->Get(histname.c_str());
    int binsx = hist->GetNbinsX();
    std::cout << "hi" << std::endl;
    std::cout << F << std::endl;
    for(int b = 1; b < binsx; b++){
      
      float valx = hist->GetBin(b);
      std::cout << valx << std::endl;
      return valx;

    }
}


  int BinHistoFunc(){
      
    cout <<"hi"<< endl;
  //  BinHistoFunc m;
 float valx  =   BinHisto("triggerRatePlots_RunE.root", "HTTrigNum");
    return 0;
  }


