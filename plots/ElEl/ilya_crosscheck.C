
#define STAND_ALONE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/HistMaker.C"
#include "../../utils/root_files.h"
#include "../../utils/PlotUtils.C"



const int type = FLAG_ELECTRONS;
const int year = 2016;

double count_evts(TTree *t1, double m_low, double m_high, bool is_data = false, int year = 2016){
    TempMaker tm(t1, is_data, year);
    tm.do_electrons = true;
    tm.setup();
    int nEvents=0;
    double tot_evts = 0.;

    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        bool pass = (tm.m >= m_low && tm.m <= m_high);
        if(pass){
            tot_evts += tm.getEvtWeight();
        }
    }
    return tot_evts;
}




void ilya_crosscheck(){

    TFile * f_data = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Analysis_ntuples/ElEl16_data_mlow_nov1.root");
    TTree * t_data  = (TTree *) gDirectory->Get("T_sig");
    TFile * f_mc = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Analysis_ntuples/ElEl16_dy_mlow_nov4.root");
    TTree * t_mc  = (TTree *) gDirectory->Get("T_sig");


    double m_low = 60.;
    double m_high = 120.;

    double n_data = count_evts(t_data, m_low, m_high, true, year);
    double n_mc = count_evts(t_mc, m_low, m_high, false, year);

    printf("data events %f  MC events %f \n", n_data, n_mc);

}





