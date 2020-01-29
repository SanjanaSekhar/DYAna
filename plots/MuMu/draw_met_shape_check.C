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
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"



const int type = FLAG_MUONS;
const int year = 2017;
const bool write_out = false;
char *plot_dir = "Paper_plots/";


void make_cost_xf_hist(TTree *t1, TH1F *h_cost,  TH1F *h_xf, bool is_data=false, int flag1 = FLAG_MUONS, 
        int year = 2016, Double_t met_min = 50., Double_t m_low = 150., Double_t m_high = 9999999.){
    //read event data
        TempMaker tm(t1, is_data, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;

        tm.setup();
        int nEvents=0;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_min  && tm.has_no_bjets && tm.not_cosmic;
            //bool good_eta = (fabs(tm.el1_eta) > 1.5 && fabs(tm.el2_eta) < 1.4) || (fabs(tm.el2_eta) > 1.5 && fabs(tm.el1_eta) < 1.4);
            //pass = pass && good_eta;


            if(pass){
                nEvents++;
                tm.doCorrections();
                tm.getEvtWeight();

                h_xf->Fill(tm.xF, tm.evt_weight);
                h_cost->Fill(tm.cost, tm.evt_weight);

            }


        
        }
        printf("Selected %i events \n", nEvents);


    tm.finish();
}


void draw_met_shape_check(){

    setTDRStyle();
    init(year);
    init_indv_bkgs(year);

    int n_xf_bins1 = 5;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.2, 0.5};
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *ttbar_low_xf = new TH1F("ttbar_low_xf", "MC signal", n_xf_bins1,  xf_bins1);
    TH1F *ttbar_high_xf = new TH1F("ttbar_high_xf", "MC signal", n_xf_bins1,  xf_bins1);


    int num_cost_bins = 10;
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *ttbar_low_cost = new TH1F("ttbar_low_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *ttbar_high_cost = new TH1F("ttbar_high_cost", "TTbar Background", n_cost_bins, -1.,1.);


    double m_low = 150.;
    double m_high = 100000.;
    make_cost_xf_hist(t_mumu_mc, ttbar_low_cost, ttbar_low_xf, false, type, year, 47., m_low, m_high);
    make_cost_xf_hist(t_mumu_mc, ttbar_cost, ttbar_xf, false, type, year, 50., m_low, m_high);
    make_cost_xf_hist(t_mumu_mc, ttbar_high_cost, ttbar_high_xf, false, type, year, 53., m_low, m_high);

    make_ratio_plot("met_up_cost_compare.pdf", ttbar_high_cost, "MET up", ttbar_cost, "MET nom", "Ratio", "cos(#theta)", false, false);
    make_ratio_plot("met_down_cost_compare.pdf", ttbar_low_cost, "MET down", ttbar_cost, "MET nom", "Ratio", "cos(#theta)", false, false);


    make_ratio_plot("met_up_xf_compare.pdf", ttbar_high_xf, "MET up", ttbar_xf, "MET nom", "Ratio", "x_F", true, false);
    make_ratio_plot("met_down_xf_compare.pdf", ttbar_low_xf, "MET down", ttbar_xf, "MET nom", "Ratio", "x_F", true, false);
}
