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



void draw_roch_cmp(){
    int year =2018;

    init(year);
    setTDRStyle();




    TH1F *h1_cost = new TH1F("data_cost", "Data", 40, -1.,1.);
    TH1F *h1_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 60, 140, 2000);
    TH1F *h1_pt = new TH1F("mc_pt", "MC signal", 40, 0, 1000);
    TH1F *h1_xf = new TH1F("msdjf_pt", "MC signal", 40, 0, 1000);

    TH1F *h2_cost = new TH1F("h2_cost", "Data", 40, -1.,1.);
    TH1F *h2_m = new TH1F("h2_m", "Data Dimuon Mass Distribution", 60, 140, 2000);
    TH1F *h2_pt = new TH1F("h2_pt", "MC signal", 40, 0, 1000);
    TH1F *h2_xf = new TH1F("hx_pt", "MC signal", 40, 0, 1000);

    TH1F *h_def_ratio = new TH1F("d", "Roch default ratio", 40, 0.5, 1.5);
    TH1F *h_alt_ratio = new TH1F("a", "Roch alt ratio", 40, 0.5, 1.5);

    bool is_data = true;
    int type = FLAG_MUONS;
    make_m_cost_pt_xf_hist(t_mumu_data, h1_m, h1_cost, h1_pt, h1_xf, is_data, type, true, year);
    make_m_cost_pt_xf_hist(t_mumu_data, h2_m, h2_cost, h2_pt, h2_xf, is_data, type, false, year);

    //make_pt_cmp(t1, h1_cost, h2_cost);
    //TCanvas *c1 = new TCanvas("c1", "", 200, 10, 900, 700);
    //h_def_ratio->Draw("hist");
    //TCanvas *c2 = new TCanvas("c2", "", 200, 10, 900, 700);
    //h_alt_ratio->Draw("hist");
    //make_ratio_plot("MuMu_mc_rocc_alt_cost_cmp.pdf", h1_cost, "RC def ",h2_cost, "RC Alt", "OFF/On", "Cos(#theta)", false);

    printf("ON integral is %.2f. OFF integral is %.2f \n", h1_m->Integral(), h2_m->Integral());

    make_ratio_plot("MuMu_data_roch_cost_cmp.pdf", h2_cost, "RC OFF ",h1_cost, "RC ON", "OFF/On", "Cos(#theta)", false);
    make_ratio_plot("MuMu_data_roch_m_cmp.pdf", h2_m, "RC OFF ",h1_m, "RC ON", "OFF/On", "M (GeV)", true);
}

