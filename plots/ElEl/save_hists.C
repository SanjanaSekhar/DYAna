
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
#include "../../utils/Colors.h"

const int type = FLAG_ELECTRONS;
const int year = 2016;
const bool write_out = true;

char *fout_name = "ElEl/LQ_saved_hists.root";

void save_hists(){

    char year_str[80];
    sprintf(year_str, "y%i", year);

    setTDRStyle();
    init(year);
    init_indv_bkgs(year);

    float yLQ = 200.0;
    float m_LQ = 2000.;
    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300., 700. };

    TH1F *dy_pt = new TH1F("dy_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *dy_tautau_pt = new TH1F("dy_tautau_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *data_pt = new TH1F("data_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *diboson_pt = new TH1F("diboson_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *wt_pt = new TH1F("wt_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *QCD_pt = new TH1F("QCD_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *gg_pt = new TH1F("gg_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *LQu_pt = new TH1F("LQu_pt", "S-eu signal", n_pt_bins1, pt_bins1);
    TH1F *LQu_vec_pt = new TH1F("LQu_vec_pt", "V-eu signal", n_pt_bins1, pt_bins1);

    int n_xf_bins1 = 5;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.2, 0.5};
    TH1F *dy_xf = new TH1F("dy_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *dy_tautau_xf = new TH1F("dy_tautau_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *data_xf = new TH1F("data_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *diboson_xf = new TH1F("diboson_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *wt_xf = new TH1F("wt_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *QCD_xf = new TH1F("QCD_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *gg_xf = new TH1F("gg_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *LQu_xf = new TH1F("LQu_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *LQu_vec_xf = new TH1F("LQu_vec_xf", "dy signal", n_xf_bins1,  xf_bins1);

    int n_m_bins = 10;
    float mbin_base = 10.;
    //Float_t mbins1[] = {170.,200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1800., 2400.};
    
    Float_t mbins1[] = {500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1800., 2400.};
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, mbins1);
    TH1F *dy_m = new TH1F("dy_m", "dy Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *dy_tautau_m = new TH1F("dy_tautau_m", "dy no signal (qq, gluglu qbarqbar)", n_m_bins, mbins1);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, mbins1);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, mbins1);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, mbins1);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, mbins1);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, mbins1);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, mbins1);
    TH1F *LQu_m = new TH1F("LQu_m", "S-eu", n_m_bins, mbins1);
    TH1F *LQu_vec_m = new TH1F("LQu_vec_m", "V-eu", n_m_bins, mbins1);

    int n_cost_bins = 10;
    float cost_bin_size = 2./n_cost_bins;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *dy_cost = new TH1F("dy_cost", "dy Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1,1);
    TH1F *dy_tautau_cost = new TH1F("dy_tautau_cost", "dy no signal (qq, gluglu qbarqbar)", n_cost_bins, -1.,1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", n_cost_bins, -1,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1,1);
    TH1F *gg_cost = new TH1F("gg_cost", "QCD", n_cost_bins, -1,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1,1);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", n_cost_bins, -1,1);
    TH1F *LQu_cost = new TH1F("LQu_cost", "S-eu", n_cost_bins, -1,1);
    TH1F *LQu_vec_cost = new TH1F("LQu_vec_cost", "V-eu", n_cost_bins, -1,1);

    int n_phi_bins = 20;
    TH1F *data_phi = new TH1F("data_phi", "Data", n_phi_bins, -4.,4.);
    TH1F *dy_phi = new TH1F("dy_phi", "dy Signal (qqbar, qglu, qbarglu)", n_phi_bins, -4,4);
    TH1F *dy_tautau_phi = new TH1F("dy_tautau_phi", "dy no signal (qq, gluglu qbarqbar)", n_phi_bins, -4.,4.);
    TH1F *ttbar_phi = new TH1F("ttbar_phi", "TTbar Background", n_phi_bins, -4.,4.);
    TH1F *diboson_phi = new TH1F("diboson_phi", "DiBoson (WW, WZ,ZZ)", n_phi_bins, -4,4);
    TH1F *QCD_phi = new TH1F("QCD_phi", "QCD", n_phi_bins, -4,4);
    TH1F *gg_phi = new TH1F("gg_phi", "QCD", n_phi_bins, -4,4);
    TH1F *WJets_phi = new TH1F("WJets_phi", "WJets", n_phi_bins, -4,4);
    TH1F *wt_phi = new TH1F("wt_phi", "tw + #bar{t}w", n_phi_bins, -4,4);
    TH1F *LQu_phi = new TH1F("LQu_phi", "S-eu", n_phi_bins, -4,4);
    TH1F *LQu_vec_phi = new TH1F("LQu_vec_phi", "V-eu", n_phi_bins, -4,4);

    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *dy_rap = new TH1F("dy_rap", "dy Signal (qqbar, qglu, qbarglu)", n_rap_bins, -2.5,2.5);
    TH1F *dy_tautau_rap = new TH1F("dy_tautau_rap", "dy no signal (qq, gluglu qbarqbar)", n_rap_bins, -2.5,2.5);
    TH1F *ttbar_rap = new TH1F("ttbar_rap", "TTbar Background", n_rap_bins, -2.5,2.5);
    TH1F *diboson_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *QCD_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *gg_rap = new TH1F("gg_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *WJets_rap = new TH1F("WJets_rap", "WJets", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "tw + #bar{t}w", n_rap_bins, -2.5,2.5);
    TH1F *LQu_rap = new TH1F("LQu_rap", "S-eu", n_rap_bins, -2.5,2.5);
    TH1F *LQu_vec_rap = new TH1F("LQu_vec_rap", "V-eu", n_rap_bins, -2.5,2.5);

    float m_low = 500.;
    float m_high = 10000.;

    make_m_cost_pt_xf_hist(t_elel_data, data_m, data_cost, data_pt, data_xf, data_phi, data_rap, true, type,  year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_mc, dy_m, dy_cost, dy_pt, dy_xf, dy_phi, dy_rap, false, type, year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_tautau, dy_tautau_m, dy_tautau_cost, dy_tautau_pt, dy_tautau_xf, dy_tautau_phi, dy_tautau_rap, false, type,  year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, ttbar_phi, ttbar_rap, false, type,  year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_wt, wt_m, wt_cost, wt_pt, wt_xf, wt_phi, wt_rap, false, type,  year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_gamgam, gg_m, gg_cost, gg_pt, gg_xf, gg_phi, gg_rap, false, type,  year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, diboson_rap, false, type,   year, m_low, m_high, false, 0);
    make_m_cost_pt_xf_hist(t_elel_mc, LQu_m, LQu_cost, LQu_pt, LQu_xf, LQu_phi, LQu_rap, false, type, year, m_low, m_high, false, 1, yLQ, m_LQ);
    make_m_cost_pt_xf_hist(t_elel_mc, LQu_vec_m, LQu_vec_cost, LQu_vec_pt, LQu_vec_xf, LQu_vec_phi, LQu_vec_rap, false, type, year, m_low, m_high, false, 2, yLQ, m_LQ);
    symmetrize1d(gg_cost);

    //gg_cost->Scale(0.);
    //gg_m->Scale(0.);

    bool ss_qcd = false;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, QCD_phi, QCD_rap, 
            type, year, m_low, m_high, ss_qcd);


    /*
    QCD_m->Scale(0.6);
    QCD_cost->Scale(0.6);
    QCD_pt->Scale(0.6);
    QCD_xf->Scale(0.6);
    */


    printf("Data integral is %.2f \n", data_cost->Integral());
    printf("DY integral is %.2f \n", dy_cost->Integral());
    printf("ttbar integral is %.2f \n", ttbar_cost->Integral());







    TFile *fout = new TFile(fout_name, "UPDATE");
    fout->cd();
    fout->mkdir(year_str);
    fout->cd(year_str);

    data_m->Write();
    diboson_m->Write();
    QCD_m->Write();
    ttbar_m->Write();
    gg_m->Write();
    dy_tautau_m->Write();
    dy_m->Write();
    wt_m->Write();
    LQu_m->Write();
    LQu_vec_m->Write();

    data_cost->Write();
    diboson_cost->Write();
    QCD_cost->Write();
    ttbar_cost->Write();
    gg_cost->Write();
    dy_tautau_cost->Write();
    dy_cost->Write();
    wt_cost->Write();
    LQu_cost->Write();
    LQu_vec_cost->Write();


    data_rap->Write();
    diboson_rap->Write();
    QCD_rap->Write();
    ttbar_rap->Write();
    gg_rap->Write();
    dy_tautau_rap->Write();
    dy_rap->Write();
    wt_rap->Write();
    LQu_rap->Write();
    LQu_vec_rap->Write();

}

    
    
