
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
const int year = 2018;
const bool write_out = true;

char *fout_name = "ElEl/saved_hists.root";

void save_hists(){

    char year_str[80];
    sprintf(year_str, "y%i", year);

    setTDRStyle();
    init(year);
    init_indv_bkgs(year);

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

    int n_m_bins = 17;
    float mbin_base = 10.;
    Float_t mbins1[] = {170.,200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1800., 2400.};
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, mbins1);
    TH1F *dy_m = new TH1F("dy_m", "dy Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *dy_tautau_m = new TH1F("dy_tautau_m", "dy no signal (qq, gluglu qbarqbar)", n_m_bins, mbins1);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, mbins1);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, mbins1);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, mbins1);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, mbins1);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, mbins1);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, mbins1);

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

    dy_tautau_cost->SetFillColor(tautau_c);
    dy_tautau_m->SetFillColor(tautau_c);
    dy_tautau_pt->SetFillColor(tautau_c);
    dy_tautau_xf->SetFillColor(tautau_c);
    dy_tautau_phi->SetFillColor(tautau_c);
    dy_tautau_rap->SetFillColor(tautau_c);

    dy_cost->SetFillColor(DY_c);
    dy_m->SetFillColor(DY_c);
    dy_pt->SetFillColor(DY_c);
    dy_xf->SetFillColor(DY_c);
    dy_phi->SetFillColor(DY_c);
    dy_rap->SetFillColor(DY_c);

    ttbar_cost->SetFillColor(ttbar_c);
    ttbar_m->SetFillColor(ttbar_c);
    ttbar_pt->SetFillColor(ttbar_c);
    ttbar_xf->SetFillColor(ttbar_c);
    ttbar_phi->SetFillColor(ttbar_c);
    ttbar_rap->SetFillColor(ttbar_c);

    wt_cost->SetFillColor(wt_c);
    wt_m->SetFillColor(wt_c);
    wt_pt->SetFillColor(wt_c);
    wt_xf->SetFillColor(wt_c);
    wt_phi->SetFillColor(wt_c);
    wt_rap->SetFillColor(wt_c);

    diboson_cost->SetFillColor(diboson_c);
    diboson_m->SetFillColor(diboson_c);
    diboson_pt->SetFillColor(diboson_c);
    diboson_xf->SetFillColor(diboson_c);
    diboson_phi->SetFillColor(diboson_c);
    diboson_rap->SetFillColor(diboson_c);

    QCD_cost->SetFillColor(qcd_c);
    QCD_m->SetFillColor(qcd_c);
    QCD_pt->SetFillColor(qcd_c);
    QCD_xf->SetFillColor(qcd_c);
    QCD_phi->SetFillColor(qcd_c);
    QCD_rap->SetFillColor(qcd_c);

    gg_cost->SetFillColor(gamgam_c);
    gg_m->SetFillColor(gamgam_c);
    gg_pt->SetFillColor(gamgam_c);
    gg_xf->SetFillColor(gamgam_c);
    gg_phi->SetFillColor(gamgam_c);
    gg_rap->SetFillColor(gamgam_c);

    float m_low = 170.;
    float m_high = 10000.;

    make_m_cost_pt_xf_hist(t_elel_data, data_m, data_cost, data_pt, data_xf, data_phi, data_rap, true, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_mc, dy_m, dy_cost, dy_pt, dy_xf, dy_phi, dy_rap, false, type,   year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_tautau, dy_tautau_m, dy_tautau_cost, dy_tautau_pt, dy_tautau_xf, dy_tautau_phi, dy_tautau_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, ttbar_phi, ttbar_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_wt, wt_m, wt_cost, wt_pt, wt_xf, wt_phi, wt_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_gamgam, gg_m, gg_cost, gg_pt, gg_xf, gg_phi, gg_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_elel_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, diboson_rap, false, type,   year, m_low, m_high);


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






    setHistError(QCD_m, qcd_sys_unc);
    setHistError(QCD_cost, qcd_sys_unc);
    setHistError(QCD_pt, qcd_sys_unc);
    setHistError(QCD_xf, qcd_sys_unc);
    setHistError(QCD_rap, qcd_sys_unc);
    setHistError(QCD_phi, qcd_sys_unc);

    setHistError(dy_cost, dy_sys_unc);
    setHistError(dy_pt, dy_sys_unc);
    setHistError(dy_xf, dy_sys_unc);
    setHistError(dy_rap, dy_sys_unc);
    setHistError(dy_phi, dy_sys_unc);

    setHistError(diboson_cost, diboson_sys_unc);
    setHistError(diboson_pt, diboson_sys_unc);
    setHistError(diboson_xf, diboson_sys_unc);
    setHistError(diboson_rap, diboson_sys_unc);
    setHistError(diboson_phi, diboson_sys_unc);

    setHistError(ttbar_cost, top_sys_unc);
    setHistError(ttbar_pt, top_sys_unc);
    setHistError(ttbar_xf, top_sys_unc);
    setHistError(ttbar_rap, top_sys_unc);
    setHistError(ttbar_phi, top_sys_unc);

    setHistError(wt_cost, top_sys_unc);
    setHistError(wt_pt, top_sys_unc);
    setHistError(wt_xf, top_sys_unc);
    setHistError(wt_rap, top_sys_unc);
    setHistError(wt_phi, top_sys_unc);


    setHistError(gg_cost, gam_sys_unc);
    setHistError(gg_pt, gam_sys_unc);
    setHistError(gg_xf, gam_sys_unc);
    setHistError(gg_rap, gam_sys_unc);
    setHistError(gg_phi, gam_sys_unc);



    setHistMassDepError(dy_m);
    setHistMassDepError(diboson_m);
    setHistMassDepError(ttbar_m);
    setHistMassDepError(wt_m);
    setHistMassDepError(wt_m);
    setHistMassDepError(gg_m);


    binwidth_normalize(data_m, mbin_base);
    binwidth_normalize(diboson_m, mbin_base);
    binwidth_normalize(QCD_m, mbin_base);
    binwidth_normalize(wt_m, mbin_base);
    binwidth_normalize(ttbar_m, mbin_base);
    binwidth_normalize(gg_m, mbin_base);
    binwidth_normalize(dy_tautau_m, mbin_base);
    binwidth_normalize(dy_m, mbin_base);

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
    
    data_cost->Write();
    diboson_cost->Write();
    QCD_cost->Write();
    ttbar_cost->Write();
    gg_cost->Write();
    dy_tautau_cost->Write();
    dy_cost->Write();
    wt_cost->Write();


    data_rap->Write();
    diboson_rap->Write();
    QCD_rap->Write();
    ttbar_rap->Write();
    gg_rap->Write();
    dy_tautau_rap->Write();
    dy_rap->Write();
    wt_rap->Write();

}

    
    
