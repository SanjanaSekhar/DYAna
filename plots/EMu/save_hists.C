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
#include "../../utils/Colors.h"

int year = 2018;
const bool write_out = true;
//char *plot_label = "e#mu Control Region";
char *plot_label = "";

char *fout_name = "EMu/saved_hists.root";




void save_hists(){

    char year_str[80];
    sprintf(year_str, "y%i", year);

    printf("Year is %i \n", year);
    setTDRStyle();

    init_emu(year);
    init_emu_indv_bkgs(year);
    setup_all_SFs(year);

    int n_m_bins = 8;
    float mbin_base = 10.;
    Float_t mbins1[] = {170.,200., 250., 300., 350., 400., 500., 700., 1000.};

    TH1F *data_m = new TH1F("data_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *ttbar_m = new TH1F("ttbar_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *diboson_m = new TH1F("diboson_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *wt_m = new TH1F("wt_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *dy_m = new TH1F("dy_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *QCD_m = new TH1F("QCD_m", "MC Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);

    float pt_bin_size  = 10.;
    TH1F *data_pt = new TH1F("data_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *wt_pt = new TH1F("wt_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *dy_pt = new TH1F("dy_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC Signal (qqbar, qglu, qbarglu)", 30, 0, 300);

    int n_cost_bins = 8;
    float cost_bin_size = 2./8;
    TH1F *data_cost = new TH1F("data_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *diboson_cost = new TH1F("diboson_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *wt_cost = new TH1F("wt_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *dy_cost = new TH1F("dy_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "MC Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1, 1);


    int n_rap_bins = 20;
	float rap_bin_low = -2.4;
    float rap_bin_high = 2.4;
    float rap_bin_size = (rap_bin_high - rap_bin_low) / n_rap_bins;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, rap_bin_low,rap_bin_high);
    TH1F *ttbar_rap = new TH1F("ttbar_rap", "TTbar Background", n_rap_bins, rap_bin_low,rap_bin_high);
    TH1F *diboson_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, rap_bin_low,rap_bin_high);
    TH1F *wt_rap = new TH1F("wt_rap", "QCD", n_rap_bins, rap_bin_low,rap_bin_high);
    TH1F *dy_rap = new TH1F("dy_rap", "QCD", n_rap_bins, rap_bin_low,rap_bin_high);
    TH1F *QCD_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, rap_bin_low,rap_bin_high);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100.);

    Double_t m_low = 170;
    Double_t m_high = 10000;

    bool ss = false;

    bool do_emu_cost_rw = false;
    make_emu_m_cost_pt_rap_hist(t_emu_data, data_m, data_cost, data_pt, data_rap, true,  year, m_low, m_high, ss);
    make_emu_m_cost_pt_rap_hist(t_emu_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_diboson, diboson_m, diboson_cost, diboson_pt, diboson_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_wt, wt_m, wt_cost, wt_pt, wt_rap, false,  year, m_low, m_high, ss, do_emu_cost_rw);
    make_emu_m_cost_pt_rap_hist(t_emu_dy, dy_m, dy_cost, dy_pt, dy_rap, false,  year, m_low, m_high, ss);



    bool fakes_reweight = true;
    bool fakes_sys_errors = false;
    Fakerate_est_emu(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam, t_emu_QCD_contam, QCD_m,  QCD_cost, QCD_pt, QCD_rap, FLAG_MUONS, year, m_low, m_high, fakes_reweight, fakes_sys_errors);


    symmetrize1d(QCD_cost);
    symmetrize1d(ttbar_cost);
    //symmetrize1d(diboson_cost);
    //symmetrize1d(wt_cost);
    //symmetrize1d(dy_cost);


    Double_t data_count = data_m->Integral();
    Double_t fake_count = QCD_m->Integral();
    Double_t mc_count = ttbar_m->Integral() + diboson_m->Integral() + wt_m->Integral() + dy_m->Integral();


    Double_t fake_unc = qcd_sys_unc * fake_count;
    Double_t mc_unc = (ttbar_m->Integral() + wt_m->Integral()) * top_sys_unc + diboson_m->Integral() * diboson_sys_unc + dy_m->Integral()*dy_sys_unc;

    printf("Data count %.0f +/- %.0f \n", data_count, sqrt(data_count));
    printf("MC count %.0f +/- %0.f \n", mc_count, mc_unc);
    printf("TTbar count %.0f  (frac %.2f) \n", ttbar_m->Integral(), ttbar_m->Integral()/(mc_count + fake_count));
    printf("WT count %.0f  (frac %.2f) \n", wt_m->Integral(), wt_m->Integral()/(mc_count + fake_count));
    printf("Diboson count %.0f (frac %.2f ) \n", diboson_m->Integral(), diboson_m->Integral()/(mc_count + fake_count));
    printf("Fake count %.0f +/- %.0f (frac %.2f) \n", fake_count, fake_unc, fake_count / (mc_count + fake_count));
    Double_t ratio = data_count / (mc_count + fake_count);
    Double_t unc = sqrt( (data_count/(mc_count + fake_count)/(mc_count + fake_count)) +  
                          pow(data_count/(mc_count + fake_count)/(mc_count + fake_count), 2) * (fake_unc*fake_unc + mc_unc*mc_unc));
    printf("Ratio is %1.3f +/- %1.3f \n", ratio, unc);

    TH1F *total_mc = (TH1F *) diboson_cost->Clone("total_mc");
    total_mc->Add(QCD_cost);
    total_mc->Add(dy_cost);
    total_mc->Add(wt_cost);
    total_mc->Add(ttbar_cost);


    Double_t data_B = data_cost->Integral(1,n_cost_bins/2);
    Double_t data_F = data_cost->Integral(n_cost_bins/2 + 1, n_cost_bins);
    Double_t data_AFB = (data_F - data_B)/(data_F+data_B);
   
    Double_t data_dAFB = AFB_counting_unc(data_F, data_B, sqrt(data_F), sqrt(data_B));
    printf("F %.0f, B %.0f \n", data_F, data_B);
    printf("AFB %.3f +/- %.3f \n", data_AFB, data_dAFB);



    Double_t total_mc_dF, total_mc_dB;
    Double_t total_mc_B = total_mc->IntegralAndError(1,n_cost_bins/2, total_mc_dB);
    Double_t total_mc_F = total_mc->IntegralAndError(n_cost_bins/2 + 1, n_cost_bins, total_mc_dF);
    Double_t total_mc_AFB = (total_mc_F - total_mc_B)/(total_mc_F+total_mc_B);
   
    Double_t total_mc_dAFB = AFB_counting_unc(total_mc_F, total_mc_B, total_mc_dF, total_mc_dB);
    printf("total_mc F %.0f, B %.0f \n", total_mc_F, total_mc_B);
    printf("total_mc AFB %.3f +/- %.3f \n", total_mc_AFB, total_mc_dAFB);


    Double_t diboson_dF, diboson_dB;
    Double_t diboson_B = diboson_cost->IntegralAndError(1,n_cost_bins/2, diboson_dB);
    Double_t diboson_F = diboson_cost->IntegralAndError(n_cost_bins/2 + 1, n_cost_bins, diboson_dF);
    Double_t diboson_AFB = (diboson_F - diboson_B)/(diboson_F+diboson_B);
   
    Double_t diboson_dAFB = AFB_counting_unc(diboson_F, diboson_B, diboson_dF, diboson_dB);
    printf("Diboson F %.0f, B %.0f \n", diboson_F, diboson_B);
    printf("Diboson AFB %.3f +/- %.3f \n", diboson_AFB, diboson_dAFB);




    TFile *fout = new TFile(fout_name, "UPDATE");
    fout->cd();
    fout->mkdir(year_str);
    fout->cd(year_str);

    data_m->Write();
    diboson_m->Write();
    QCD_m->Write();
    ttbar_m->Write();
    dy_m->Write();
    wt_m->Write();
    
    data_cost->Write();
    diboson_cost->Write();
    QCD_cost->Write();
    ttbar_cost->Write();
    dy_cost->Write();
    wt_cost->Write();


    data_rap->Write();
    diboson_rap->Write();
    QCD_rap->Write();
    ttbar_rap->Write();
    dy_rap->Write();
    wt_rap->Write();


}








