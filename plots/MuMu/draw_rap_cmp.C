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



const int type = FLAG_MUONS;
const int year = 2017;
const bool write_out = true;
char *plot_dir = "Misc_plots/rap_check/";




void draw_rap_cmp(){


    
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);

    gROOT->SetBatch(1);
    TH1F *dummy = new TH1F("dummy", "dummy", 100, 0, 100);

    int n_rap_bins = 10;
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

    dy_tautau_rap->SetFillColor(tautau_c);

    dy_rap->SetFillColor(DY_c);

    ttbar_rap->SetFillColor(ttbar_c);


    wt_rap->SetFillColor(wt_c);

    diboson_rap->SetFillColor(diboson_c);

    QCD_rap->SetFillColor(qcd_c);

    gg_rap->SetFillColor(kOrange);


    for(int i=0; i<n_m_bins; i++){

        data_rap->Reset();
        dy_rap->Reset();
        dy_tautau_rap->Reset();
        ttbar_rap->Reset();
        wt_rap->Reset();
        gg_rap->Reset();
        diboson_rap->Reset();
        QCD_rap->Reset();




        float m_low = m_bins[i];
        float m_high = m_bins[i+1];;
        bool ss = false;

        make_m_cost_pt_xf_hist(t_mumu_data, dummy, dummy, dummy, dummy, dummy, data_rap, true, type,  year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_mc, dummy, dummy, dummy, dummy, dummy,  dy_rap,              false, type,   year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_tautau, dummy, dummy, dummy, dummy, dummy,  dy_tautau_rap, false, type,  year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_ttbar, dummy, dummy, dummy, dummy, dummy, ttbar_rap, false, type,  year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_wt, dummy, dummy, dummy, dummy, dummy, wt_rap, false, type,  year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_gamgam, dummy, dummy, dummy, dummy, dummy, gg_rap, false, type,  year, m_low, m_high);
        make_m_cost_pt_xf_hist(t_mumu_diboson, dummy, dummy, dummy, dummy, dummy, diboson_rap, false, type,   year, m_low, m_high);



        bool ss_qcd = false;
        make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, dummy, dummy, dummy, dummy, dummy, QCD_rap, type, year, m_low, m_high, ss_qcd);





        setHistError(QCD_rap, qcd_sys_unc);

        setHistError(dy_rap, dy_sys_unc);

        setHistError(diboson_rap, diboson_sys_unc);

        setHistError(ttbar_rap, top_sys_unc);

        setHistError(wt_rap, top_sys_unc);


        setHistError(gg_rap, gam_sys_unc);



        THStack *rap_stack = new THStack("rap_stack", "DiElectron Rapidity Distribution: Data vs MC; y");
        rap_stack->Add(diboson_rap);
        rap_stack->Add(QCD_rap);
        rap_stack->Add(wt_rap);
        rap_stack->Add(ttbar_rap);
        rap_stack->Add(gg_rap);
        rap_stack->Add(dy_tautau_rap);
        rap_stack->Add(dy_rap);



        gStyle->SetLegendBorderSize(0);
        float x_size = 0.3;
        float y_size = 0.3;
        TLegend *leg1 = new TLegend(x_size, y_size);
        char leg_title[100];
        sprintf(leg_title, "Muons %i, Mass %.0f - %.0f GeV", year, m_low, m_high);
        leg1->SetHeader(leg_title, "C");
        leg1->AddEntry(data_rap, "data", "p");
        leg1->AddEntry(dy_rap, "DY Signal", "f");
        leg1->AddEntry(dy_tautau_rap, "DY #rightarrow #tau#tau", "f");
        leg1->AddEntry(gg_rap, "#gamma#gamma #rightarrow #mu#mu", "f");
        leg1->AddEntry(ttbar_rap, "t#bar{t}", "f");
        leg1->AddEntry(wt_rap, "tW + #bar{t}W", "f");
        leg1->AddEntry(QCD_rap, "QCD + WJets", "f");
        leg1->AddEntry(diboson_rap, "WW + WZ + ZZ", "f");

     
        //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
        TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
        TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
        int iPeriod = 4; 
        writeExtraText = true;
        char plt_file[100], y_ax_label[100];



        bool logy = true;

        bool logx = false;
        bool draw_sys_uncs = true;
        float ratio_range = 0.5;

        char plot_label[100];

        sprintf(plot_label, "Muons %.0f-%.0f GeV",m_low, m_high);


        sprintf(y_ax_label, "Events/%.2f", rap_bin_size);
        std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg1, "rap", "dimuon Y",y_ax_label, plot_label, -1., logy, logx, draw_sys_uncs, ratio_range);
        CMS_lumi(p_rap, year, 33);
        sprintf(plt_file, "%sMuMu%i_mbin%i_rap_cmp.png", plot_dir, year % 2000, i);
        if(write_out) c_rap->Print(plt_file);

    }

}

    
    
