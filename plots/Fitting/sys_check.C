#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/TemplateMaker_systematics.C"
#include "../../utils/root_files.h"
#include "utils.C"

TH2F *h_elel_asym, *h_elel_sym, *h_elel_alpha, *h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd, *h_elel_gam;
TH2F *h_elel_mc_count, *h_elel_sym_count;
TH2F *h_mumu_asym, *h_mumu_sym, *h_mumu_alpha, *h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd, *h_mumu_gam;



void sys_check(){
        gStyle->SetOptStat(0);
    
        int year = 2016;
        init(year);
        char *plot_dir = "template_plots";
        char *sys = "RENORM";
        int i = 0;
        //setup_all_SFs(year);

        string sys_up = string(sys) + string("Up");
        string sys_down = string(sys) + string("Down");

        Double_t alpha_denom = amc_alpha[i];
        double m_low = m_bins[i];
        double m_high = m_bins[i+1];


        char mu_fname1[100],  el_fname1[100];

        sprintf(mu_fname1, "%s/MuMu%i_M%.0f_%s_chk.png", plot_dir, year, m_low, sys);
        sprintf(el_fname1, "%s/ElEl%i_M%.0f_%s_chk.png", plot_dir, year, m_low, sys);


        TH2F * h_elel_plain = new TH2F("elel_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_elel_sys_up = new TH2F("elel_up", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_elel_sys_down = new TH2F("elel_down", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);

        TH2F * h_mumu_plain = new TH2F("mumu_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_mumu_sys_up = new TH2F("mumu_up", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_mumu_sys_down = new TH2F("mumu_down", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);



        bool ss = false;
        bool do_RC = true;




        TTree *elel_ts[1] = {t_elel_mc};
        TTree *mumu_ts[1] = {t_mumu_mc};
        printf("Making elel temps \n");
        gen_combined_background_template(1, elel_ts, h_elel_plain, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss, "");
        gen_combined_background_template(1, elel_ts, h_elel_sys_up, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss, sys_up);
        gen_combined_background_template(1, elel_ts, h_elel_sys_down, year, m_low, m_high, FLAG_ELECTRONS,  do_RC, ss, sys_down);

        printf("Making mumu temps \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_plain, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, "");
        gen_combined_background_template(1, mumu_ts, h_mumu_sys_up, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_up);
        gen_combined_background_template(1, mumu_ts, h_mumu_sys_down, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_down);


        TH1F *h1_mumu_plain = convert2d(h_mumu_plain);
        TH1F *h1_mumu_sys_up = convert2d(h_mumu_sys_up);
        TH1F *h1_mumu_sys_down = convert2d(h_mumu_sys_down);


        TH1F *h1_elel_plain = convert2d(h_elel_plain);
        TH1F *h1_elel_sys_up = convert2d(h_elel_sys_up);
        TH1F *h1_elel_sys_down = convert2d(h_elel_sys_down);

        h1_mumu_plain->SetLineColor(kBlack);
        h1_mumu_plain->SetLineWidth(2);

        h1_mumu_sys_up->SetLineColor(kBlue);
        h1_mumu_sys_up->SetLineWidth(2);
        h1_mumu_sys_down->SetLineColor(kGreen+3);
        h1_mumu_sys_down->SetLineWidth(2);


        h1_elel_plain->SetLineColor(kBlack);
        h1_elel_plain->SetLineWidth(2);

        h1_elel_sys_up->SetLineColor(kBlue);
        h1_elel_sys_up->SetLineWidth(2);
        h1_elel_sys_down->SetLineColor(kGreen+3);
        h1_elel_sys_down->SetLineWidth(2);


        char mu_title[100], el_title[100];

        sprintf(mu_title, "Muons: %s", sys);
        sprintf(el_title, "Electrons: %s", sys);

        TCanvas *c_mumu1 = new TCanvas("c_mumu", "Muons", 200, 10, 900, 700);
        h1_mumu_plain->SetTitle(mu_title);
        h1_mumu_plain->Draw("hist");
        h1_mumu_sys_up->Draw("hist same");
        h1_mumu_sys_down->Draw("hist same");

        TLegend *leg1 = new TLegend(0.15, 0.15);
        leg1->AddEntry(h1_mumu_plain, "Nominal Template", "l");
        leg1->AddEntry(h1_mumu_sys_up, "Sys Up Template", "l");
        leg1->AddEntry(h1_mumu_sys_down, "Sys Down Template", "l");
        leg1->Draw();

        c_mumu1->Print(mu_fname1);




        TCanvas *c_elel1 = new TCanvas("c_elel", "Electrons", 200, 10, 900, 700);
        h1_elel_plain->SetTitle(mu_title);
        h1_elel_plain->Draw("hist");
        h1_elel_sys_up->Draw("hist same");
        h1_elel_sys_down->Draw("hist same");

        TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
        leg2->Draw();


        c_elel1->Print(el_fname1);

}




