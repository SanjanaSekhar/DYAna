#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"

TH3F *h_elel_asym, *h_elel_sym, *h_elel_alpha, *h_elel_LQpure,*h_elel_LQint,*h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd, *h_elel_gam;
TH3F *h_elel_mc_count, *h_elel_sym_count;
TH3F *h_mumu_asym, *h_mumu_sym, *h_mumu_alpha, *h_mumu_LQpure,*h_mumu_LQint,*h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd, *h_mumu_gam;

Double_t m_LQ;

void LQ_sys_check(){
        std::cout << "enter m_LQ:";        
        std::cin >> m_LQ;
        gStyle->SetOptStat(0);
        gROOT->SetBatch(1);
    
        int year = 2018;
        init(year);
        char *plot_dir = "Misc_plots";
        char *sys = "_REFAC";
        bool do_bkg = true;
        bool do_electrons = true;
        bool do_muons = false;
        int i = 5;
        setup_all_SFs(year);

        string sys_up = string(sys) + string("Up");
        string sys_down = string(sys) + string("Down");

        Double_t alpha_denom = amc_alpha[i];
        double m_low = m_bins[i];
        double m_high = m_bins[i+1];
        double afb = 0.6;


        char mu_fname1[100],  el_fname1[100];

        sprintf(mu_fname1, "%s/MuMu%i_LQ%s_chk.png", plot_dir, year, sys);
        sprintf(el_fname1, "%s/ElEl%i_LQ%s_chk.png", plot_dir, year, sys);

        bool use_xf = false;

        int n_var1_bins = n_y_bins;
        float *var1_bins = y_bins;
        if(use_xf){
            n_var1_bins = n_xf_bins;
            var1_bins = xf_bins;
        }

        TH3F * h_elel_bkg = new TH3F("elel_bkg", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_bkg_up = new TH3F("elel_bkg_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_bkg_down = new TH3F("elel_bkg_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_plain = new TH3F("elel_plain", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_sys_up = new TH3F("elel_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_sys_down = new TH3F("elel_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);

        TH3F * h_mumu_bkg = new TH3F("mumu_bkg", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_bkg_up = new TH3F("mumu_bkg_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_bkg_down = new TH3F("mumu_bkg_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_plain = new TH3F("mumu_plain", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_sys_up = new TH3F("mumu_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_sys_down = new TH3F("mumu_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);



        bool ss = false;




        TTree *elel_ts[3] = {t_elel_ttbar, t_elel_wt, t_elel_diboson};
        TTree *mumu_ts[3] = {t_mumu_ttbar, t_mumu_wt, t_mumu_diboson};
        char mu_title[100], el_title[100];

        TH1F *h1_elel_bkg, *h1_mumu_bkg, *h1_elel_bkg_up, *h1_elel_bkg_down, *h1_mumu_bkg_up, *h1_mumu_bkg_down;

        if(do_muons){
            printf("Making mumu temps \n");
            one_mc_template(t_mumu_mc, afb, h_mumu_plain, year, m_LQ, FLAG_MUONS, use_xf, "");
            one_mc_template(t_mumu_mc, afb, h_mumu_sys_up, year, m_LQ, FLAG_MUONS, use_xf, sys_up);
            one_mc_template(t_mumu_mc,  afb, h_mumu_sys_down, year, m_LQ, FLAG_MUONS, use_xf, sys_down);
            TH1F *h1_mumu_plain = convert3d(h_mumu_plain);
            TH1F *h1_mumu_sys_up = convert3d(h_mumu_sys_up);
            TH1F *h1_mumu_sys_down = convert3d(h_mumu_sys_down);

            h1_mumu_plain->SetLineColor(kBlack);
            h1_mumu_plain->SetLineWidth(2);


            h1_mumu_sys_up->SetLineColor(kBlue);
            h1_mumu_sys_up->SetLineWidth(2);
            h1_mumu_sys_down->SetLineColor(kGreen+3);
            h1_mumu_sys_down->SetLineWidth(2);
            printf("MuMu: nom %.0f, up %.0f, down %.0f \n", h_mumu_plain->Integral(), h_mumu_sys_up->Integral(), h_mumu_sys_down->Integral());

            if(do_bkg){

                gen_combined_background_template(3, mumu_ts, h_mumu_bkg, year, FLAG_MUONS,  ss, use_xf,  "");
                gen_combined_background_template(3, mumu_ts, h_mumu_bkg_up, year, FLAG_MUONS,  ss, use_xf,  sys_up);
                gen_combined_background_template(3, mumu_ts, h_mumu_bkg_down, year, FLAG_MUONS,  ss, use_xf, sys_down);
                h1_mumu_bkg = convert3d(h_mumu_bkg);
                h1_mumu_bkg_up = convert3d(h_mumu_bkg_up);
                h1_mumu_bkg_down = convert3d(h_mumu_bkg_down);

                h1_mumu_bkg->SetLineColor(kRed);
                h1_mumu_bkg->SetLineWidth(2);
                h1_mumu_bkg_up->SetLineColor(kMagenta);
                h1_mumu_bkg_up->SetLineWidth(2);
                h1_mumu_bkg_down->SetLineColor(kRed-7);
                h1_mumu_bkg_down->SetLineWidth(2);
                printf("MuMu Bkg: nom %.0f, up %.0f, down %.0f \n", h_mumu_bkg->Integral(), h_mumu_bkg_up->Integral(), h_mumu_bkg_down->Integral());
            }



            sprintf(mu_title, "Muons: %s", sys);
            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Muons", 200, 10, 900, 700);
            h1_mumu_plain->SetTitle(mu_title);
            h1_mumu_plain->Draw("hist");
            h1_mumu_sys_up->Draw("hist same");
            h1_mumu_sys_down->Draw("hist same");
            if(do_bkg){
                h1_mumu_bkg->Draw("hist same");
                h1_mumu_bkg_up->Draw("hist same");
                h1_mumu_bkg_down->Draw("hist same");
            }

            TLegend *leg1 = new TLegend(0.15, 0.15);
            leg1->AddEntry(h1_mumu_plain, "Nominal Template", "l");
            leg1->AddEntry(h1_mumu_sys_up, "Sys Up Template", "l");
            leg1->AddEntry(h1_mumu_sys_down, "Sys Down Template", "l");
            if(do_bkg){
                leg1->AddEntry(h1_mumu_bkg, "Nominal Bkg Template", "l");
                leg1->AddEntry(h1_mumu_bkg_up, "Sys Up Bkg Template", "l");
                leg1->AddEntry(h1_mumu_bkg_down, "Sys Down Bkg Template", "l");
            }
            leg1->Draw();

            c_mumu1->Print(mu_fname1);
        }


        if(do_electrons){
            printf("Making elel temps \n");

            one_mc_template(t_elel_mc,  afb, h_elel_plain, year, m_LQ, FLAG_ELECTRONS, use_xf, "");
            one_mc_template(t_elel_mc, afb, h_elel_sys_up, year, m_LQ, FLAG_ELECTRONS, use_xf, sys_up);
            one_mc_template(t_elel_mc, afb, h_elel_sys_down, year, m_LQ, FLAG_ELECTRONS, use_xf, sys_down);
            TH1F *h1_elel_plain = convert3d(h_elel_plain);
            TH1F *h1_elel_sys_up = convert3d(h_elel_sys_up);
            TH1F *h1_elel_sys_down = convert3d(h_elel_sys_down);

            h1_elel_plain->SetLineColor(kBlack);
            h1_elel_plain->SetLineWidth(2);


            h1_elel_sys_up->SetLineColor(kBlue);
            h1_elel_sys_up->SetLineWidth(2);
            h1_elel_sys_down->SetLineColor(kGreen+3);
            h1_elel_sys_down->SetLineWidth(2);
            printf("elel: nom %.0f, up %.0f, down %.0f \n", h_elel_plain->Integral(), h_elel_sys_up->Integral(), h_elel_sys_down->Integral());

            if(do_bkg){

                gen_combined_background_template(3, elel_ts, h_elel_bkg, year, FLAG_ELECTRONS,  ss, use_xf,  "");
                gen_combined_background_template(3, elel_ts, h_elel_bkg_up, year, FLAG_ELECTRONS,  ss, use_xf,  sys_up);
                gen_combined_background_template(3, elel_ts, h_elel_bkg_down, year, FLAG_ELECTRONS,  ss, use_xf, sys_down);
                h1_elel_bkg = convert3d(h_elel_bkg);
                h1_elel_bkg_up = convert3d(h_elel_bkg_up);
                h1_elel_bkg_down = convert3d(h_elel_bkg_down);

                h1_elel_bkg->SetLineColor(kRed);
                h1_elel_bkg->SetLineWidth(2);
                h1_elel_bkg_up->SetLineColor(kMagenta);
                h1_elel_bkg_up->SetLineWidth(2);
                h1_elel_bkg_down->SetLineColor(kRed-7);
                h1_elel_bkg_down->SetLineWidth(2);
                printf("elel Bkg: nom %.0f, up %.0f, down %.0f \n", h_elel_bkg->Integral(), h_elel_bkg_up->Integral(), h_elel_bkg_down->Integral());
            }



            sprintf(el_title, "Electrons: %s", sys);
            TCanvas *c_elel1 = new TCanvas("c_elel", "Electrons", 200, 10, 900, 700);
            h1_elel_plain->SetTitle(el_title);
            h1_elel_plain->Draw("hist");
            h1_elel_sys_up->Draw("hist same");
            h1_elel_sys_down->Draw("hist same");
            if(do_bkg){
                h1_elel_bkg->Draw("hist same");
                h1_elel_bkg_up->Draw("hist same");
                h1_elel_bkg_down->Draw("hist same");
            }

            TLegend *leg1 = new TLegend(0.15, 0.15);
            leg1->AddEntry(h1_elel_plain, "Nominal Template", "l");
            leg1->AddEntry(h1_elel_sys_up, "Sys Up Template", "l");
            leg1->AddEntry(h1_elel_sys_down, "Sys Down Template", "l");
            if(do_bkg){
                leg1->AddEntry(h1_elel_bkg, "Nominal Bkg Template", "l");
                leg1->AddEntry(h1_elel_bkg_up, "Sys Up Bkg Template", "l");
                leg1->AddEntry(h1_elel_bkg_down, "Sys Down Bkg Template", "l");
            }
            leg1->Draw();

            c_elel1->Print(el_fname1);
        
        }






}


