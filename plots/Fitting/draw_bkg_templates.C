#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "utils.C"
//#include "../../utils/TemplateMaker_systematics.C"
//#include "../../utils/root_files.h"
#include "../../analyze/combine/TemplateUtils.h"



void draw_bkg_templates(){
        gStyle->SetOptStat(0);
    
        int year = 2018;
        init(year);
        char *plot_dir = "Paper_plots/template_plots";
        //setup_all_SFs(year);
        string sys_label = "";

        for(int i=0; i<n_m_bins; i++){
            Double_t alpha_denom = amc_alpha[i];
            double m_low = m_bins[i];
            double m_high = m_bins[i+1];

            auto h_mumu_dy_gg = new TH2F("h_mumu_dy_gg", "dy_gg",
                    n_xf_bins, xf_bins, n_cost_bins, cost_bins);
            h_mumu_dy_gg->SetDirectory(0);
            auto h_mumu_back = new TH2F("h_mumu_bkg", "dy_gg",
                    n_xf_bins, xf_bins, n_cost_bins, cost_bins);
            h_mumu_back->SetDirectory(0);
            auto h_elel_dy_gg = new TH2F("h_elel_bkg", "dy_gg",
                    n_xf_bins, xf_bins, n_cost_bins, cost_bins);
            h_elel_dy_gg->SetDirectory(0);
            auto h_elel_back = new TH2F("h_elel_bkg", "dy_gg",
                    n_xf_bins, xf_bins, n_cost_bins, cost_bins);
            h_elel_back->SetDirectory(0);

            bool ss = false;
            bool do_RC = true;

            TTree *mumu_ts[1] = {t_mumu_back};
            gen_combined_background_template(1, mumu_ts, h_mumu_back, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_label);
            mumu_ts[0] = t_mumu_nosig;
            gen_combined_background_template(1, mumu_ts, h_mumu_dy_gg, year, m_low, m_high, FLAG_MUONS,  do_RC, ss, sys_label);


            TTree *elel_ts[1] = {t_elel_back};
            gen_combined_background_template(1, elel_ts, h_elel_back, year, m_low, m_high, FLAG_ELECTRONS, do_RC,ss, sys_label);
            elel_ts[0] = t_elel_nosig;
            gen_combined_background_template(1, elel_ts, h_elel_dy_gg, year, m_low, m_high, FLAG_ELECTRONS, do_RC,ss, sys_label);

            h_elel_dy_gg->Print("all");
            h_mumu_dy_gg->Print("all");


            symmetrize2d(h_elel_back);
            symmetrize2d(h_mumu_back);

            symmetrize2d(h_elel_dy_gg);
            symmetrize2d(h_mumu_dy_gg);





            char mu_title[100], el_title[100];
            sprintf(mu_title, "Muons Mass %.0f to %.0f", m_low, m_high);
            sprintf(el_title, "Electrons Mass %.0f to %.0f", m_low, m_high);

            char mu_fname1[100], mu_fname2[100], el_fname1[100], el_fname2[100];

            sprintf(mu_fname1, "%s/MuMu%i_M%.0f_bkg_temps.png", plot_dir, year, m_low);
            sprintf(el_fname1, "%s/ElEl%i_M%.0f_bkg_temps.png", plot_dir, year, m_low);


            auto h1_mumu_back = convert2d(h_mumu_back);
            auto h1_mumu_dy_gg = convert2d(h_mumu_dy_gg);
            


            h1_mumu_back->SetLineColor(kGreen +3);
            h1_mumu_dy_gg->SetLineColor(kRed +1);

            h1_mumu_back->SetLineWidth(2);
            h1_mumu_dy_gg->SetLineWidth(2);

            h1_mumu_dy_gg->SetMaximum(std::max(h1_mumu_dy_gg->GetMaximum(), h1_mumu_back->GetMaximum()) *1.2);

            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_dy_gg->SetTitle(mu_title); 
            h1_mumu_dy_gg->Draw("hist");
            h1_mumu_back->Draw("hist same ");

            TLegend *leg1 = new TLegend(0.15, 0.15);
            leg1->AddEntry(h1_mumu_dy_gg, "DY No Asym Template", "l");
            leg1->AddEntry(h1_mumu_back, "Backgrounds Template", "l");
            leg1->Draw();
            c_mumu1->Print(mu_fname1);

            auto h1_elel_back = convert2d(h_elel_back);
            auto h1_elel_dy_gg = convert2d(h_elel_dy_gg);
            


            h1_elel_back->SetLineColor(kGreen +3);
            h1_elel_dy_gg->SetLineColor(kRed +1);

            h1_elel_back->SetLineWidth(2);
            h1_elel_dy_gg->SetLineWidth(2);

            h1_elel_dy_gg->SetMaximum(std::max(h1_elel_dy_gg->GetMaximum(), h1_elel_back->GetMaximum()) *1.2);

            TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
            h1_elel_dy_gg->SetTitle(el_title); 
            h1_elel_dy_gg->Draw("hist");
            h1_elel_back->Draw("hist same ");

            TLegend *leg2 = new TLegend(0.15, 0.15);
            leg2->AddEntry(h1_elel_dy_gg, "DY No Asym Template", "l");
            leg2->AddEntry(h1_elel_back, "Backgrounds Template", "l");
            leg2->Draw();
            c_elel1->Print(el_fname1);


        }
}




