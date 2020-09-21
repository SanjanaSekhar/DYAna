#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"
#include <iostream>


void LQ_consistency_check(){

    
        bool ss = false;
        bool use_xF =false;
        //bool use_LQ_denom=true;
        bool draw_muons = false;
        bool draw_electrons = true;
        const string sys_label = "";
        
        //char *plot_dir = "Paper_plots/template_plots";
        char *plot_dir = "Misc_plots/LQ_consistency_check";
        int i_start = 1; 
        int i_end = 4;
        float x_start = 0.75;
        float x_end = 0.9;
        float y_start = 0.75;
        float y_end = 0.9;
       
       // std::cout << "enter m_LQ:";        
       // std::cin >> m_LQ; 
        gStyle->SetOptStat(0);
        gROOT->SetBatch(1);

         Double_t m_LQ=1500.;
        printf("=========================\n m_LQ = %f, draw_muons = %d, draw_electrons = %d \n=========================\n",m_LQ,draw_muons,draw_electrons );
         for (int year = 2016; year<=2018; year++)
        {

        init(year);
        setup_all_SFs(year);
            
            char title[100];

            char mu_title[100], el_title[100];
            char mu_fname1[100], mu_fname2[100], mu_fname3[100], el_fname1[100], el_fname2[100], el_fname3[100];

            sprintf(mu_title, "Muons, m_LQ=%0.1f TeV, year=%i",m_LQ/1000,year);
            sprintf(el_title, "Electrons, m_LQ=%0.1f TeV, year=%i",m_LQ/1000,year);
             sprintf(mu_fname1, "%s/Mu%i_dy_old_new.png", plot_dir, year%2000);

            auto h_mumu_sym = new TH3F(title, "Symmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_sym->SetDirectory(0);
            
            sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
            auto h_mumu_asym = new TH3F(title, "Asymmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_asym->SetDirectory(0);
            //--------------------------------------
            sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
            auto h_mumu_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_alpha->SetDirectory(0);
             sprintf(title, "mumu%i_LQpure_u%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_u = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_u->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_u%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_u = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_u->SetDirectory(0);
            sprintf(title, "mumu%i_LQpure_d%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_d = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_d->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_d%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_d = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_d->SetDirectory(0);
            //------------------------------------------
            
            bool old = true;
            gen_mc_template(t_mumu_mc, h_mumu_sym, h_mumu_asym, h_mumu_alpha,h_mumu_LQpure_u, h_mumu_LQint_u,h_mumu_LQpure_d, h_mumu_LQint_d, 
                year, m_LQ, FLAG_MUONS, use_xF, old, "");

            auto h1_mumu_sym = convert3d(h_mumu_sym);
            auto h1_mumu_asym = convert3d(h_mumu_asym);

            int n_1d_bins = n_lq_m_bins*(std::round(std::ceil(n_y_bins/2.) * n_cost_bins + std::floor(n_y_bins/2.) * (n_cost_bins-2)));

            sprintf(title, "mumu%i_fpl%s", year%2000, sys_label.c_str());
            auto h1_mumu_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
            h1_mumu_pl->SetDirectory(0);
            sprintf(title, "mumu%i_fmn%s", year%2000, sys_label.c_str());
            auto h1_mumu_mn = new TH1F(title, "minus template of DY", n_1d_bins, 0, n_1d_bins);
            h1_mumu_mn->SetDirectory(0);
            sprintf(title, "mumu%i_dy%s", year%2000, sys_label.c_str());
            auto h1_mumu_dy = new TH1F(title, "template of DY", n_1d_bins, 0, n_1d_bins);
            h1_mumu_dy->SetDirectory(0);
            
            make_pl_mn_templates(h1_mumu_sym, h1_mumu_asym, h1_mumu_pl, h1_mumu_mn);

            double norm = 3./8.;
            double afb = 0.6;

            h1_mumu_dy->Add(h1_mumu_pl,h1_mumu_mn,(norm+afb),(norm-afb));

            auto h_mumu_dy_new = new TH3F(title, "new dy template",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_dy_new->SetDirectory(0);

            old = false;
            gen_mc_template(t_mumu_mc, h_mumu_dy_new, h_mumu_asym, h_mumu_alpha,h_mumu_LQpure_u, h_mumu_LQint_u,h_mumu_LQpure_d, h_mumu_LQint_d, 
                year, m_LQ, FLAG_MUONS, use_xF, old, "");

            auto h1_mumu_dy_new = convert3d(h_mumu_dy_new);

             h1_mumu_dy->SetLineColor(kRed);
             h1_mumu_dy->SetLineWidth(2);
             h1_mumu_dy_new->SetLineColor(kBlue);
             h1_mumu_dy_new->SetLineWidth(2);

             TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_dy_new->SetTitle(mu_title); 
            h1_mumu_dy_new->Draw("hist");
            h1_mumu_dy->Draw("hist same");

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_mumu_dy, "DY Template - old", "l");
            leg1->AddEntry(h1_mumu_dy_new, "DY Template - new", "l");
            leg1->Draw();
            c_mumu1->Print(mu_fname1);
            delete c_mumu1;
        }
    }