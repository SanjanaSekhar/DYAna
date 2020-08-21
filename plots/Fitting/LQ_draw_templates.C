#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"
#include <iostream>


void LQ_draw_templates(){

    
        bool ss = false;
        bool use_xF =false;
        //bool use_LQ_denom=true;
        bool draw_muons = true;
        bool draw_electrons = false;
        const string sys_label = "";
        
        //char *plot_dir = "Paper_plots/template_plots";
        char *plot_dir = "Misc_plots/template_plots2";
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

        TCanvas *c16_el_lqpall = new TCanvas("c_el_lqpall", "Histograms", 200, 10, 900, 700);
        TLegend *leg_lqpall = new TLegend(x_start, y_start, x_end, y_end);
        TCanvas *c16_el_lqiall = new TCanvas("c_el_lqiall", "Histograms", 200, 10, 900, 700);
        TLegend *leg_lqiall = new TLegend(x_start, y_start, x_end, y_end);
        TCanvas *c16_mu_lqpall = new TCanvas("c_mu_lqpall", "Histograms", 200, 10, 900, 700);
        //TLegend *leg_lqpall = new TLegend(x_start, y_start, x_end, y_end);
        TCanvas *c16_mu_lqiall = new TCanvas("c_mu_lqiall", "Histograms", 200, 10, 900, 700);
        //TLegend *leg_lqiall = new TLegend(x_start, y_start, x_end, y_end);

        char el_all_lqp16[100], el_all_lqi16[100];
        sprintf(el_all_lqp16,"%s/Electrons_LQpure_all16.png",plot_dir);
        sprintf(el_all_lqi16,"%s/Electrons_LQint_all16.png",plot_dir);
        char el_title_all_1[100], el_title_all_2[100];
        sprintf(el_title_all_1, "Electrons: LQpure, year=2016");
        sprintf(el_title_all_2, "Electrons: LQint");

        char mu_all_lqp16[100], mu_all_lqi16[100];
        sprintf(mu_all_lqp16,"%s/Muons_LQpure_all16.png",plot_dir);
        sprintf(mu_all_lqi16,"%s/Muons_LQint_all16.png",plot_dir);
        char mu_title_all_1[100], mu_title_all_2[100];
        sprintf(mu_title_all_1, "Muons: LQpure, year=2016");
        sprintf(mu_title_all_2, "Muons: LQint");

     // for(int i=i_start;i<=i_end;i++)
      // {
        Double_t m_LQ=1500.;
        printf("=========================\n m_LQ = %f, draw_muons = %d, draw_electrons = %d \n=========================\n",m_LQ,draw_muons,draw_electrons );

        TH1F *h16_mumu_LQpure_u, *h16_mumu_LQint_u, *h16_mumu_LQpure_d, *h16_mumu_LQint_d;
        TH1F *h17_mumu_LQpure_u, *h17_mumu_LQint_u, *h17_mumu_LQpure_d, *h17_mumu_LQint_d;
        TH1F *h18_mumu_LQpure_u, *h18_mumu_LQint_u, *h18_mumu_LQpure_d, *h18_mumu_LQint_d;

        TH1F *h16_elel_LQpure_u, *h16_elel_LQint_u, *h16_elel_LQpure_d, *h16_elel_LQint_d;
        TH1F *h17_elel_LQpure_u, *h17_elel_LQint_u, *h17_elel_LQpure_d, *h17_elel_LQint_d;
        TH1F *h18_elel_LQpure_u, *h18_elel_LQint_u, *h18_elel_LQpure_d, *h18_elel_LQint_d;

        
        //int year = 2017;
        for (int year = 2016; year<=2018; year++)
        {

        init(year);
        setup_all_SFs(year);
            
            char title[100];

            char mu_title[100], el_title[100];
            char mu_fname1[100], mu_fname2[100], mu_fname3[100], el_fname1[100], el_fname2[100], el_fname3[100];

            sprintf(mu_title, "Muons, m_LQ=%0.1f TeV, year=%i",m_LQ/1000,year);
            sprintf(el_title, "Electrons, m_LQ=%0.1f TeV, year=%i",m_LQ/1000,year);

            if(draw_muons){

            auto h_mumu_sym = new TH3F(title, "Symmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_sym->SetDirectory(0);
            sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
            auto h_mumu_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_alpha->SetDirectory(0);
            sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
            auto h_mumu_asym = new TH3F(title, "Asymmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_asym->SetDirectory(0);
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
            
            

            gen_mc_template(t_mumu_mc, h_mumu_sym, h_mumu_asym, h_mumu_alpha,h_mumu_LQpure_u, h_mumu_LQint_u,h_mumu_LQpure_d, h_mumu_LQint_d, 
                year, m_LQ, FLAG_MUONS, use_xF,"");

            sprintf(mu_fname1, "%s/Mu%i_MC_SM_m%i.png", plot_dir, year%2000,int(m_LQ));
            sprintf(mu_fname2, "%s/Mu%i_MC_LQpure_m%i.png", plot_dir, year%2000, int(m_LQ));
            sprintf(mu_fname3, "%s/Mu%i_MC_LQint_m%i.png", plot_dir, year%2000, int(m_LQ));

            auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
            auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
            h_mumu_pl.Scale(0.5);
            h_mumu_mn.Scale(0.5);
        //    Double_t alpha= 0.05;
          //  double norm = 3./4./(2.+alpha);
           // h_mumu_alpha->Scale(norm);
           
            auto h1_mumu_pl = convert3d(&h_mumu_pl);
            auto h1_mumu_mn = convert3d(&h_mumu_mn);
            auto h1_mumu_alpha = convert3d(h_mumu_alpha);
            auto h1_mumu_sym = convert3d(h_mumu_sym);
            auto h1_mumu_asym = convert3d(h_mumu_asym);
            auto h1_mumu_LQpure_u = convert3d(h_mumu_LQpure_u);
            auto h1_mumu_LQint_u = convert3d(h_mumu_LQint_u);
            auto h1_mumu_LQpure_d = convert3d(h_mumu_LQpure_d);
            auto h1_mumu_LQint_d = convert3d(h_mumu_LQint_d);


            h1_mumu_alpha->SetLineColor(kGreen +3);
            h1_mumu_sym->SetLineColor(kBlue);
            h1_mumu_asym->SetLineColor(kRed+1);
            h1_mumu_pl->SetLineColor(kOrange +7);
            h1_mumu_mn->SetLineColor(kRed + 1);
            h1_mumu_LQpure_u->SetLineColor(kOrange +9);
            h1_mumu_LQint_u->SetLineColor(kGreen +3);
            h1_mumu_LQpure_d->SetLineColor(kBlue);
            h1_mumu_LQint_d->SetLineColor(kRed );

            h1_mumu_alpha->SetLineWidth(2);
            h1_mumu_sym->SetLineWidth(2);
            h1_mumu_asym->SetLineWidth(2);
            h1_mumu_pl->SetLineWidth(2);
            h1_mumu_mn->SetLineWidth(2);
            h1_mumu_LQpure_u->SetLineWidth(2);
            h1_mumu_LQint_u->SetLineWidth(2);
            h1_mumu_LQpure_d->SetLineWidth(2);
            h1_mumu_LQint_d->SetLineWidth(2);


            if(year==2016)
            {
                h16_mumu_LQpure_u = h1_mumu_LQpure_u;
                h16_mumu_LQint_u = h1_mumu_LQint_u;
                h16_mumu_LQpure_d = h1_mumu_LQpure_d;
                h16_mumu_LQint_d = h1_mumu_LQint_d;
            }
            if(year==2017)
            {
                h17_mumu_LQpure_u = h1_mumu_LQpure_u;
                h17_mumu_LQint_u = h1_mumu_LQint_u;
                h17_mumu_LQpure_d = h1_mumu_LQpure_d;
                h17_mumu_LQint_d = h1_mumu_LQint_d;
            }
            if(year==2018)
            {
                h18_mumu_LQpure_u = h1_mumu_LQpure_u;
                h18_mumu_LQint_u = h1_mumu_LQint_u;
                h18_mumu_LQpure_d = h1_mumu_LQpure_d;
                h18_mumu_LQint_d = h1_mumu_LQint_d;
            }

           // h1_mumu_asym->SetMaximum(h1_mumu_sym->GetMaximum()*1.2);
            //h1_mumu_LQpure->SetMaximum(h1_mumu_LQint->GetMaximum()*1.2);
            
            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_alpha->SetTitle(mu_title); 
            //h1_mumu_asym->Draw("hist");
            //h1_mumu_sym->Draw("hist same ");
            //h1_mumu_pl->Draw("hist");
            h1_mumu_alpha->Draw("hist");
            //h1_mumu_mn->Draw("hist same");
            

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            //leg1->AddEntry(h1_mumu_asym, "Asym Template", "l");
            //leg1->AddEntry(h1_mumu_sym, "Sym Template", "l");
            leg1->AddEntry(h1_mumu_pl, "Plus Template", "l");
            leg1->AddEntry(h1_mumu_mn, "Minus Template", "l");
            leg1->AddEntry(h1_mumu_alpha, "alpha Template", "l");
            leg1->Draw();
            //c_mumu1->Print(mu_fname1);
            delete c_mumu1;
        
            TCanvas *c_mumu2 = new TCanvas("c_mumu2", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure_u->SetTitle(mu_title);
            h1_mumu_LQpure_u->Draw("hist");
            h1_mumu_LQpure_d->Draw("hist same");

            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_mumu_LQpure_u,"u-LQpure Template","l");
            leg2->AddEntry(h1_mumu_LQpure_d,"d-LQpure Template","l");
            leg2->Draw();

            c_mumu2->Print(mu_fname2);
            delete c_mumu2;

            TCanvas *c_mumu3 = new TCanvas("c_mumu3", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQint_u->SetTitle(mu_title);
            h1_mumu_LQint_u->Draw("hist");
            h1_mumu_LQint_d->Draw("hist same");

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h1_mumu_LQint_u,"u-LQint Template","l");
            leg3->AddEntry(h1_mumu_LQint_d,"d-LQint Template","l");
            leg3->Draw();
            
           c_mumu3->Print(mu_fname3);
            delete c_mumu3;

        }

        if(draw_electrons){

            auto h_elel_sym = new TH3F(title, "Symmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_sym->SetDirectory(0);
            sprintf(title, "elel%i_alpha%s", year %2000, sys_label.c_str());
            auto h_elel_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_alpha->SetDirectory(0);
            sprintf(title, "elel%i_asym%s", year %2000, sys_label.c_str());
            auto h_elel_asym = new TH3F(title, "Asymmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_asym->SetDirectory(0);
            sprintf(title, "elel%i_LQpure_u%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_u = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_u->SetDirectory(0);
            sprintf(title, "elel%i_LQint_u%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_u = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_u->SetDirectory(0);
            sprintf(title, "elel%i_LQpure_d%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_d = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_d->SetDirectory(0);
            sprintf(title, "elel%i_LQint_d%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_d = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_d->SetDirectory(0);
            

            gen_mc_template(t_elel_mc, h_elel_sym, h_elel_asym, h_elel_alpha,h_elel_LQpure_u, h_elel_LQint_u,h_elel_LQpure_d, h_elel_LQint_d, 
                year, m_LQ, FLAG_ELECTRONS, use_xF, "");

            sprintf(el_fname1, "%s/El%i_MC_SM_m%i.png", plot_dir, year%2000,int(m_LQ));
            sprintf(el_fname2, "%s/El%i_MC_LQpure_m%i.png", plot_dir, year%2000,int(m_LQ));
            sprintf(el_fname3, "%s/El%i_MC_LQint_m%i.png", plot_dir, year%2000,int(m_LQ));

            auto h_elel_pl = *h_elel_sym + *h_elel_asym;
            auto h_elel_mn = *h_elel_sym - *h_elel_asym;
            
            h_elel_pl.Scale(0.5);
            h_elel_mn.Scale(0.5);
            Double_t alpha= 0.05;
            double norm = 3./4./(2.+alpha);
            h_elel_alpha->Scale(norm);

            auto h1_elel_pl = convert3d(&h_elel_pl);
            auto h1_elel_mn = convert3d(&h_elel_mn);
            auto h1_elel_alpha = convert3d(h_elel_alpha);
            auto h1_elel_sym = convert3d(h_elel_sym);
            auto h1_elel_asym = convert3d(h_elel_asym);
            auto h1_elel_LQpure_u = convert3d(h_elel_LQpure_u);
            auto h1_elel_LQint_u = convert3d(h_elel_LQint_u);
            auto h1_elel_LQpure_d = convert3d(h_elel_LQpure_d);
            auto h1_elel_LQint_d = convert3d(h_elel_LQint_d);

            h1_elel_alpha->SetLineColor(kGreen +3);
            h1_elel_sym->SetLineColor(kBlue);
            h1_elel_asym->SetLineColor(kRed+1);
            h1_elel_pl->SetLineColor(kOrange +7);
            h1_elel_mn->SetLineColor(kRed + 1);
            h1_elel_LQpure_u->SetLineColor(kOrange +9);
            h1_elel_LQint_u->SetLineColor(kGreen +3);
            h1_elel_LQpure_d->SetLineColor(kBlue);
            h1_elel_LQint_d->SetLineColor(kRed);

            h1_elel_alpha->SetLineWidth(2);
            h1_elel_sym->SetLineWidth(2);
            h1_elel_asym->SetLineWidth(2);
            h1_elel_pl->SetLineWidth(2);
            h1_elel_mn->SetLineWidth(2);
            h1_elel_LQpure_u->SetLineWidth(2);
            h1_elel_LQint_u->SetLineWidth(2);
            h1_elel_LQpure_d->SetLineWidth(2);
            h1_elel_LQint_d->SetLineWidth(2);

            if(year==2016)
            {
                h16_elel_LQpure_u = h1_elel_LQpure_u;
                h16_elel_LQint_u = h1_elel_LQint_u;
               h16_elel_LQpure_d = h1_elel_LQpure_d;
               h16_elel_LQint_d = h1_elel_LQint_d;
            }
            if(year==2017)
            {
                h17_elel_LQpure_u = h1_elel_LQpure_u;
                h17_elel_LQint_u = h1_elel_LQint_u;
                h17_elel_LQpure_d = h1_elel_LQpure_d;
               h17_elel_LQint_d = h1_elel_LQint_d;
            }
            if(year==2018)
            {
                h18_elel_LQpure_u = h1_elel_LQpure_u;
                h18_elel_LQint_u = h1_elel_LQint_u;
                h18_elel_LQpure_d = h1_elel_LQpure_d;
                h18_elel_LQint_d = h1_elel_LQint_d;
            }
            // h1_elel_asym->SetMaxielm(h1_elel_sym->GetMaxielm()*1.2);
            //h1_elel_LQpure->SetMaxielm(h1_elel_LQint->GetMaxielm()*1.2);
            
            TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
            h1_elel_pl->SetTitle(el_title); 
            //h1_elel_asym->Draw("hist");
            //h1_elel_sym->Draw("hist same ");
            h1_elel_pl->Draw("hist");
            h1_elel_alpha->Draw("hist same");
            h1_elel_mn->Draw("hist same");
            

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            //leg1->AddEntry(h1_elel_asym, "Asym Template", "l");
            //leg1->AddEntry(h1_elel_sym, "Sym Template", "l");
            leg1->AddEntry(h1_elel_pl, "Plus Template", "l");
            leg1->AddEntry(h1_elel_mn, "Minus Template", "l");
            leg1->AddEntry(h1_elel_alpha, "alpha Template", "l");
            leg1->Draw();
            c_elel1->Print(el_fname1);
            delete c_elel1;
        
            TCanvas *c_elel2 = new TCanvas("c_elel2", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure_u->SetTitle(el_title);
            h1_elel_LQpure_u->Draw("hist");
            h1_elel_LQpure_d->Draw("hist same");

            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_elel_LQpure_u,"u-LQpure Template","l");
            leg2->AddEntry(h1_elel_LQpure_d,"d-LQpure Template","l");
            leg2->Draw();

            c_elel2->Print(el_fname2);
            delete c_elel2;

            TCanvas *c_elel3 = new TCanvas("c_elel3", "Histograms", 200, 10, 900, 700);
            h1_elel_LQint_u->SetTitle(el_title);
            h1_elel_LQint_u->Draw("hist");
            h1_elel_LQint_d->Draw("hist same");

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h1_elel_LQint_u,"u-LQint Template","l");
            leg3->AddEntry(h1_elel_LQint_d,"d-LQint Template","l");
            leg3->Draw();
            
            c_elel3->Print(el_fname3);
            delete c_elel3;
         
        }
    
    }
    char mu_title_1[100], el_title_1[100];
    sprintf(mu_title_1, "Muons: m_LQ = %i",int(m_LQ));
    sprintf(el_title_1, "Electrons: m_LQ = %i",int(m_LQ));
    char mu_1[100], mu_2[100], el_1[100], el_2[100];
    sprintf(mu_1, "%s/Muons_LQpure_m%i.png",plot_dir, int(m_LQ));
    sprintf(mu_2, "%s/Muons_LQint_m%i.png",plot_dir,int(m_LQ));
    sprintf(el_1, "%s/Elecs_LQpure_m%i.png",plot_dir,int(m_LQ));
    sprintf(el_2, "%s/Elecs_LQint_m%i.png",plot_dir,int(m_LQ));

    if(draw_muons){
    //print lq pure weights for all years
    h17_mumu_LQpure_u->SetLineColor(kBlue);
    h18_mumu_LQpure_u->SetLineColor(kGreen);

    TCanvas *c_mu_pwt = new TCanvas("c_mu_pwt", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQpure_u->SetTitle("Muons, m_LQ = 1.0 TeV");
    h18_mumu_LQpure_u->Draw("hist");
    h17_mumu_LQpure_u->Draw("hist same");
    h16_mumu_LQpure_u->Draw("hist same");
            

    TLegend *leg0 = new TLegend(x_start, y_start, x_end, y_end);
    leg0->AddEntry(h16_mumu_LQpure_u, "u-LQpure16", "l");
    leg0->AddEntry(h17_mumu_LQpure_u, "u-LQpure17", "l");
    leg0->AddEntry(h18_mumu_LQpure_u, "u-LQpure18", "l");
    leg0->Draw();
            
    //c_mu_pwt->Print("Misc_plots/template_plots/Muons_LQpure_u.png");
    delete c_mu_pwt;
    //print lq int weights for all years
    h17_mumu_LQint_u->SetLineColor(kRed);
    h18_mumu_LQint_u->SetLineColor(kGreen);

    TCanvas *c_mu_iwt = new TCanvas("c_mu_iwt", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQint_u->SetTitle("Muons: u-LQint");
    h18_mumu_LQint_u->Draw("hist");
    h17_mumu_LQint_u->Draw("hist same");
    h16_mumu_LQint_u->Draw("hist same");
            

    TLegend *leg = new TLegend(x_start, y_start, x_end, y_end);
    leg->AddEntry(h16_mumu_LQint_u, "u-LQint16", "l");
    leg->AddEntry(h17_mumu_LQint_u, "u-LQint17", "l");
    leg->AddEntry(h18_mumu_LQint_u, "u-LQint18", "l");
    leg->Draw();
            
    //c_mu_iwt->Print("Misc_plots/template_plots/Muons_LQint_u.png");
    delete c_mu_iwt;

      //print lq pure for all years
    h17_mumu_LQpure_d->SetLineColor(kBlue);
    h18_mumu_LQpure_d->SetLineColor(kRed);

    TCanvas *c_mu_lqp = new TCanvas("c_mu_lqp", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQpure_d->SetTitle("Muons: d-LQpure");
    h18_mumu_LQpure_d->Draw("hist");
    h17_mumu_LQpure_d->Draw("hist same");
    h16_mumu_LQpure_d->Draw("hist same");
            

    TLegend *leg_2 = new TLegend(x_start, y_start, x_end, y_end);
    leg_2->AddEntry(h16_mumu_LQpure_d, "LQpure16", "l");
    leg_2->AddEntry(h17_mumu_LQpure_d, "LQpure17", "l");
    leg_2->AddEntry(h18_mumu_LQpure_d, "LQpure18", "l");
    leg_2->Draw();
            
    //c_mu_lqp->Print("Misc_plots/template_plots/Muons_LQpure_d.png");
    delete c_mu_lqp;
   
    //print lq int for all years
    h17_mumu_LQint_d->SetLineColor(kBlue);
    h18_mumu_LQint_d->SetLineColor(kRed);

    TCanvas *c_mu_lqi = new TCanvas("c_mu_lqi", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQint_d->SetTitle("Muons: d-LQint");
    h18_mumu_LQint_d->Draw("hist");
    h17_mumu_LQint_d->Draw("hist same");
    h16_mumu_LQint_d->Draw("hist same");
            

    TLegend *leg_1 = new TLegend(x_start, y_start, x_end, y_end);
    leg_1->AddEntry(h16_mumu_LQint_d, "d-LQint16", "l");
    leg_1->AddEntry(h17_mumu_LQint_d, "d-LQint17", "l");
    leg_1->AddEntry(h18_mumu_LQint_d, "d-LQint18", "l");
    leg_1->Draw();
            
    //c_mu_lqi->Print("Misc_plots/template_plots/Muons_LQint_d.png");
    delete c_mu_lqi;
    }

    if(draw_electrons){
    //=============================================================================================
    //print lq pure weights for all years
    h17_elel_LQpure_u->SetLineColor(kBlue);
    h18_elel_LQpure_u->SetLineColor(kGreen);

    TCanvas *c_el_pwt = new TCanvas("c_el_pwt", "Histograms", 200, 10, 900, 700);
    h18_elel_LQpure_u->SetTitle("Electrons, m_LQ = 1.0 TeV");
    h18_elel_LQpure_u->Draw("hist");
    h17_elel_LQpure_u->Draw("hist same");
    h16_elel_LQpure_u->Draw("hist same");
    
    TLegend *leg0 = new TLegend(x_start, y_start, x_end, y_end);
    leg0->AddEntry(h16_elel_LQpure_u, "u-LQpure16", "l");
    leg0->AddEntry(h17_elel_LQpure_u, "u-LQpure17", "l");
    leg0->AddEntry(h18_elel_LQpure_u, "u-LQpure18", "l");
    leg0->Draw();
            
   // c_el_pwt->Print("Misc_plots/template_plots/Elecs_LQpure_u.png");
    delete c_el_pwt;
    //print lq int weights for all years
    h17_elel_LQint_u->SetLineColor(kRed);
    h18_elel_LQint_u->SetLineColor(kGreen);

    TCanvas *c_el_iwt = new TCanvas("c_el_iwt", "Histograms", 200, 10, 900, 700);
    h18_elel_LQint_u->SetTitle("Elecs: u-LQint");
    h18_elel_LQint_u->Draw("hist");
    h17_elel_LQint_u->Draw("hist same");
    h16_elel_LQint_u->Draw("hist same");

    TLegend *leg = new TLegend(x_start, y_start, x_end, y_end);
    leg->AddEntry(h16_elel_LQint_u, "u-LQint16", "l");
    leg->AddEntry(h17_elel_LQint_u, "u-LQint17", "l");
    leg->AddEntry(h18_elel_LQint_u, "u-LQint18", "l");
    leg->Draw();
            
    //c_el_iwt->Print("Misc_plots/template_plots/Elecs_LQint_u.png");
    delete c_el_iwt;

      //print lq pure for all years
    h17_elel_LQpure_d->SetLineColor(kBlue);
    h18_elel_LQpure_d->SetLineColor(kGreen);

    TCanvas *c_el_lqp = new TCanvas("c_el_lqp", "Histograms", 200, 10, 900, 700);
    h18_elel_LQpure_d->SetTitle("Elecs: d-LQpure");
    h18_elel_LQpure_d->Draw("hist");
    h17_elel_LQpure_d->Draw("hist same");
    h16_elel_LQpure_d->Draw("hist same");

    TLegend *leg_2 = new TLegend(x_start, y_start, x_end, y_end);
    leg_2->AddEntry(h16_elel_LQpure_d, "LQpure16", "l");
    leg_2->AddEntry(h17_elel_LQpure_d, "LQpure17", "l");
    leg_2->AddEntry(h18_elel_LQpure_d, "LQpure18", "l");
    leg_2->Draw();
            
    //c_el_lqp->Print("Misc_plots/template_plots/Elecs_LQpure_d.png");
    delete c_el_lqp;
   
    //print lq int for all years
    h17_elel_LQint_d->SetLineColor(kBlue);
    h18_elel_LQint_d->SetLineColor(kGreen);

    TCanvas *c_el_lqi = new TCanvas("c_el_lqi", "Histograms", 200, 10, 900, 700);
    h18_elel_LQint_d->SetTitle("Elecs: d-LQint");
    h18_elel_LQint_d->Draw("hist");
    h17_elel_LQint_d->Draw("hist same");
    h16_elel_LQint_d->Draw("hist same");

    TLegend *leg_1 = new TLegend(x_start, y_start, x_end, y_end);
    leg_1->AddEntry(h16_elel_LQint_d, "d-LQint16", "l");
    leg_1->AddEntry(h17_elel_LQint_d, "d-LQint17", "l");
    leg_1->AddEntry(h18_elel_LQint_d, "d-LQint18", "l");
    leg_1->Draw();
            
    //c_el_lqi->Print("Misc_plots/template_plots/Elecs_LQint_d.png");
    delete c_el_lqi;
    }
 /*   
//print el lqp for all masses in 2016
    if(draw_electrons){
    char leg_entry[50];
    sprintf(leg_entry,"m = %i TeV",int(m_LQ/1000));
    c16_el_lqpall->cd();
    if(i==i_start)
    {
        h16_elel_LQpure_u->SetLineColor(kRed);
        h16_elel_LQpure_u->SetTitle(el_title_all_1);
        h16_elel_LQpure_u->Draw("hist");
        leg_lqpall->AddEntry(h16_elel_LQpure_u,leg_entry);
    }
    else 
    {
        h16_elel_LQpure_u->SetLineColor(kGreen+((i-2)*11));
        h16_elel_LQpure_u->Draw("hist same");
        leg_lqpall->AddEntry(h16_elel_LQpure_u,leg_entry);
    }
//print el lqp for all masses in 2016

    c16_el_lqiall->cd();
    if(i==i_start)
    {
        h16_elel_LQint_u->SetLineColor(kBlue);
        h16_elel_LQint_u->SetTitle(el_title_all_2);
        h16_elel_LQint_u->Draw("hist");
        leg_lqiall->AddEntry(h16_elel_LQint_u,leg_entry);
    }
    else 
    {
        h16_elel_LQint_u->SetLineColor(kBlue+((i)*10));
        h16_elel_LQint_u->Draw("hist same");
        leg_lqiall->AddEntry(h16_elel_LQint_u,leg_entry);
    }
}
//print mu lqp for all masses in 2016
    if(draw_muons){
     char leg_entry[50];
    sprintf(leg_entry,"m = %i TeV",int(m_LQ/1000));
    c16_mu_lqpall->cd();
    if(i==i_start)
    {
        h16_mumu_LQpure_u->SetLineColor(kRed);
        h16_mumu_LQpure_u->SetTitle(mu_title_all_1);
        h16_mumu_LQpure_u->Draw("hist");
        leg_lqpall->AddEntry(h16_mumu_LQpure_u,leg_entry);
    }
    else 
    {
        h16_mumu_LQpure_u->SetLineColor(kGreen+((i-2)*11));
        h16_mumu_LQpure_u->Draw("hist same");
       leg_lqpall->AddEntry(h16_mumu_LQpure_u,leg_entry);
    }
//print mu lqp for all masses in 2016

    c16_mu_lqiall->cd();
    if(i==i_start)
    {
        h16_mumu_LQint_u->SetLineColor(kBlue);
        h16_mumu_LQint_u->SetTitle(mu_title_all_2);
        h16_mumu_LQint_u->Draw("hist");
       // leg_lqiall->AddEntry(h16_mumu_LQint_u,"m = %i"+int(m_LQ));
    }
    else 
    {
        h16_mumu_LQint_u->SetLineColor(kBlue+((i)*10));
        h16_mumu_LQint_u->Draw("hist same");
      //  leg_lqiall->AddEntry(h16_mumu_LQint_u,"m = %i"+int(m_LQ));
    }
}
}    
    

    if(draw_electrons){
    c16_el_lqpall->cd();
    leg_lqpall->Draw();
            
    c16_el_lqpall->Print(el_all_lqp16);
    delete c16_el_lqpall;

    c16_el_lqiall->cd();
    leg_lqiall->Draw();
            
   // c16_el_lqiall->Print(el_all_lqi16);
    delete c16_el_lqiall;
}
    if(draw_muons){
    c16_mu_lqpall->cd();
    leg_lqpall->Draw();
            
    c16_mu_lqpall->Print(mu_all_lqp16);
    delete c16_mu_lqpall;

    c16_mu_lqiall->cd();
    leg_lqiall->Draw();
            
   // c16_mu_lqiall->Print(mu_all_lqi16);
    delete c16_mu_lqiall;
  }
  */
  }





