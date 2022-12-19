#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"
#include <iostream>

TLine *line1, *line2, *line3, *line4, *line5, *line6, *line7, *line8;

void set_line_attributes(float ymin, float ymax){
	
	    line1 = new TLine(8,ymin,8,ymax);
        line2 = new TLine(14,ymin,14,ymax);
        line3 = new TLine(20,ymin,20,ymax);
        line4 = new TLine(28,ymin,28,ymax);
        line5 = new TLine(34,ymin,34,ymax);
        line6 = new TLine(40,ymin,40,ymax);
        line7 = new TLine(48,ymin,48,ymax);
        line8 = new TLine(54,ymin,54,ymax);	

	    line1->SetLineColor(kGreen);
        line2->SetLineColor(kGreen);
        line3->SetLineColor(kOrange);
        line4->SetLineColor(kGreen);
        line5->SetLineColor(kGreen);
        line6->SetLineColor(kOrange);
        line7->SetLineColor(kGreen);
        line8->SetLineColor(kGreen);

        line1->SetLineStyle(9);
        line2->SetLineStyle(9);
        line3->SetLineStyle(9);
        line4->SetLineStyle(9);
        line5->SetLineStyle(9);
        line6->SetLineStyle(9);
        line7->SetLineStyle(9);
        line8->SetLineStyle(9);

        line1->SetLineWidth(2);
        line2->SetLineWidth(2);
        line3->SetLineWidth(2);
        line4->SetLineWidth(2);
        line5->SetLineWidth(2);
        line6->SetLineWidth(2);
        line7->SetLineWidth(2);
        line8->SetLineWidth(2);

        line1->Draw("same");
        line2->Draw("same");
        line3->Draw("same");
        line4->Draw("same");
        line5->Draw("same");
        line6->Draw("same");
        line7->Draw("same");
        line8->Draw("same");
}
	

void LQ_draw_templates(){

    
        bool ss = false;
        bool use_xF =false;
        //bool use_LQ_denom=true;
        bool draw_muons = true;
        bool draw_electrons = true;
        const string sys_label = "";
        
        //char *plot_dir = "Paper_plots/template_plots";
        char *plot_dir = "AN_plots/LQ_templates/";
        int i_start = 1; 
        int i_end = 4;
        float x_start = 0.1;
        float x_end = 0.3;
        float y_start = 0.75;
        float y_end = 0.9;
       
       // std::cout << "enter m_LQ:";        
       // std::cin >> m_LQ; 
        gStyle->SetOptStat(0);
        gROOT->SetBatch(1);

        float ymin = 0, ymax = 0;

    
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
        Double_t m_LQ=2000.;
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
            char mu_fname1[100], mu_fname2[100], mu_fname3[100], mu_fname4[100], mu_fname5[100], mu_fname6[100], el_fname1[100], el_fname2[100], el_fname3[100], el_fname4[100], el_fname5[100], el_fname6[100];

            

            if(draw_muons){

           //scalar u
            sprintf(title, "mumu%i_LQpure_u%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_u = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_u->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_u%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_u = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_u->SetDirectory(0);
            //scalar d
            sprintf(title, "mumu%i_LQpure_d%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_d = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_d->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_d%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_d = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_d->SetDirectory(0);
            //vector u
             sprintf(title, "mumu%i_LQpure_u_vec%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_u_vec = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_u_vec->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_u_vec%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_u_vec = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_u_vec->SetDirectory(0);
            //vector d
            sprintf(title, "mumu%i_LQpure_d_vec%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure_d_vec = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure_d_vec->SetDirectory(0);
            sprintf(title, "mumu%i_LQint_d_vec%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint_d_vec = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_mumu_LQint_d_vec->SetDirectory(0);

            sprintf(mu_fname1, "%s/%i/Mu%i_Scm_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(mu_fname2, "%s/%i/Mu%i_Ssm_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(mu_fname3, "%s/%i/Mu%i_Vcm_pure_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(mu_fname4, "%s/%i/Mu%i_Vcm_int_m%i.png", plot_dir,year, year%2000, int(m_LQ));
            sprintf(mu_fname5, "%s/%i/Mu%i_Vsm_pure_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(mu_fname6, "%s/%i/Mu%i_Vsm_int_m%i.png", plot_dir, year,year%2000, int(m_LQ));

            bool make_ud = false;

            //gen_mc_SM_template(t_elel_mc,  h_elel_sym, h_elel_asym, h_elel_alpha, year, FLAG_ELECTRONS, use_xF, sys_label );
            gen_mc_LQ_template(t_mumu_mc,  h_mumu_LQpure_u, h_mumu_LQint_u, h_mumu_LQpure_d, h_mumu_LQint_d, h_mumu_LQpure_u_vec, h_mumu_LQint_u_vec, h_mumu_LQpure_d_vec, h_mumu_LQint_d_vec, year, m_LQ, FLAG_MUONS, make_ud, false,  use_xF, sys_label );

            
           auto h1_mumu_LQpure_u = convert3d(h_mumu_LQpure_u);
           auto h1_mumu_LQint_u = convert3d(h_mumu_LQint_u);
           auto h1_mumu_LQpure_d = convert3d(h_mumu_LQpure_d);
           auto h1_mumu_LQint_d = convert3d(h_mumu_LQint_d);
           auto h1_mumu_LQpure_u_vec = convert3d(h_mumu_LQpure_u_vec);
           auto h1_mumu_LQint_u_vec = convert3d(h_mumu_LQint_u_vec);
           auto h1_mumu_LQpure_d_vec = convert3d(h_mumu_LQpure_d_vec);
           auto h1_mumu_LQint_d_vec = convert3d(h_mumu_LQint_d_vec);


            h1_mumu_LQpure_u->SetLineColor(kBlue);
            h1_mumu_LQint_u->SetLineColor(kRed);
            h1_mumu_LQpure_d->SetLineColor(kBlue);
            h1_mumu_LQint_d->SetLineColor(kRed );
            h1_mumu_LQpure_u_vec->SetLineColor(kBlue);
            h1_mumu_LQint_u_vec->SetLineColor(kRed);
            h1_mumu_LQpure_d_vec->SetLineColor(kBlue);
            h1_mumu_LQint_d_vec->SetLineColor(kRed );


            
            h1_mumu_LQpure_u->SetLineWidth(2);
            h1_mumu_LQint_u->SetLineWidth(2);
            h1_mumu_LQpure_d->SetLineWidth(2);
            h1_mumu_LQint_d->SetLineWidth(2);
            h1_mumu_LQpure_u_vec->SetLineWidth(2);
            h1_mumu_LQint_u_vec->SetLineWidth(2);
            h1_mumu_LQpure_d_vec->SetLineWidth(2);
            h1_mumu_LQint_d_vec->SetLineWidth(2);
         
            sprintf(mu_title, "Channel : Muons, %.1f TeV S_{#mu c}",m_LQ/1000,year);
            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure_u->SetTitle(mu_title); 
           h1_mumu_LQpure_u->Draw("hist");
           h1_mumu_LQint_u->Draw("hist same");
          // h1_mumu_pl->Draw("hist");
          //  h1_mumu_alpha->Draw("hist same");
            //h1_mumu_mn->Draw("hist same");
            ymax = h1_mumu_LQpure_u->GetMaximum();
            ymin = h1_mumu_LQpure_u->GetMinimum();
	    
	    
            set_line_attributes(ymin, ymax);
            
            

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_mumu_LQpure_u, "Pure LQ Template", "l");
            leg1->AddEntry(h1_mumu_LQint_u, "Intereference LQ Template", "l");
           //leg1->AddEntry(h1_mumu_pl, "Plus Template", "l");
           //leg1->AddEntry(h1_mumu_mn, "Minus Template", "l");
            //leg1->AddEntry(h1_mumu_alpha, "alpha Template", "l");
            leg1->Draw();
            c_mumu1->Print(mu_fname1);
            delete c_mumu1;
            
            sprintf(mu_title, "Channel : Muons, %.1f TeV S_{#mu s}",m_LQ/1000,year);
            TCanvas *c_mumu2 = new TCanvas("c_mumu2", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure_d->SetTitle(mu_title);
            h1_mumu_LQpure_d->Draw("hist");
            h1_mumu_LQint_d->Draw("hist same");

            ymax = h1_mumu_LQpure_d->GetMaximum();
            ymin = h1_mumu_LQpure_d->GetMinimum();

	    set_line_attributes(ymin, ymax);
            

            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_mumu_LQpure_d,"Pure LQ Template","l");
            leg2->AddEntry(h1_mumu_LQint_d,"Interference LQ Template","l");
            leg2->Draw();

            c_mumu2->Print(mu_fname2);
            delete c_mumu2;

            sprintf(mu_title, "Channel : Muons, %.1f TeV V_{#mu c}",m_LQ/1000,year);
            TCanvas *c_mumu3 = new TCanvas("c_mumu3", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure_u_vec->SetTitle(mu_title);
            h1_mumu_LQpure_u_vec->Draw("hist");
            //h1_mumu_LQint_u->Draw("hist same");

            ymax = h1_mumu_LQpure_u_vec->GetMaximum();
            ymin = h1_mumu_LQpure_u_vec->GetMinimum();

	    set_line_attributes(ymin, ymax);
            

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h1_mumu_LQpure_u_vec,"Pure LQ Template","l");
            //leg3->AddEntry(h1_mumu_LQint_d,"d-LQint Template","l");
            leg3->Draw();
            
            c_mumu3->Print(mu_fname3);
            delete c_mumu3;

            sprintf(mu_title, "Channel : Muons, %.1f TeV V_{#mu c}",m_LQ/1000,year);
            TCanvas *c_mumu4 = new TCanvas("c_mumu4", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQint_u_vec->SetTitle(mu_title);
            h1_mumu_LQint_u_vec->Draw("hist");
            //h1_mumu_LQint_u->Draw("hist same");

            ymax = h1_mumu_LQint_u_vec->GetMaximum();
            ymin = h1_mumu_LQint_u_vec->GetMinimum();

	    set_line_attributes(ymin, ymax);
            

            TLegend *leg4 = new TLegend(x_start, y_start, x_end, y_end);
            leg4->AddEntry(h1_mumu_LQint_u_vec,"Interference LQ Template","l");
            //leg3->AddEntry(h1_mumu_LQint_d,"d-LQint Template","l");
            leg4->Draw();
            
            c_mumu4->Print(mu_fname4);
            delete c_mumu4;

            sprintf(mu_title, "Channel : Muons, %.1f TeV V_{#mu s}",m_LQ/1000,year);
            TCanvas *c_mumu5 = new TCanvas("c_mumu5", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure_d_vec->SetTitle(mu_title);
            h1_mumu_LQpure_d_vec->Draw("hist");
            //h1_mumu_LQint_u->Draw("hist same");

            ymax = h1_mumu_LQpure_d_vec->GetMaximum();
            ymin = h1_mumu_LQpure_d_vec->GetMinimum();

	    set_line_attributes(ymin, ymax);
            

            TLegend *leg5 = new TLegend(x_start, y_start, x_end, y_end);
            leg5->AddEntry(h1_mumu_LQpure_d_vec,"Pure LQ Template","l");
            //leg3->AddEntry(h1_mumu_LQint_d,"d-LQint Template","l");
            leg5->Draw();
            
            c_mumu5->Print(mu_fname5);
            delete c_mumu5;

            sprintf(mu_title, "Channel : Muons, %.1f TeV V_{#mu s}",m_LQ/1000,year);
            TCanvas *c_mumu6 = new TCanvas("c_mumu6", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQint_d_vec->SetTitle(mu_title);
            h1_mumu_LQint_d_vec->Draw("hist");
            //h1_mumu_LQint_u->Draw("hist same");

            ymax = h1_mumu_LQint_d_vec->GetMaximum();
            ymin = h1_mumu_LQint_d_vec->GetMinimum();
	    set_line_attributes(ymin, ymax);

            

            TLegend *leg6 = new TLegend(x_start, y_start, x_end, y_end);
            leg6->AddEntry(h1_mumu_LQint_d_vec,"Interference LQ Template","l");
            //leg3->AddEntry(h1_mumu_LQint_d,"d-LQint Template","l");
            leg6->Draw();
            
            c_mumu6->Print(mu_fname6);
            delete c_mumu6;

        }

        if(draw_electrons){

           //scalar u
            sprintf(title, "elel%i_LQpure_u%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_u = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_u->SetDirectory(0);
            sprintf(title, "elel%i_LQint_u%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_u = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_u->SetDirectory(0);
            //scalar d
            sprintf(title, "elel%i_LQpure_d%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_d = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_d->SetDirectory(0);
            sprintf(title, "elel%i_LQint_d%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_d = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_d->SetDirectory(0);
            //vector u
             sprintf(title, "elel%i_LQpure_u_vec%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_u_vec = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_u_vec->SetDirectory(0);
            sprintf(title, "elel%i_LQint_u_vec%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_u_vec = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_u_vec->SetDirectory(0);
            //vector d
            sprintf(title, "elel%i_LQpure_d_vec%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure_d_vec = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQpure_d_vec->SetDirectory(0);
            sprintf(title, "elel%i_LQint_d_vec%s", year %2000, sys_label.c_str());
            auto h_elel_LQint_d_vec = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
            h_elel_LQint_d_vec->SetDirectory(0);

            sprintf(el_fname1, "%s/%i/El%i_Sce_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(el_fname2, "%s/%i/El%i_Sse_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(el_fname3, "%s/%i/El%i_Vce_pure_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(el_fname4, "%s/%i/El%i_Vce_int_m%i.png", plot_dir,year, year%2000, int(m_LQ));
            sprintf(el_fname5, "%s/%i/El%i_Vse_pure_m%i.png", plot_dir, year,year%2000, int(m_LQ));
            sprintf(el_fname6, "%s/%i/El%i_Vse_int_m%i.png", plot_dir, year,year%2000, int(m_LQ));

            bool make_ud = false;

            //gen_mc_SM_template(t_elel_mc,  h_elel_sym, h_elel_asym, h_elel_alpha, year, FLAG_ELECTRONS, use_xF, sys_label );
            gen_mc_LQ_template(t_elel_mc,  h_elel_LQpure_u, h_elel_LQint_u, h_elel_LQpure_d, h_elel_LQint_d, h_elel_LQpure_u_vec, h_elel_LQint_u_vec, h_elel_LQpure_d_vec, h_elel_LQint_d_vec, year, m_LQ, FLAG_ELECTRONS, make_ud,false, use_xF, sys_label );

            
           auto h1_elel_LQpure_u = convert3d(h_elel_LQpure_u);
           auto h1_elel_LQint_u = convert3d(h_elel_LQint_u);
           auto h1_elel_LQpure_d = convert3d(h_elel_LQpure_d);
           auto h1_elel_LQint_d = convert3d(h_elel_LQint_d);
           auto h1_elel_LQpure_u_vec = convert3d(h_elel_LQpure_u_vec);
           auto h1_elel_LQint_u_vec = convert3d(h_elel_LQint_u_vec);
           auto h1_elel_LQpure_d_vec = convert3d(h_elel_LQpure_d_vec);
           auto h1_elel_LQint_d_vec = convert3d(h_elel_LQint_d_vec);


            h1_elel_LQpure_u->SetLineColor(kBlue);
            h1_elel_LQint_u->SetLineColor(kRed);
            h1_elel_LQpure_d->SetLineColor(kBlue);
            h1_elel_LQint_d->SetLineColor(kRed );
            h1_elel_LQpure_u_vec->SetLineColor(kBlue);
            h1_elel_LQint_u_vec->SetLineColor(kRed);
            h1_elel_LQpure_d_vec->SetLineColor(kBlue);
            h1_elel_LQint_d_vec->SetLineColor(kRed );

	    
            
            h1_elel_LQpure_u->SetLineWidth(2);
            h1_elel_LQint_u->SetLineWidth(2);
            h1_elel_LQpure_d->SetLineWidth(2);
            h1_elel_LQint_d->SetLineWidth(2);
            h1_elel_LQpure_u_vec->SetLineWidth(2);
            h1_elel_LQint_u_vec->SetLineWidth(2);
            h1_elel_LQpure_d_vec->SetLineWidth(2);
            h1_elel_LQint_d_vec->SetLineWidth(2);
         
            sprintf(el_title, "Channel : Electrons, %.1f TeV S_{ec}",m_LQ/1000,year);
            TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure_u->SetTitle(el_title); 
           h1_elel_LQpure_u->Draw("hist");
           h1_elel_LQint_u->Draw("hist same");
          // h1_elel_pl->Draw("hist");
          //  h1_elel_alpha->Draw("hist same");
            //h1_elel_mn->Draw("hist same");

           ymax = h1_elel_LQpure_u->GetMaximum();
            ymin = h1_elel_LQpure_u->GetMinimum();
	    set_line_attributes(ymin, ymax);

           
            

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_elel_LQpure_u, "Pure LQ Template", "l");
            leg1->AddEntry(h1_elel_LQint_u, "Intereference LQ Template", "l");
           //leg1->AddEntry(h1_elel_pl, "Plus Template", "l");
           //leg1->AddEntry(h1_elel_mn, "Minus Template", "l");
            //leg1->AddEntry(h1_elel_alpha, "alpha Template", "l");
            leg1->Draw();
            c_elel1->Print(el_fname1);
            delete c_elel1;
            
            sprintf(el_title, "Channel : Electrons, %.1f TeV S_{es}",m_LQ/1000,year);
            TCanvas *c_elel2 = new TCanvas("c_elel2", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure_d->SetTitle(el_title);
            h1_elel_LQpure_d->Draw("hist");
            h1_elel_LQint_d->Draw("hist same");

            ymax = h1_elel_LQpure_d->GetMaximum();
            ymin = h1_elel_LQpure_d->GetMinimum();
	    set_line_attributes(ymin, ymax);

            

            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_elel_LQpure_d,"Pure LQ Template","l");
            leg2->AddEntry(h1_elel_LQint_d,"Interference LQ Template","l");
            leg2->Draw();

            c_elel2->Print(el_fname2);
            delete c_elel2;

            sprintf(el_title, "Channel : Electrons, %.1f TeV V_{ec}",m_LQ/1000,year);
            TCanvas *c_elel3 = new TCanvas("c_elel3", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure_u_vec->SetTitle(el_title);
            h1_elel_LQpure_u_vec->Draw("hist");
            //h1_elel_LQint_u->Draw("hist same");

            ymax = h1_elel_LQpure_u_vec->GetMaximum();
            ymin = h1_elel_LQpure_u_vec->GetMinimum();
	    set_line_attributes(ymin, ymax);

            

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h1_elel_LQpure_u_vec,"Pure LQ Template","l");
            //leg3->AddEntry(h1_elel_LQint_d,"d-LQint Template","l");
            leg3->Draw();
            
            c_elel3->Print(el_fname3);
            delete c_elel3;

            sprintf(el_title, "Channel : Electrons, %.1f TeV V_{ec}",m_LQ/1000,year);
            TCanvas *c_elel4 = new TCanvas("c_elel4", "Histograms", 200, 10, 900, 700);
            h1_elel_LQint_u_vec->SetTitle(el_title);
            h1_elel_LQint_u_vec->Draw("hist");
            //h1_elel_LQint_u->Draw("hist same");

            ymax = h1_elel_LQint_u_vec->GetMaximum();
            ymin = h1_elel_LQint_u_vec->GetMinimum();
	    set_line_attributes(ymin, ymax);

            

            TLegend *leg4 = new TLegend(x_start, y_start, x_end, y_end);
            leg4->AddEntry(h1_elel_LQint_u_vec,"Interference LQ Template","l");
            //leg3->AddEntry(h1_elel_LQint_d,"d-LQint Template","l");
            leg4->Draw();
            
            c_elel4->Print(el_fname4);
            delete c_elel4;

            sprintf(el_title, "Channel : Electrons, %.1f TeV V_{es}",m_LQ/1000,year);
            TCanvas *c_elel5 = new TCanvas("c_elel5", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure_d_vec->SetTitle(el_title);
            h1_elel_LQpure_d_vec->Draw("hist");
            //h1_elel_LQint_u->Draw("hist same");

            ymax = h1_elel_LQpure_d_vec->GetMaximum();
            ymin = h1_elel_LQpure_d_vec->GetMinimum();

	    set_line_attributes(ymin, ymax);
            

            TLegend *leg5 = new TLegend(x_start, y_start, x_end, y_end);
            leg5->AddEntry(h1_elel_LQpure_d_vec,"Pure LQ Template","l");
            //leg3->AddEntry(h1_elel_LQint_d,"d-LQint Template","l");
            leg5->Draw();
            
            c_elel5->Print(el_fname5);
            delete c_elel5;

            sprintf(el_title, "Channel : Electrons, %.1f TeV V_{es}",m_LQ/1000,year);
            TCanvas *c_elel6 = new TCanvas("c_elel6", "Histograms", 200, 10, 900, 700);
            h1_elel_LQint_d_vec->SetTitle(el_title);
            h1_elel_LQint_d_vec->Draw("hist");
            //h1_elel_LQint_u->Draw("hist same");

            ymax = h1_elel_LQint_d_vec->GetMaximum();
            ymin = h1_elel_LQint_d_vec->GetMinimum();
	    set_line_attributes(ymin, ymax);	    
	    

            TLegend *leg6 = new TLegend(x_start, y_start, x_end, y_end);
            leg6->AddEntry(h1_elel_LQint_d_vec,"Interference LQ Template","l");
            //leg3->AddEntry(h1_elel_LQint_d,"d-LQint Template","l");
            leg6->Draw();
            
            c_elel6->Print(el_fname6);
            delete c_elel6;

        }

    }
    /*
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





