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
        bool draw_electrons = true;
       const string sys_label = "";
       const int n_rap_bins = 4;
        float rap_bins[] = {0., 0.6, 1., 1.5,  2.4};
        //char *plot_dir = "Paper_plots/template_plots";
        char *plot_dir = "Misc_plots/template_plots2";

       for(int i=1;i<=10;i++){
        Double_t m_LQ=1000.*i;
        printf("=========================\n m_LQ = %f, draw_muons = %d, draw_electrons = %d \n=========================\n",m_LQ,draw_muons,draw_electrons );

        TH1F *h16_mumu_pl, *h16_mumu_mn, *h16_mumu_LQpure, *h16_mumu_LQint; TH1D *h16_mumu_LQpure_wt, *h16_mumu_LQint_wt;
        TH1F *h17_mumu_pl, *h17_mumu_mn, *h17_mumu_LQpure, *h17_mumu_LQint; TH1D *h17_mumu_LQpure_wt, *h17_mumu_LQint_wt;
        TH1F *h18_mumu_pl, *h18_mumu_mn, *h18_mumu_LQpure, *h18_mumu_LQint; TH1D *h18_mumu_LQpure_wt, *h18_mumu_LQint_wt;

        TH1F *h16_elel_pl, *h16_elel_mn, *h16_elel_LQpure, *h16_elel_LQint; TH1D *h16_elel_LQpure_wt, *h16_elel_LQint_wt;
        TH1F *h17_elel_pl, *h17_elel_mn, *h17_elel_LQpure, *h17_elel_LQint; TH1D *h17_elel_LQpure_wt, *h17_elel_LQint_wt;
        TH1F *h18_elel_pl, *h18_elel_mn, *h18_elel_LQpure, *h18_elel_LQint; TH1D *h18_elel_LQpure_wt, *h18_elel_LQint_wt;

        float x_start = 0.75;
        float x_end = 0.9;

        float y_start = 0.75;
        float y_end = 0.9;

       // std::cout << "enter m_LQ:";        
       // std::cin >> m_LQ; 
        gStyle->SetOptStat(0);
        gROOT->SetBatch(1);
        //int year = 2017;
        for (int year = 2016; year<=2018; year++)
        {

        init(year);
        setup_all_SFs(year);
            
            char title[100];

            char mu_title[100], el_title[100];
            char mu_fname1[100], mu_fname2[100], mu_fname3[100], el_fname1[100], el_fname2[100], el_fname3[100];

            sprintf(mu_title, "Muons, m_LQ=%0.1f, year=%i",m_LQ,year);
            sprintf(el_title, "Electrons, m_LQ=%0.1f, year=%i",m_LQ,year);

            if(draw_muons){

            auto h_mumu_sym = new TH3F(title, "Symmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_sym->SetDirectory(0);
            sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
            auto h_mumu_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_alpha->SetDirectory(0);
            sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
            auto h_mumu_asym = new TH3F(title, "Asymmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_asym->SetDirectory(0);
            sprintf(title, "mumu%i_LQpure%s", year %2000, sys_label.c_str());
            auto h_mumu_LQpure = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_LQpure->SetDirectory(0);
            sprintf(title, "mumu%i_LQint%s", year %2000, sys_label.c_str());
            auto h_mumu_LQint = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_LQint->SetDirectory(0);
            auto h_mu_LQpure_wt = new TH1D(title, "Weights distribution of LQpure", 100, 0., .2);
            h_mu_LQpure_wt->SetDirectory(0);
            auto h_mu_LQint_wt = new TH1D(title, "Weights distribution of LQint", 100, 0., .2);
            h_mu_LQint_wt->SetDirectory(0);
            

            gen_mc_template(t_mumu_mc, h_mumu_sym, h_mumu_asym, h_mumu_alpha,h_mumu_LQpure, h_mumu_LQint, 
                year, m_LQ, FLAG_MUONS, use_xF,"");

            sprintf(mu_fname1, "%s/Mu%i_MC_SM_m%i.png", plot_dir, year%2000,int(m_LQ));
            //sprintf(mu_fname2, "%s/MuMu%i_M_fit_temps.png", plot_dir, year);
            sprintf(mu_fname2, "%s/Mu%i_MC_LQ_m%i.png", plot_dir, year%2000, int(m_LQ));
            sprintf(mu_fname3, "%s/Mu%i_MC_LQwts_m%i.png", plot_dir, year%2000, int(m_LQ));

            auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
            auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
            /*
            h_mumu_pl.Scale(0.5);
            h_mumu_mn.Scale(0.5);
            double norm = 3./4./(2.+alpha_denom);
            h_mumu_alpha->Scale(norm);
            h_mumu_LQpure->Scale(norm);
            h_mumu_LQint->Scale(norm);
           */
           // h_mumu_pl.Print("range");
            auto h1_mumu_pl = convert3d(&h_mumu_pl);
           // h1_mumu_pl->Print("range");
            auto h1_mumu_mn = convert3d(&h_mumu_mn);
            auto h1_mumu_alpha = convert3d(h_mumu_alpha);
            auto h1_mumu_sym = convert3d(h_mumu_sym);
            auto h1_mumu_asym = convert3d(h_mumu_asym);
            auto h1_mumu_LQpure = convert3d(h_mumu_LQpure);
            auto h1_mumu_LQint = convert3d(h_mumu_LQint);


            h1_mumu_alpha->SetLineColor(kGreen +3);
            h1_mumu_sym->SetLineColor(kBlue);
            h1_mumu_asym->SetLineColor(kRed+1);
            h1_mumu_pl->SetLineColor(kOrange +7);
            h1_mumu_mn->SetLineColor(kRed + 1);
            h1_mumu_LQpure->SetLineColor(kOrange +9);
            h1_mumu_LQint->SetLineColor(kYellow + 3);
          //  h_mu_LQpure_wt->SetLineColor(kRed);
           // h_mu_LQint_wt->SetLineColor(kGreen+3);

            h1_mumu_alpha->SetLineWidth(2);
            h1_mumu_sym->SetLineWidth(2);
            h1_mumu_asym->SetLineWidth(2);
            h1_mumu_pl->SetLineWidth(2);
            h1_mumu_mn->SetLineWidth(2);
            h1_mumu_LQpure->SetLineWidth(2);
            h1_mumu_LQint->SetLineWidth(2);
           // h_mu_LQpure_wt->SetLineWidth(2);
           // h_mu_LQint_wt->SetLineWidth(2);



            if(year==2016)
            {
                h16_mumu_pl = h1_mumu_pl;
                h16_mumu_mn = h1_mumu_mn;
                h16_mumu_LQpure = h1_mumu_LQpure;
                h16_mumu_LQint = h1_mumu_LQint;
            //    h16_mumu_LQpure_wt = h_mu_LQpure_wt;
             //   h16_mumu_LQint_wt = h_mu_LQint_wt;
            }
            if(year==2017)
            {
                h17_mumu_pl = h1_mumu_pl;
                h17_mumu_mn = h1_mumu_mn;
                h17_mumu_LQpure = h1_mumu_LQpure;
                h17_mumu_LQint = h1_mumu_LQint;
            //    h17_mumu_LQpure_wt = h_mu_LQpure_wt;
             //   h17_mumu_LQint_wt = h_mu_LQint_wt;
            }
            if(year==2018)
            {
                h18_mumu_pl = h1_mumu_pl;
                h18_mumu_mn = h1_mumu_mn;
                h18_mumu_LQpure = h1_mumu_LQpure;
                h18_mumu_LQint = h1_mumu_LQint;
              //  h18_mumu_LQpure_wt = h_mu_LQpure_wt;
               // h18_mumu_LQint_wt = h_mu_LQint_wt;
            }

            h1_mumu_asym->SetMaximum(h1_mumu_sym->GetMaximum()*1.2);
           // h1_mumu_LQpure->SetMaximum(h1_mumu_LQint->GetMaximum()*1.2);
            /*
            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_asym->SetTitle(mu_title); 
            h1_mumu_asym->Draw("hist");
            h1_mumu_sym->Draw("hist same ");

            

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_mumu_asym, "Asym Template", "l");
            leg1->AddEntry(h1_mumu_sym, "Sym Template", "l");
            leg1->Draw();
         //   c_mumu1->Print(mu_fname1);
            delete c_mumu1;
        */
            TCanvas *c_mumu2 = new TCanvas("c_mumu2", "Histograms", 200, 10, 900, 700);
            h1_mumu_LQpure->SetTitle(mu_title);
            h1_mumu_LQpure->Draw("hist");
           // h1_mumu_mn->Draw("hist same");
           // h1_mumu_alpha->Draw("hist same");
            //h1_mumu_LQpure->SetTitle(mu_title);
            //h1_mumu_pl->Draw("hist same ");
            h1_mumu_LQint->Draw("hist same");


           // c_mumu2->cd();
            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
           // leg2->AddEntry(h1_mumu_pl, "Plus Template", "l");
           //leg2->AddEntry(h1_mumu_mn, "Minus Template", "l");
           //leg2->AddEntry(h1_mumu_alpha, "#alpha Template", "l");
            leg2->AddEntry(h1_mumu_LQpure,"LQpure Template","l");
            leg2->AddEntry(h1_mumu_LQint,"LQint Template","l");
            leg2->Draw();

            c_mumu2->Print(mu_fname2);
            delete c_mumu2;
/*
            TCanvas *c_mumu3 = new TCanvas("c_mumu3", "Histograms", 200, 10, 900, 700);
            h_mu_LQint_wt->SetTitle(mu_title);
            h_mu_LQint_wt->Draw("hist");
            h_mu_LQpure_wt->Draw("hist same ");

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h_mu_LQint_wt, "LQint weights*evtwt", "l");
            leg3->AddEntry(h_mu_LQpure_wt, "LQpure weights*evtwt", "l");
            leg3->Draw();
            

         //   c_mumu3->Print(mu_fname3);
            delete c_mumu3;
*/
        }

        if(draw_electrons){

            auto h_elel_sym = new TH3F(title, "Symmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_sym->SetDirectory(0);
            sprintf(title, "elel%i_alpha%s", year %2000, sys_label.c_str());
            auto h_elel_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_alpha->SetDirectory(0);
            sprintf(title, "elel%i_asym%s", year %2000, sys_label.c_str());
            auto h_elel_asym = new TH3F(title, "Asymmetric template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_asym->SetDirectory(0);
            sprintf(title, "elel%i_LQpure%s", year %2000, sys_label.c_str());
            auto h_elel_LQpure = new TH3F(title, "LQpure template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_LQpure->SetDirectory(0);
            sprintf(title, "elel%i_LQint%s", year %2000, sys_label.c_str());
            auto h_elel_LQint = new TH3F(title, "LQint template of mc",
                    n_lq_m_bins, lq_m_bins, n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_LQint->SetDirectory(0);
            auto h_el_LQpure_wt = new TH1D(title, "Weights distribution of LQpure", 100, 0., .2);
            h_el_LQpure_wt->SetDirectory(0);
            auto h_el_LQint_wt = new TH1D(title, "Weights distribution of LQint", 100, 0., .2);
            h_el_LQint_wt->SetDirectory(0);

            gen_mc_template(t_elel_mc, h_elel_sym, h_elel_asym, h_elel_alpha,h_elel_LQpure, h_elel_LQint, 
                year, m_LQ, FLAG_ELECTRONS, use_xF, "");

            sprintf(el_fname1, "%s/El%i_MC_SM_m%i.png", plot_dir, year%2000,int(m_LQ));
            //sprintf(el_fname2, "%s/ElEl%i_M_fit_temps.png", plot_dir, year);
            sprintf(el_fname2, "%s/El%i_MC_LQ_m%i.png", plot_dir, year%2000,int(m_LQ));
            sprintf(el_fname3, "%s/El%i_MC_LQwts_m%i.png", plot_dir, year%2000,int(m_LQ));

            auto h_elel_pl = *h_elel_sym + *h_elel_asym;
            auto h_elel_mn = *h_elel_sym - *h_elel_asym;
            /*
            h_elel_pl.Scale(0.5);
            h_elel_mn.Scale(0.5);
            h_elel_alpha->Scale(norm);
            h_elel_LQpure->Scale(norm);
            h_elel_LQint->Scale(norm);
            
*/
           
            auto h1_elel_pl = convert3d(&h_elel_pl);
            auto h1_elel_mn = convert3d(&h_elel_mn);
            auto h1_elel_alpha = convert3d(h_elel_alpha);
            auto h1_elel_sym = convert3d(h_elel_sym);
            auto h1_elel_asym = convert3d(h_elel_asym);
            auto h1_elel_LQpure = convert3d(h_elel_LQpure);
            auto h1_elel_LQint = convert3d(h_elel_LQint);

            h1_elel_alpha->SetLineColor(kGreen +3);
            h1_elel_sym->SetLineColor(kBlue);
            h1_elel_asym->SetLineColor(kRed+1);
            h1_elel_pl->SetLineColor(kOrange +7);
            h1_elel_mn->SetLineColor(kRed + 1);
            h1_elel_LQpure->SetLineColor(kOrange +9);
            h1_elel_LQint->SetLineColor(kYellow + 3);
           // h_el_LQpure_wt->SetLineColor(kRed);
            //h_el_LQint_wt->SetLineColor(kGreen+3);

            h1_elel_alpha->SetLineWidth(2);
            h1_elel_sym->SetLineWidth(2);
            h1_elel_asym->SetLineWidth(2);
            h1_elel_pl->SetLineWidth(2);
            h1_elel_mn->SetLineWidth(2);
            h1_elel_LQpure->SetLineWidth(2);
            h1_elel_LQint->SetLineWidth(2);
          //  h_el_LQpure_wt->SetLineWidth(2);
           // h_el_LQint_wt->SetLineWidth(2);

            if(year==2016)
            {
                h16_elel_pl = h1_elel_pl;
                h16_elel_mn = h1_elel_mn;
                h16_elel_LQpure = h1_elel_LQpure;
                h16_elel_LQint = h1_elel_LQint;
             //   h16_elel_LQpure_wt = h_el_LQpure_wt;
              //  h16_elel_LQint_wt = h_el_LQint_wt;
            }
            if(year==2017)
            {
                h17_elel_pl = h1_elel_pl;
                h17_elel_mn = h1_elel_mn;
                h17_elel_LQpure = h1_elel_LQpure;
                h17_elel_LQint = h1_elel_LQint;
              //  h17_elel_LQpure_wt = h_el_LQpure_wt;
               // h17_elel_LQint_wt = h_el_LQint_wt;
            }
            if(year==2018)
            {
                h18_elel_pl = h1_elel_pl;
                h18_elel_mn = h1_elel_mn;
                h18_elel_LQpure = h1_elel_LQpure;
                h18_elel_LQint = h1_elel_LQint;
               // h18_elel_LQpure_wt = h_el_LQpure_wt;
               // h18_elel_LQint_wt = h_el_LQint_wt;
            }


            h1_elel_asym->SetMaximum(h1_elel_sym->GetMaximum()*1.2);
            //h1_elel_LQpure->SetMaximum(h1_elel_LQint->GetMaximum()*2);

            TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
            h1_elel_asym->SetTitle(el_title);
            h1_elel_asym->Draw("hist");
            h1_elel_sym->Draw("hist same ");

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_elel_asym, "Asym Template", "l");
            leg1->AddEntry(h1_elel_sym, "Sym Template", "l");

            leg1->Draw();

          //  c_elel1->Print(el_fname1);
            delete c_elel1;

            TCanvas *c_elel2 = new TCanvas("c_elel2", "Histograms", 200, 10, 900, 700);
            h1_elel_LQpure->SetTitle(el_title);
           // h1_elel_LQpure->Draw("hist");
           // h1_elel_mn->Draw("hist same");
            //h1_elel_alpha->Draw("hist same");
            //h1_elel_LQpure->SetTitle(el_title);
           // h1_elel_pl->Draw("hist same ");
            h1_elel_LQint->Draw("hist same");

             TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
           // leg2->AddEntry(h1_elel_pl, "Plus Template", "l");
           // leg2->AddEntry(h1_elel_mn, "Minus Template", "l");
           // leg2->AddEntry(h1_elel_alpha, "#alpha Template", "l");
            leg2->AddEntry(h1_elel_LQpure,"LQpure Template","l");
            leg2->AddEntry(h1_elel_LQint,"LQint Template","l");
           

            leg2->Draw();

            c_elel2->Print(el_fname2);
            delete c_elel2;
/*
            TCanvas *c_elel3 = new TCanvas("c_elel3", "Histograms", 200, 10, 900, 700);
            h_el_LQint_wt->SetTitle(el_title);
            h_el_LQint_wt->Draw("hist");
           h_el_LQpure_wt->Draw("hist same ");

            TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
           leg3->AddEntry(h_el_LQint_wt, "LQint weights*evtwt", "l");
            leg3->AddEntry(h_el_LQpure_wt, "LQpure weights*evtwt", "l");
            leg3->Draw();

         //   c_elel3->Print(el_fname3);
            delete c_elel3;
 */           
        }
    }
    char mu_1[100], mu_2[100], el_1[100], el_2[100];
    sprintf(mu_1, "Muons_LQpure_m%i.png",int(m_LQ));
    sprintf(mu_2, "Muons_LQint_m%i.png",int(m_LQ));
    sprintf(el_1, "Elecs_LQpure_m%i.png",int(m_LQ));
    sprintf(el_2, "Elecs_LQ_m%iint.png",int(m_LQ));
    //print lq pure weights for all years
  /*  h17_mumu_LQpure_wt->SetLineColor(kBlue);
    h18_mumu_LQpure_wt->SetLineColor(kGreen);

    TCanvas *c_mu_pwt = new TCanvas("c_mu_pwt", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQpure_wt->SetTitle("Muons: LQpure weights");
    h18_mumu_LQpure_wt->Draw("hist");
    h17_mumu_LQpure_wt->Draw("hist same");
    h16_mumu_LQpure_wt->Draw("hist same");
            

    TLegend *leg0 = new TLegend(x_start, y_start, x_end, y_end);
    leg0->AddEntry(h16_mumu_LQpure_wt, "LQpure16 weights*evtwt", "l");
    leg0->AddEntry(h17_mumu_LQpure_wt, "LQpure17 weights*evtwt", "l");
    leg0->AddEntry(h18_mumu_LQpure_wt, "LQpure18 weights*evtwt", "l");
    leg0->Draw();
            
    c_mu_pwt->Print("Misc_plots/template_plots2/Muons_LQpure_weights.png");
    delete c_mu_pwt;
    //print lq int weights for all years
    h17_mumu_LQint_wt->SetLineColor(kRed);
    h18_mumu_LQint_wt->SetLineColor(kGreen);

    TCanvas *c_mu_iwt = new TCanvas("c_mu_iwt", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQint_wt->SetTitle("Muons: LQint weights");
    h18_mumu_LQint_wt->Draw("hist");
    h17_mumu_LQint_wt->Draw("hist same");
    h16_mumu_LQint_wt->Draw("hist same");
            

    TLegend *leg = new TLegend(x_start, y_start, x_end, y_end);
    leg->AddEntry(h16_mumu_LQint_wt, "LQint16 weights*evtwt", "l");
    leg->AddEntry(h17_mumu_LQint_wt, "LQint17 weights*evtwt", "l");
    leg->AddEntry(h18_mumu_LQint_wt, "LQint18 weights*evtwt", "l");
    leg->Draw();
            
    c_mu_iwt->Print("Misc_plots/template_plots2/Muons_LQint_weights.png");
    delete c_mu_iwt;
   //print plus minus temps for all years
    
    h17_mumu_pl->SetLineColor(kBlue);
    h18_mumu_pl->SetLineColor(kGreen);
    h17_mumu_mn->SetLineColor(kBlue+5);
    h18_mumu_mn->SetLineColor(kGreen+5);

    TCanvas *c_mu_plmn = new TCanvas("c_mu_plmn", "Histograms", 200, 10, 900, 700);
    h16_mumu_pl->SetTitle("Muons: Plus and Minus templates");
    h16_mumu_pl->Draw("hist");
    h17_mumu_pl->Draw("hist same");
    h18_mumu_pl->Draw("hist same");
    h16_mumu_mn->Draw("hist same");
    h17_mumu_mn->Draw("hist same");
    h18_mumu_mn->Draw("hist same");

    TLegend *leg_ = new TLegend(x_start, y_start, x_end, y_end);
    leg_->AddEntry(h16_mumu_pl, "plus16", "l");
    leg_->AddEntry(h17_mumu_pl, "plus17", "l");
    leg_->AddEntry(h18_mumu_pl, "plus18", "l");
    leg_->AddEntry(h16_mumu_mn, "minus16", "l");
    leg_->AddEntry(h17_mumu_mn, "minus17", "l");
    leg_->AddEntry(h18_mumu_mn, "minus18", "l");
    leg_->Draw();
            
    c_mu_plmn->Print("Muons_PlMn.png");
    delete c_mu_plmn;
    */
    //print lq int for all years
    h17_mumu_LQint->SetLineColor(kBlue);
    h18_mumu_LQint->SetLineColor(kGreen+3);

    TCanvas *c_mu_lqi = new TCanvas("c_mu_lqi", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQint->SetTitle("Muons: LQint templates, m_LQ = %i",int(m_LQ));
    h18_mumu_LQint->Draw("hist");
    h17_mumu_LQint->Draw("hist same");
    h16_mumu_LQint->Draw("hist same");
            

    TLegend *leg_1 = new TLegend(x_start, y_start, x_end, y_end);
    leg_1->AddEntry(h16_mumu_LQint, "LQint16", "l");
    leg_1->AddEntry(h17_mumu_LQint, "LQint17", "l");
    leg_1->AddEntry(h18_mumu_LQint, "LQint18", "l");
    leg_1->Draw();
            
    c_mu_lqi->Print(mu_2);
    delete c_mu_lqi;

    //print lq pure for all years
    h17_mumu_LQpure->SetLineColor(kBlue);
    h18_mumu_LQpure->SetLineColor(kGreen+3);

    TCanvas *c_mu_lqp = new TCanvas("c_mu_lqp", "Histograms", 200, 10, 900, 700);
    h18_mumu_LQpure->SetTitle("Muons: LQpure templates, m_LQ = %i",int(m_LQ));
    h18_mumu_LQpure->Draw("hist");
    h17_mumu_LQpure->Draw("hist same");
    h16_mumu_LQpure->Draw("hist same");
            

    TLegend *leg_2 = new TLegend(x_start, y_start, x_end, y_end);
    leg_2->AddEntry(h16_mumu_LQpure, "LQpure16", "l");
    leg_2->AddEntry(h17_mumu_LQpure, "LQpure17", "l");
    leg_2->AddEntry(h18_mumu_LQpure, "LQpure18", "l");
    leg_2->Draw();
            
    c_mu_lqp->Print(mu_1);
    delete c_mu_lqp;

    //=============================================================================================
    //print lq pure weights for all years
/*
    h17_elel_LQpure_wt->SetLineColor(kBlue);
    h18_elel_LQpure_wt->SetLineColor(kGreen);

    TCanvas *c_el_pwt = new TCanvas("c_el_pwt", "Histograms", 200, 10, 900, 700);
    h18_elel_LQpure_wt->SetTitle("elecs: LQpure weights");
    h18_elel_LQpure_wt->Draw("hist");
    h17_elel_LQpure_wt->Draw("hist same");
    h16_elel_LQpure_wt->Draw("hist same");
            
    leg0->Draw();
            
    c_el_pwt->Print("Misc_plots/template_plots2/elecs_LQpure_weights.png");
    delete c_el_pwt;
    //print lq int weights for all years
    h17_elel_LQint_wt->SetLineColor(kRed);
    h18_elel_LQint_wt->SetLineColor(kGreen);

    TCanvas *c_el_iwt = new TCanvas("c_el_iwt", "Histograms", 200, 10, 900, 700);
    h18_elel_LQint_wt->SetTitle("elecs: LQint weights");
    h18_elel_LQint_wt->Draw("hist");
    h17_elel_LQint_wt->Draw("hist same");
    h16_elel_LQint_wt->Draw("hist same");

    leg->Draw();
            
    c_el_iwt->Print("Misc_plots/template_plots2/elecs_LQint_weights.png");
    delete c_el_iwt;
    /*
    //print plus minus temps for all years
    h17_elel_pl->SetLineColor(kBlue);
    h18_elel_pl->SetLineColor(kGreen);
    h17_elel_mn->SetLineColor(kBlue+5);
    h18_elel_mn->SetLineColor(kGreen+5);

    TCanvas *c_el_plmn = new TCanvas("c_el_plmn", "Histograms", 200, 10, 900, 700);
    h16_elel_pl->SetTitle("elecs: Plus and Minus templates");
    h16_elel_pl->Draw("hist");
    h17_elel_pl->Draw("hist same");
    h18_elel_pl->Draw("hist same");
    h16_elel_mn->Draw("hist same");
    h17_elel_mn->Draw("hist same");
    h18_elel_mn->Draw("hist same");

    leg_->Draw();
            
    c_el_plmn->Print("elecs_PlMn.png");
    delete c_el_plmn;
    */
    //print lq int for all years
    h17_elel_LQint->SetLineColor(kBlue);
    h18_elel_LQint->SetLineColor(kGreen+3);

    TCanvas *c_el_lqi = new TCanvas("c_el_lqi", "Histograms", 200, 10, 900, 700);
    h18_elel_LQint->SetTitle("Elecs: LQint templates, m_LQ = %i",int(m_LQ));
    h18_elel_LQint->Draw("hist");
    h17_elel_LQint->Draw("hist same");
    h16_elel_LQint->Draw("hist same");
            
    leg_1->Draw();
            
    c_el_lqi->Print(el_2);
    delete c_el_lqi;

    //print lq pure for all years
    h17_elel_LQpure->SetLineColor(kBlue);
    h18_elel_LQpure->SetLineColor(kGreen+3);

    TCanvas *c_el_lqp = new TCanvas("c_el_lqp", "Histograms", 200, 10, 900, 700);
    h18_elel_LQpure->SetTitle("Elecs: LQpure templates, m_LQ = %i",int(m_LQ));
    h18_elel_LQpure->Draw("hist");
    h17_elel_LQpure->Draw("hist same");
    h16_elel_LQpure->Draw("hist same");
    
    leg_2->Draw();
            
    c_el_lqp->Print(el_1);
    delete c_el_lqp;

    }
}




