#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"
#include "../../utils/PlotUtils.C"

TH3F *h_elel_asym, *h_elel_sym, *h_elel_alpha, *h_elel_LQpure,*h_elel_LQint,*h_elel_back,  *h_elel_dy_gg, *h_elel_data, *h_elel_mc, *h_elel_qcd, *h_elel_gam;
TH3F *h_elel_mc_count, *h_elel_sym_count;
TH3F *h_mumu_asym, *h_mumu_sym, *h_mumu_alpha, *h_mumu_LQpure,*h_mumu_LQint,*h_mumu_back,  *h_mumu_dy_gg, *h_mumu_data, *h_mumu_mc, *h_mumu_qcd, *h_mumu_gam;



void LQ_sys_check(){




    gStyle->SetOptStat(0);
    gROOT->SetBatch(1); 
    const int num_sys = 2;
    string sys_array[num_sys] = {"_elScaleGain","_ptrw7b"};
    for(int year = 2016; year <= 2018; year++){
	for(int flag_q = 1; flag_q <=2; flag_q++){
        init(year);

        float m_LQ = 2000.;
        char *plot_dir = "Misc_plots";
	char *date;
        if(flag_q==1) date = "010923_u";
	else date = "010923_d";
        //char *sys = "_";
        bool do_bkg = false;
	bool do_qcd = false;
        bool do_electrons = false;
        bool do_muons = true;
        bool vec = false;
	bool make_ud = true;
	if(!make_ud) date = "101322_c";
        //int flag_q = 2;
        float yLQ = 1.0;

        setup_all_SFs(year);
	for(int i = 0; i< num_sys; i++){

	const char *sys = sys_array[i].c_str();
        string sys_up = string(sys) + string("Up");
        string sys_down = string(sys) + string("Down");

        Double_t alpha_denom = (amc_alpha[5]+amc_alpha[6]+amc_alpha[7])/3.;
        //double m_low = m_bins[i];
        //double m_high = m_bins[i+1];
        double afb = 0.6;


        char mu_fname1[100],  el_fname1[100];

        sprintf(mu_fname1, "%s/mumu%i_yLQ%.1f_%s_chk_%s.png", plot_dir, year, yLQ, sys, date);
        sprintf(el_fname1, "%s/ee%i_yLQ%.1f_%s_chk_%s.png", plot_dir, year, yLQ, sys, date);

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
        TH3F * h_elel_qcd = new TH3F("elel_qcd", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_qcd_up = new TH3F("elel_qcd_up", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_elel_qcd_down = new TH3F("elel_qcd_down", "",n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);

        TH3F * h_mumu_bkg = new TH3F("mumu_bkg", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_bkg_up = new TH3F("mumu_bkg_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_bkg_down = new TH3F("mumu_bkg_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_plain = new TH3F("mumu_plain", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_sys_up = new TH3F("mumu_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_sys_down = new TH3F("mumu_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_qcd = new TH3F("mumu_qcd", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_qcd_up = new TH3F("mumu_qcd_up", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH3F * h_mumu_qcd_down = new TH3F("mumu_qcd_down", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);


        bool ss = false;




        TTree *elel_ts[3] = {t_elel_ttbar, t_elel_wt, t_elel_diboson};
        TTree *mumu_ts[3] = {t_mumu_ttbar, t_mumu_wt, t_mumu_diboson};
        char mu_title[100], el_title[100];

        TH1F *h1_elel_bkg, *h1_mumu_bkg, *h1_elel_bkg_up, *h1_elel_bkg_down, *h1_mumu_bkg_up, *h1_mumu_bkg_down;
	TH1F *h1_elel_qcd, *h1_mumu_qcd, *h1_elel_qcd_up, *h1_elel_qcd_down, *h1_mumu_qcd_up, *h1_mumu_qcd_down;
        if(do_muons){
            printf("Making mumu temps \n");
            one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_plain, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS,make_ud,  use_xf, "");
            one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_sys_up, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS, make_ud, use_xf, sys_up);
            one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_sys_down, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS, make_ud, use_xf, sys_down);
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
                bool emu_reweight = false;
                gen_combined_background_template(3, mumu_ts, h_mumu_bkg, year, FLAG_MUONS,  ss, use_xf, emu_reweight, "");
                gen_combined_background_template(3, mumu_ts, h_mumu_bkg_up, year, FLAG_MUONS,  ss, use_xf, emu_reweight, sys_up);
                gen_combined_background_template(3, mumu_ts, h_mumu_bkg_down, year, FLAG_MUONS,  ss, use_xf, emu_reweight, sys_down);
                
                symmetrize3d(h_mumu_bkg);
                symmetrize3d(h_mumu_bkg_up);
                symmetrize3d(h_mumu_bkg_down);
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
             if(do_qcd){
                bool incl_ss = true;
                bool ss_binning = false;
                gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, 
                        FLAG_MUONS, incl_ss, ss_binning, use_xF, "");
                gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd_up, year,  
                        FLAG_MUONS, incl_ss, ss_binning, use_xF, sys_up);
                gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd_down, year, 
                        FLAG_MUONS, incl_ss, ss_binning, use_xF, sys_down);

                symmetrize3d(h_mumu_qcd);
                symmetrize3d(h_mumu_qcd_up);
                symmetrize3d(h_mumu_qcd_down);


                h1_mumu_qcd = convert3d(h_mumu_qcd);
                h1_mumu_qcd_up = convert3d(h_mumu_qcd_up);
                h1_mumu_qcd_down = convert3d(h_mumu_qcd_down);
                

                h1_mumu_qcd->SetLineColor(kGray);
                h1_mumu_qcd->SetLineWidth(2);
                h1_mumu_qcd_up->SetLineColor(kOrange + 3);
                h1_mumu_qcd_up->SetLineWidth(2);
                h1_mumu_qcd_down->SetLineColor(kOrange +1);
                h1_mumu_qcd_down->SetLineWidth(2);
                printf("mumu fakes: nom %.0f, up %.0f, down %.0f \n", h_mumu_qcd->Integral(), h_mumu_qcd_up->Integral(), h_mumu_qcd_down->Integral());
                h1_mumu_qcd->Print("range");
            }
		unzero_bins(h_mumu_plain);
		unzero_bins(h_mumu_sys_down);

            sprintf(mu_title, "SM DY + LQ_um + all bkgs  %s, m_LQ = %i GeV, year = %i, yLQ = %.1f ", sys, int(m_LQ), year, yLQ);
            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Muons", 200, 10, 900, 700);
            TPad *pad1 = new TPad("p1", "pad1", 0.,0.3,0.98,1.);
            pad1->SetBottomMargin(0);
            pad1->Draw();
            pad1->cd();
            h1_mumu_plain->SetTitle(mu_title);
            h1_mumu_plain->Draw("hist");
            h1_mumu_sys_up->Draw("hist same");
            h1_mumu_sys_down->Draw("hist same");

            gStyle->SetLegendBorderSize(0);
            TLegend *leg1 = new TLegend(0.7, 0.75, 0.9, 0.9);
            leg1->AddEntry(h1_mumu_plain, "Nominal Template", "l");
            leg1->AddEntry(h1_mumu_sys_up, "Sys Up Template", "l");
            leg1->AddEntry(h1_mumu_sys_down, "Sys Down Template", "l");

            c_mumu1->cd();
            TPad *pad2 = new TPad("p2", "pad2", 0.,0,.98,0.3);
            //pad2->SetTopMargin(0);
            pad2->SetBottomMargin(0.2);
            pad2->SetGridy();
            pad2->Draw();
            pad2->cd();
            TH1F *ratio_up, *ratio_down;
    
            ratio_up = (TH1F *) h1_mumu_sys_up->Clone("h_ratio_up");
            ratio_up->Sumw2();
            ratio_up->SetStats(0);
            ratio_up->Divide(h1_mumu_plain);

            ratio_down = (TH1F *) h1_mumu_sys_down->Clone("h_ratio_down");
            ratio_down->Sumw2();
            ratio_down->SetStats(0);
            ratio_down->Divide(h1_mumu_plain);

            //ratio_up->SetMarkerStyle(21);
            ratio_down->SetMinimum(0.94);
	    ratio_down->SetMaximum(1.06);
	    ratio_down->SetTitle("");
            ratio_down->SetLineColor(kGreen+3);
            ratio_down->Draw("hist");
           // ratio_down->SetMarkerStyle(21);
            ratio_up->SetLineColor(kBlue);
            ratio_up->Draw("hist same");
            
            c_mumu1->cd();

            if(do_bkg){
                h1_mumu_bkg->Draw("hist same");
                h1_mumu_bkg_up->Draw("hist same");
                h1_mumu_bkg_down->Draw("hist same");
            }
            if(do_qcd){
                h1_mumu_qcd->Draw("hist same");
                h1_mumu_qcd_up->Draw("hist same");
                h1_mumu_qcd_down->Draw("hist same");
            }
            
            if(do_bkg){
                leg1->AddEntry(h1_mumu_bkg, "Nominal Bkg Template", "l");
                leg1->AddEntry(h1_mumu_bkg_up, "Sys Up Bkg Template", "l");
                leg1->AddEntry(h1_mumu_bkg_down, "Sys Down Bkg Template", "l");
            }
             if(do_qcd){
                leg1->AddEntry(h1_mumu_qcd, "Nominal qcd Template", "l");
                leg1->AddEntry(h1_mumu_qcd_up, "Sys Up qcd Template", "l");
                leg1->AddEntry(h1_mumu_qcd_down, "Sys Down qcd Template", "l");
            }
            leg1->Draw();

            c_mumu1->Print(mu_fname1);
        }


        if(do_electrons){
            printf("Making elel temps \n");

            one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_plain, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS, make_ud, use_xf, "");
            one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_sys_up, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS, make_ud, use_xf, sys_up);
            one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_sys_down, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS, make_ud, use_xf, sys_down);
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
                bool emu_reweight = false;
                gen_combined_background_template(3, elel_ts, h_elel_bkg, year, FLAG_ELECTRONS,  ss, use_xf,  emu_reweight,"");
                gen_combined_background_template(3, elel_ts, h_elel_bkg_up, year, FLAG_ELECTRONS,  ss, use_xf, emu_reweight, sys_up);
                gen_combined_background_template(3, elel_ts, h_elel_bkg_down, year, FLAG_ELECTRONS,  ss, use_xf,emu_reweight, sys_down);
                
                symmetrize3d(h_elel_bkg);
                symmetrize3d(h_elel_bkg_up);
                symmetrize3d(h_elel_bkg_down);
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
            if(do_qcd){
                bool incl_ss = true;
                bool ss_binning = false;
                gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, 
                        FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, "");
                gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd_up, year, 
                        FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, sys_up);
                gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd_down, year, 
                        FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, sys_down);

                symmetrize3d(h_elel_qcd);
                symmetrize3d(h_elel_qcd_up);
                symmetrize3d(h_elel_qcd_down);


                h1_elel_qcd = convert3d(h_elel_qcd);
                h1_elel_qcd_up = convert3d(h_elel_qcd_up);
                h1_elel_qcd_down = convert3d(h_elel_qcd_down);

                h1_elel_qcd->SetLineColor(kGray);
                h1_elel_qcd->SetLineWidth(2);
                h1_elel_qcd_up->SetLineColor(kOrange + 3);
                h1_elel_qcd_up->SetLineWidth(2);
                h1_elel_qcd_down->SetLineColor(kOrange +1);
                h1_elel_qcd_down->SetLineWidth(2);
                printf("elel fakes: nom %.0f, up %.0f, down %.0f \n", h_elel_qcd->Integral(), h_elel_qcd_up->Integral(), h_elel_qcd_down->Integral());
            }


            sprintf(el_title, "SM DY + LQ_ue + all bkgs %s, m_LQ = %i GeV, year = %i, yLQ = %.1f ", sys, int(m_LQ), year, yLQ);
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
             if(do_qcd){
                h1_elel_qcd->Draw("hist same");
                h1_elel_qcd_up->Draw("hist same");
                h1_elel_qcd_down->Draw("hist same");
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
               if(do_qcd){
                leg1->AddEntry(h1_elel_qcd, "Nominal qcd Template", "l");
                leg1->AddEntry(h1_elel_qcd_up, "Sys Up qcd Template", "l");
                leg1->AddEntry(h1_elel_qcd_down, "Sys Down qcd Template", "l");
            }
            leg1->Draw();

            c_elel1->Print(el_fname1);

        }

    }


}

}

}
