/*
vector<string>  sys_labels_uncorr = 
       // { "_BTAGCOR", "_BTAGUNCOR", "_BTAGLIGHT",  "_METJER", "_METJEC", "_METHEM", 
       { "_muPref","_prefire", "_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear", "_muRC", "_Pu", 
            "_muHLTBAR", "_muIDBAR", "_muISOBAR",  "_muHLTEND", "_muIDEND", "_muISOEND",  "_muIDSYS", "_muISOSYS",  
            "_elHLTBARPTHIGH", "_elIDBARPTHIGH", "_elRECOBARPTHIGH", "_elHLTENDPTHIGH", "_elIDENDPTHIGH", "_elRECOENDPTHIGH",
            "_elHLTBARPTLOW", "_elIDBARPTLOW", "_elRECOBARPTLOW", "_elHLTENDPTLOW", "_elIDENDPTLOW", "_elRECOENDPTLOW",
           // "_ptrw1b", "_ptrw2b", "_ptrw3b", "_ptrw4b", "_ptrw5b", "_ptrw6b", "_ptrw7b",
            "_emucostrw1b", "_emucostrw2b", "_emucostrw3b", "_emucostrw4b",
            "_elfakesrw1b", "_elfakesrw2b", "_elfakesrw3b", "_elfakesrw4b",
            "_mufakesrw1b", "_mufakesrw2b", "_mufakesrw3b", "_mufakesrw4b",
            "_RENORM", "_FAC", "_REFAC","_alphaS",
        };
*/


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
    const int num_sys = 4;
    //string sys_array[num_sys] = {"_RENORM","_REFAC","_FAC","_muRC"};
    string sys_array[num_sys] = {"_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear"};
    string sys_l = "momentum_scale", plot_t = "Electron Momentum Scale";
    const char *sys_label = sys_l.c_str();
    const char *plot_title = plot_t.c_str();
    char *plot_dir = "AN_plots/Systematics/UpDown";
    float m_LQ = 2500.;
    Double_t alpha_denom = (amc_alpha[5]+amc_alpha[6]+amc_alpha[7])/3.;

    double afb = 0.6;

    bool use_xf = false;

    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xf){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }

    for(int flag_q = 1; flag_q <=2; flag_q++){

        bool do_bkg = false;
        bool do_qcd = false;
        bool do_electrons = true;
        bool do_muons = false;
        bool vec = false;
        bool make_ud = true;
        float yLQ = 1.0;
        char *date;
        if(flag_q==2) {date = "u"; if(vec) date = "u_vec";}
        else {date = "d"; if(vec) date = "d_vec";}
        char mu_title[100], el_title[100];
        char mu_fname1[100],  el_fname1[100];
        sprintf(mu_fname1, "%s/mumu_yLQ%.1f_%s_chk_%s.png", plot_dir, yLQ, sys_label, date);
        sprintf(el_fname1, "%s/ee_yLQ%.1f_%s_chk_%s.png", plot_dir, yLQ, sys_label, date);
        

        

        
        TH1F *h1_elel_plain_comb, *h1_elel_sys_up_comb, *h1_elel_sys_down_comb;
        TH1F *h1_mumu_plain_comb, *h1_mumu_sys_up_comb, *h1_mumu_sys_down_comb;
        TH1F *h1_elel_bkg, *h1_mumu_bkg, *h1_elel_bkg_up, *h1_elel_bkg_down, *h1_mumu_bkg_up, *h1_mumu_bkg_down;
        TH1F *h1_elel_qcd, *h1_mumu_qcd, *h1_elel_qcd_up, *h1_elel_qcd_down, *h1_mumu_qcd_up, *h1_mumu_qcd_down;
        
        for(int year = 2016; year <= 2016; year++){

            init(year);
            setup_all_SFs(year);


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

            TH3F * h_mumu_qcd = new TH3F("mumu_qcd", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
            TH3F * h_mumu_qcd_up = new TH3F("mumu_qcd_up", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
            TH3F * h_mumu_qcd_down = new TH3F("mumu_qcd_down", "", n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);


            bool ss = false;




            TTree *elel_ts[3] = {t_elel_ttbar, t_elel_wt, t_elel_diboson};
            TTree *mumu_ts[3] = {t_mumu_ttbar, t_mumu_wt, t_mumu_diboson};
            //char mu_title[100], el_title[100];


            if(do_muons){
                printf("Making mumu temps \n");
                one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_plain, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS,make_ud,  use_xf, "");
                TH1F* h1_mumu_plain = convert3d(h_mumu_plain);
                h1_mumu_plain->SetLineColor(kBlack);
                h1_mumu_plain->SetLineWidth(2);

                //just for initializing TH1F
                TH1F* h1_mumu_sys_up = convert3d(h_mumu_plain);
                TH1F* h1_mumu_sys_down = convert3d(h_mumu_plain);
                int nbins = h1_mumu_plain->GetNbinsX();
                float up_diff[60] = {0.};
                float down_diff[60] = {0.};
                float up = 0, down = 0;

                for(int i = 0; i< num_sys; i++){

                    TH3F * h_mumu_sys_up = new TH3F("mumu_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
                    TH3F * h_mumu_sys_down = new TH3F("mumu_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
                    const char *sys = sys_array[i].c_str();
                    string sys_up = string(sys) + string("Up");
                    string sys_down = string(sys) + string("Down");
                    one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_sys_up, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS, make_ud, use_xf, sys_up);
                    one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_sys_down, year, m_LQ, yLQ , flag_q, vec, FLAG_MUONS, make_ud, use_xf, sys_down);

                    h1_mumu_sys_up = convert3d(h_mumu_sys_up);
                    h1_mumu_sys_down = convert3d(h_mumu_sys_down);

                    for(int j = 1; j <= h1_mumu_sys_up->GetNbinsX(); j++){
                        up_diff[j] += h1_mumu_sys_up->GetBinContent(j) - h1_mumu_plain->GetBinContent(j);
                        down_diff[j] += h1_mumu_plain->GetBinContent(j) - h1_mumu_sys_down->GetBinContent(j);
                    }
                }
                for(int i = 1; i <= h1_mumu_plain->GetNbinsX(); i++){

                    up = h1_mumu_plain->GetBinContent(i) + up_diff[i];
                    down = h1_mumu_plain->GetBinContent(i) - down_diff[i];
                    h1_mumu_sys_up->SetBinContent(i, up);
                    h1_mumu_sys_down->SetBinContent(i, down);
                }
                
                h1_mumu_sys_up->SetLineColor(kRed);
                h1_mumu_sys_up->SetLineWidth(2);
                h1_mumu_sys_down->SetLineColor(kBlue);
                h1_mumu_sys_down->SetLineWidth(2);
                
                if(year==2016){
                   h1_mumu_plain_comb = (TH1F*)h1_mumu_plain->Clone();
                   h1_mumu_sys_up_comb = (TH1F*)h1_mumu_sys_up->Clone();
                   h1_mumu_sys_down_comb = (TH1F*)h1_mumu_sys_down->Clone(); 
               } else {

                h1_mumu_plain_comb->Add(h1_mumu_plain);
                h1_mumu_sys_up_comb->Add(h1_mumu_sys_up);
                h1_mumu_sys_down_comb->Add(h1_mumu_sys_down);
                
            }



            printf("MuMu: nom %.0f, up %.0f, down %.0f \n", h1_mumu_plain->Integral(), h1_mumu_sys_up->Integral(), h1_mumu_sys_down->Integral());
            /*
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
            */


        }


        if(do_electrons){
            printf("Making elel temps \n");
            one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_plain, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS,make_ud,  use_xf, "");
            TH1F* h1_elel_plain = convert3d(h_elel_plain);
            h1_elel_plain->SetLineColor(kBlack);
            h1_elel_plain->SetLineWidth(2);

            //just for initializing TH1F
            TH1F* h1_elel_sys_up = convert3d(h_elel_plain);
            TH1F* h1_elel_sys_down = convert3d(h_elel_plain);
            int nbins = h1_elel_plain->GetNbinsX();
            float up_diff[60] = {0.};
            float down_diff[60] = {0.};
            float up = 0, down = 0;

            for(int i = 0; i< num_sys; i++){

                TH3F * h_elel_sys_up = new TH3F("elel_up", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
                TH3F * h_elel_sys_down = new TH3F("elel_down", "", n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
                const char *sys = sys_array[i].c_str();
                string sys_up = string(sys) + string("Up");
                string sys_down = string(sys) + string("Down");
                one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_sys_up, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS, make_ud, use_xf, sys_up);
                one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_sys_down, year, m_LQ, yLQ , flag_q, vec, FLAG_ELECTRONS, make_ud, use_xf, sys_down);

                TH1F* h1_elel_sys_up = convert3d(h_elel_sys_up);
                TH1F* h1_elel_sys_down = convert3d(h_elel_sys_down);

                for(int j = 1; j <= h1_elel_sys_up->GetNbinsX(); j++){
                    up_diff[j] += h1_elel_sys_up->GetBinContent(j) - h1_elel_plain->GetBinContent(j);
                    down_diff[j] += h1_elel_plain->GetBinContent(j) - h1_elel_sys_down->GetBinContent(j);
                }
            }
            for(int i = 1; i <= h1_elel_plain->GetNbinsX(); i++){

                up = h1_elel_plain->GetBinContent(i) + up_diff[i];
                down = h1_elel_plain->GetBinContent(i) - down_diff[i];
                h1_elel_sys_up->SetBinContent(i, up);
                h1_elel_sys_down->SetBinContent(i, down);
            }

            h1_elel_sys_up->SetLineColor(kRed);
            h1_elel_sys_up->SetLineWidth(2);
            h1_elel_sys_down->SetLineColor(kBlue);
            h1_elel_sys_down->SetLineWidth(2);

            if(year==2016){
               h1_elel_plain_comb = (TH1F*)h1_elel_plain->Clone();
               h1_elel_sys_up_comb = (TH1F*)h1_elel_sys_up->Clone();
               h1_elel_sys_down_comb = (TH1F*)h1_elel_sys_down->Clone(); 
           } else {

            h1_elel_plain_comb->Add(h1_elel_plain);
            h1_elel_sys_up_comb->Add(h1_elel_sys_up);
            h1_elel_sys_down_comb->Add(h1_elel_sys_down);

        }



        printf("MuMu: nom %.0f, up %.0f, down %.0f \n", h1_elel_plain->Integral(), h1_elel_sys_up->Integral(), h1_elel_sys_down->Integral());
        /*
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
        */

    }

}



if(do_muons){

    if(flag_q == 2) {
        sprintf(mu_title, "%s, m_{LQ} = %.1f TeV, y_{#mu u} = %.1f, Channel: #mu-u 2016,2017,2018", plot_title, m_LQ/1000, yLQ);
        if(vec) sprintf(mu_title, "%s, m_{LQ} = %.1f TeV, g_{#mu u} = %.1f, Channel: #mu-u-vec 2016,2017,2018", plot_title, m_LQ/1000, yLQ);
    }
    else {
        sprintf(mu_title, "%s, m_{LQ} = %.1f TeV, y_{#mu d} = %.1f, Channel: #mu-d : 2016,2017,2018 ", plot_title, m_LQ/1000, yLQ);
        if(vec) sprintf(mu_title, "%s, m_{LQ} = %.1f TeV, g_{#mu d} = %.1f, Channel: #mu-d-vec : 2016,2017,2018 ", plot_title, m_LQ/1000, yLQ);
    }
    TCanvas *c_mumu1 = new TCanvas("c_mumu", "Muons", 200, 10, 900, 700);
    TPad *pad1 = new TPad("p1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    h1_mumu_sys_up_comb->SetMaximum(h1_mumu_sys_up_comb->GetMaximum()*1.5);
    h1_mumu_sys_up_comb->SetTitle(mu_title);
    h1_mumu_sys_up_comb->GetYaxis()->SetLabelSize(0.035);
    h1_mumu_sys_up_comb->Draw("hist");
    h1_mumu_plain_comb->Draw("hist same");
    h1_mumu_sys_down_comb->Draw("hist same");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h1_mumu_plain_comb, "Nominal LQ Signal template", "l");
    leg1->AddEntry(h1_mumu_sys_up_comb, "Sys Up Template", "l");
    leg1->AddEntry(h1_mumu_sys_down_comb, "Sys Down Template", "l");

    c_mumu1->cd();
    TPad *pad2 = new TPad("p2", "pad2", 0.,0,.98,0.3);
            //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    TH1F *ratio_up, *ratio_down;

    ratio_up = (TH1F *) h1_mumu_sys_up_comb->Clone("h_ratio_up");
    ratio_up->Sumw2();
    ratio_up->SetStats(0);
    ratio_up->Divide(h1_mumu_plain_comb);

    ratio_down = (TH1F *) h1_mumu_sys_down_comb->Clone("h_ratio_down");
    ratio_down->Sumw2();
    ratio_down->SetStats(0);
    ratio_down->Divide(h1_mumu_plain_comb);

            //ratio_up->SetMarkerStyle(21);
    ratio_down->SetMinimum(0.8);
    ratio_down->SetMaximum(1.2);
    ratio_down->SetTitle("");
    ratio_down->SetLineColor(kBlue);
    ratio_down->GetXaxis()->SetLabelSize(0.1);
    ratio_down->GetYaxis()->SetLabelSize(0.1);
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
    if(flag_q == 2) {
        sprintf(el_title, "%s, m_{LQ} = %.1f TeV, y_{eu} = %.1f, Channel: e-u 2016,2017,2018", plot_title, m_LQ/1000, yLQ);
        if(vec) sprintf(el_title, "%s, m_{LQ} = %.1f TeV, g_{eu} = %.1f, Channel: e-u-vec 2016,2017,2018", plot_title, m_LQ/1000, yLQ);
    }
    else {
        sprintf(el_title, "%s, m_{LQ} = %.1f TeV, y_{ed} = %.1f, Channel: e-d : 2016,2017,2018 ", plot_title, m_LQ/1000, yLQ);
        if(vec) sprintf(el_title, "%s, m_{LQ} = %.1f TeV, g_{ed} = %.1f, Channel: e-d-vec : 2016,2017,2018 ", plot_title, m_LQ/1000, yLQ);
    }
    TCanvas *c_elel1 = new TCanvas("c_mumu", "Muons", 200, 10, 900, 700);
    TPad *pad1 = new TPad("p1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    h1_elel_sys_up_comb->SetMaximum(h1_elel_sys_up_comb->GetMaximum()*1.5);
    h1_elel_sys_up_comb->SetTitle(el_title);
    h1_elel_sys_up_comb->GetYaxis()->SetLabelSize(0.035);
    h1_elel_sys_up_comb->Draw("hist");
    h1_elel_plain_comb->Draw("hist same");
    h1_elel_sys_down_comb->Draw("hist same");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h1_elel_plain_comb, "Nominal LQ Signal template", "l");
    leg1->AddEntry(h1_elel_sys_up_comb, "Sys Up Template", "l");
    leg1->AddEntry(h1_elel_sys_down_comb, "Sys Down Template", "l");

    c_elel1->cd();
    TPad *pad2 = new TPad("p2", "pad2", 0.,0,.98,0.3);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    TH1F *ratio_up, *ratio_down;

    ratio_up = (TH1F *) h1_elel_sys_up_comb->Clone("h_ratio_up");
    ratio_up->Sumw2();
    ratio_up->SetStats(0);
    ratio_up->Divide(h1_elel_plain_comb);

    ratio_down = (TH1F *) h1_elel_sys_down_comb->Clone("h_ratio_down");
    ratio_down->Sumw2();
    ratio_down->SetStats(0);
    ratio_down->Divide(h1_elel_plain_comb);
    ratio_down->SetMinimum(0.8);
    ratio_down->SetMaximum(1.2);
    ratio_down->SetTitle("");
    ratio_down->SetLineColor(kBlue);
    ratio_down->GetXaxis()->SetLabelSize(0.1);
    ratio_down->GetYaxis()->SetLabelSize(0.1);
    ratio_down->Draw("hist");
           // ratio_down->SetMarkerStyle(21);
    ratio_up->SetLineColor(kBlue);
    ratio_up->Draw("hist same");

    c_elel1->cd(); 


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
