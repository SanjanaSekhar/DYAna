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
#include "results.h"


void draw_graph_ratio(TGraphErrors *g_data, TGraphErrors *g_mc, TGraphErrors *g_ratio, char label[10], std::string title){

    g_mc->SetMarkerColor(kBlue);
    g_mc->SetLineColor(kBlue);
    g_mc->SetLineWidth(4);


    g_data->SetMarkerColor(kBlack);
    g_data->SetLineColor(kBlack);

    g_data->SetMarkerStyle(kFullSquare);

    g_data->SetLineWidth(2);




    TCanvas *c_m = new TCanvas(label, label, 200, 10, 1000, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0.012);
    pad1->Draw();
    pad1->cd();

    g_mc->Draw("ALP");
    g_mc->GetYaxis()->SetRangeUser(0.0, 0.75);
    g_mc->GetXaxis()->SetLimits(100., 1400.);
    g_mc->Draw("ALP");
    g_data->Draw("PE same");

    g_mc->SetTitle(title.c_str());


    g_mc->GetYaxis()->SetTitle("Forward Backward Asymmetry");
    g_mc->GetYaxis()->SetNdivisions(505);
    g_mc->GetYaxis()->SetTitleSize(20);
    g_mc->GetYaxis()->SetTitleFont(43);
    //g_mc->GetYaxis()->SetTitleOffset(1.2);
    g_mc->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_mc->GetYaxis()->SetLabelSize(15);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.2, 0.2);
    leg1->AddEntry(g_mc, "aMC@NLO", "l");
    leg1->AddEntry(g_data, "Fitted AFB * Dilution", "p");
    leg1->Draw();

    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    g_ratio->Draw("APE");
    g_ratio->GetYaxis()->SetRangeUser(0.8, 1.20);
    g_ratio->GetXaxis()->SetLimits(100., 1400.);
    g_ratio->Draw("APE");


    g_ratio->GetYaxis()->SetTitle("Ratio");
    g_ratio->GetYaxis()->SetNdivisions(505);
    g_ratio->GetYaxis()->SetTitleSize(20);
    g_ratio->GetYaxis()->SetTitleFont(43);
    g_ratio->GetYaxis()->SetTitleOffset(1.2);
    g_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetYaxis()->SetLabelSize(15);
    // X axis g_ratio plot settings
    g_ratio->GetXaxis()->SetTitle("M_{ll} (GeV)");
    g_ratio->GetXaxis()->SetTitleSize(20);
    g_ratio->GetXaxis()->SetTitleFont(43);
    g_ratio->GetXaxis()->SetTitleOffset(3.);
    g_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetXaxis()->SetLabelSize(20);
    int iPeriod = -1; 
    writeExtraText = true;
    extraTextFont = cmsTextFont;
    extraText = TString(title);
    CMS_lumi(pad1, iPeriod, 33 );
    c_m->Update();
}




void fill_fb_dilu(TTree *t_dy, TH2D *h_for, TH2D *h_back, TH2D *h_for_st, TH2D *h_back_st, 
        TH2D *h_f_corr, TH2D *h_f_inc, TH2D *h_b_corr, TH2D *h_b_inc){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    double gen_weight, m, cost, cost_st;
    t_dy->SetBranchAddress("gen_p", &gen_lep_p);
    t_dy->SetBranchAddress("gen_m", &gen_lep_m);
    t_dy->SetBranchAddress("m", &m);
    t_dy->SetBranchAddress("cost", &cost);
    t_dy->SetBranchAddress("cost_st", &cost_st);
    t_dy->SetBranchAddress("gen_weight", &gen_weight);

    double m_low = 150.;




    int nEvents=0;

    for (int i=0; i<t_dy->GetEntries(); i++) {
        t_dy->GetEntry(i);
        if(m > m_low){
            cm = *gen_lep_p + *gen_lep_m;
            nEvents++;

            // if same sign > 0
            double mult = cost * cost_st;

            Double_t y = abs(cm.Rapidity());

            if(cost > 0.){
                h_for->Fill(y, m, gen_weight);
            }
            else{
                h_back->Fill(y, m, gen_weight);
            }

            if(cost_st > 0.){
                if(mult > 0.) h_f_corr->Fill(y, m, gen_weight);
                else h_f_inc->Fill(y, m, gen_weight);
                h_for_st->Fill(y, m, gen_weight);
            }
            else{
                if(mult > 0.) h_b_corr->Fill(y, m, gen_weight);
                else h_b_inc->Fill(y, m, gen_weight);
                h_back_st->Fill(y, m, gen_weight);
            }

        }
    }
    printf("selected %i events \n", nEvents);

}


void plot_counting_cmp(){
    setTDRStyle();
    TFile *f_gen = TFile::Open("Results/DY_gen_level_nov8.root");

    TTree *t_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_el = (TTree *) f_gen->Get("T_gen_el");


    const int n_eta_bins = 4;
    double eta_bins[] = {0., 1.0, 1.25, 1.5, 2.4};

    double m_bins[] = {150, 171, 200,  250, 320, 510, 700, 1000, 14000};

    TH2D *h_for = new TH2D("h_for", "Foward", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_back = new TH2D("h_back", "Backward", n_eta_bins, eta_bins, n_m_bins, m_bins);

    TH2D *h_for_st = new TH2D("h_for_st", "Foward", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_back_st = new TH2D("h_back_st", "Backward", n_eta_bins, eta_bins, n_m_bins, m_bins);

    TH2D *h_f_corr = new TH2D("h_f_corr", "Correct sign", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_f_inc = new TH2D("h_f_inc", "Incorrect sign", n_eta_bins, eta_bins, n_m_bins, m_bins);

    TH2D *h_b_corr = new TH2D("h_b_corr", "Correct sign", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_b_inc = new TH2D("h_b_inc", "Incorrect sign", n_eta_bins, eta_bins, n_m_bins, m_bins);

    std::vector<TGraphErrors *> g_mc, g_fitted, g_ratio;
    fill_fb_dilu(t_mu, h_for, h_back, h_for_st, h_back_st, h_f_corr, h_f_inc, h_b_corr, h_b_inc);
    fill_fb_dilu(t_el, h_for, h_back, h_for_st, h_back_st, h_f_corr, h_f_inc, h_b_corr, h_b_inc);

    double afb_mc[8], afb_fit[8], afb_fit_unc[8], ratio[8], ratio_unc[8];
    for(int i=1; i<= n_eta_bins; i++){
        for(int j=1; j<= n_m_bins; j++){
            double n_for = h_for->GetBinContent(i,j);
            double n_back = h_back->GetBinContent(i,j);

            double n_for_st =  h_for_st->GetBinContent(i,j);
            double n_back_st = h_back_st->GetBinContent(i,j);
            
            double n_f_corr = h_f_corr->GetBinContent(i,j);
            double n_f_inc = h_f_inc->GetBinContent(i,j);

            double n_b_corr = h_b_corr->GetBinContent(i,j);
            double n_b_inc = h_b_inc->GetBinContent(i,j);
            
            double dilu_f = std::min((n_f_corr - n_f_inc)/ (n_f_corr + n_f_inc), 1.);
            double dilu_b = std::min((n_b_corr - n_b_inc)/ (n_b_corr + n_b_inc), 1.);

            double e_for = h_for->GetBinError(i,j);
            double e_back = h_back->GetBinError(i,j);

            afb_mc[j-1] = (n_for - n_back) / (n_for + n_back);
            double afb_mc_st = (n_for_st - n_back_st) / (n_for_st + n_back_st);

            double afb_meas = y_amc[j-1];
            //double afb_meas = y_comb[j-1];
            
            double f_for = (1 + afb_meas)/2.;
            double f_back = (1 - afb_meas)/2.;

            afb_fit[j-1] = dilu_f * f_for - dilu_b * f_back;
            afb_fit_unc[j-1] = ((dilu_f + dilu_b) / 2.) * y_comb_errs[j-1]; //rough estimate (too lazy for full calc)

            ratio[j-1] = afb_fit[j-1] / afb_mc[j-1];
            ratio_unc[j-1] = afb_fit_unc[j-1] / afb_mc[j-1];

            printf("%i %i Raw %.3f, true %.3f, corrected %.3f ratio %.3f  \n", i, j, afb_mc[j-1], afb_mc_st, afb_fit[j-1], ratio[j-1]);



        }
        g_mc.push_back(new TGraphErrors(n_m_bins, m, afb_mc, nullptr, nullptr));
        g_fitted.push_back( new TGraphErrors(n_m_bins, m, afb_fit, m_err, afb_fit_unc));
        g_ratio.push_back(new TGraphErrors(n_m_bins, m, ratio, m_err, ratio_unc));

    }


    for(int i=0; i<n_eta_bins; i++){
        char label[10], title[40];
        sprintf(label, "eta%i", i);
        sprintf(title, " %.1f < |y| < %.1f", eta_bins[i], eta_bins[i+1]);

        draw_graph_ratio(g_fitted[i], g_mc[i], g_ratio[i], label, title);
    }

}
