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
#include "../../utils/root_files.h"
#include "../../utils/HistUtils.C"
#include "../../utils/ScaleFactors.C"
#include "../../utils/PlotUtils.C"

int make_amc_gen_cost(TTree *t_gen, TH1F *h_cost_st, TH1F *h_cost_r, TH1F *h_pt, TH1F *h_xf,  
        float m_low, float m_high, float pt_low, float pt_high, float rap_low, float rap_high, 
        bool do_ptrw = false, int year = 2016){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    Bool_t sig_event(1);
    t_gen->SetBranchAddress("gen_p", &gen_lep_p);
    t_gen->SetBranchAddress("gen_m", &gen_lep_m);
    //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
    //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
    t_gen->SetBranchAddress("m", &m);
    t_gen->SetBranchAddress("cost", &cost);
    t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);
    t_gen->SetBranchAddress("sig_event", &sig_event);



    ptrw_helper ptrw_SFs; 
    setup_ptrw_helper(&ptrw_SFs, year);


    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        //bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 30.;
        if(m >= m_low && m <= m_high){
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */
            float my_cost = cost_st;
            float rap = abs(cm.Rapidity());
            if(pt >= pt_low && pt <= pt_high && rap >= rap_low && rap <= rap_high){
                if(gen_weight >0) nEvents++;
                else  nEvents--;

                if(do_ptrw){
                    float ptrw = get_ptrw_SF(ptrw_SFs, m, pt, 0); 
                    gen_weight *= ptrw;
                }

                h_cost_st->Fill(my_cost, gen_weight);
                //h_cost_st->Fill(-cost_st, gen_weight);
                h_cost_r->Fill(cost, gen_weight);


                double xf = abs(2.*cm.Pz()/13000.);

                h_pt->Fill(pt, gen_weight);
                h_xf->Fill(xf, gen_weight);
            }

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}




void fit_amc_gen_cost(){

    bool write_out = false;
    int year = 2016;
    bool do_ptrw = true;
    char *out_file = "../analyze/SFs/2017/a0_fits.root";
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_nov13.root");
    gROOT->SetBatch(1);

    TFile * f_out;
    if(write_out)
        f_out = TFile::Open(out_file, "RECREATE");
    

    //TFile *f_gen = TFile::Open("../MuMu17_dy_gen.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/A0_fits/";

    if(do_ptrw) printf("PT reweighting is ON! \n");





    int n_bins = 40;

    TH1F *h_cost1 = new TH1F("h_cost1", "", n_bins, -1., 1.);
    TH1F *h_cost2 = new TH1F("h_cost2", "", n_bins, -1., 1.);
    TH1F *h_fit = new TH1F("h_fit", "", n_bins, -1., 1.);
    TH1F *h_pt = new TH1F("h_pt", "", 40, 0., 400.);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, -1., 1.);
    float bin_size = 2./n_bins;

    int nEvents = 0;
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.47);
    gStyle->SetStatY(0.9);

    //TF1 *func = new TF1("func", "(1 + x*x + [1]*(1-x*x) + (4./3.)*(2. + [1])*[0]*x) /(8./3. + 4.*[1]/3.)", -1., 1.);
    TF1 *func = new TF1("func", "3./8.*(1 + x*x + ([1]/2.)*(1-3*x*x)) + [0]*x", -1., 1.);
    func->SetParameter(0,0.6);
    func->SetParameter(1,0.1);
    func->SetParName(0, "AFB");
    func->SetParName(1, "A0");
    //func->SetParLimits(1, 0., 1.);


    for(int m_idx = 0; m_idx < n_m_bins; m_idx++){
        float m_low = m_bins[m_idx];
        float m_high = m_bins[m_idx+1];
        //float m_low = 150.;
        //float m_high = 200.;
        //float m_mid = 0.5 * (m_low + m_high);
        printf("Mass range from %.0f to %.0f \n", m_low, m_high);
        h_pt->Reset();
        h_cost1->Reset();
        float pt_low = 0.;
        float pt_high = 100000.;
        float rap_low = 0.;
        float rap_high = 10000.;
        char title[100];
        TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
        nEvents = make_amc_gen_cost(t_gen_mu,  h_cost1, h_dummy, h_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, year);
        float norm = h_cost1->Integral();
        nEvents += make_amc_gen_cost(t_gen_el,  h_cost1, h_dummy, h_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, year);
        sprintf(title, "M %.0f-%.0f GeV, dilepton pT; pT", m_low, m_high);
        h_pt->SetTitle(title);
        h_pt->Scale(1./h_pt->Integral());
        h_pt->Draw();
        c1->SetLogy();
        sprintf(title, "%sy%i_m%.0f_pt.png", plot_dir, year -2000, m_low);
        if(write_out) c1->Print(title);

        Double_t dnB, dnF;
        Double_t nB = h_cost1->IntegralAndError(1,n_bins/2, dnB);
        Double_t nF = h_cost1->IntegralAndError(n_bins/2 + 1,n_bins, dnF);
        Double_t n_tot = nB + nF;
        
        Double_t AFB = ((nF - nB))/((nF+nB));
        Double_t dAFB_v2 = sqrt( pow(dnB * 2. * nF / (n_tot*n_tot),2) + pow(dnF * 2. * nB / (n_tot*n_tot),2));



        c1->SetLogy(false);
        h_cost1->Scale(1./h_cost1->Integral() / bin_size);
        h_cost1->Draw();
        h_cost1->Fit(func);
        sprintf(title, "M %.0f-%.0f GeV; cos(#theta_{*})", m_low, m_high);
        h_cost1->SetTitle(title);

        sprintf(title, "%sy%i_m%.0f_cost.png", plot_dir, year-2000,  m_low);
        c1->Print(title);
        //
        printf("Mass range from %.0f to %.0f, nEvents %i \n", m_low, m_high, nEvents);
        printf("AFB: %.4f +/- %.4f \n", func->GetParameter(0), func->GetParError(0));
        printf("A0: %.3f +/- %.3f \n", func->GetParameter(1), func->GetParError(1));
        printf("Counting AFB %.4f +/- %.4f \n", AFB, dAFB_v2);
        printf("Counting: Total NF %.0f NB %.0f \n", nF, nB);

        continue;
        //exit(1);
        //


        for(int pt_idx = 0; pt_idx < n_pt_bins; pt_idx++){
            for(int rap_idx =0; rap_idx < n_rap_bins; rap_idx++){
                pt_low = pt_bins[pt_idx];
                pt_high = pt_bins[pt_idx+1];
                rap_low = rap_bins[rap_idx];
                rap_high = rap_bins[rap_idx+1];
                //pt_low = pt_bins[pt_idx] * m_low;
                //pt_high = pt_bins[pt_idx+1] * m_low;
                printf("pt low, high  %.1f %.1f \n", pt_low, pt_high);
                printf("rap low, high  %.1f %.1f \n", rap_low, rap_high);

                //float pt_low = 0.;
                //float pt_high = 10000.;

                h_cost1->Reset();
                h_pt->Reset();
                nEvents = make_amc_gen_cost(t_gen_mu,  h_cost1, h_dummy, h_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, year);
                nEvents += make_amc_gen_cost(t_gen_el,  h_cost1, h_dummy, h_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, year);



                
                //Normalize hist & divide by bin size
                h_cost1->Scale(1./h_cost1->Integral() / bin_size);
                h_cost1->Fit(func);
                sprintf(title, "a0_y%i_m%i_pt%i_rap%i", year -2000, m_idx, pt_idx, rap_idx);
                RooRealVar *A0_fit = new RooRealVar(title, "", func->GetParameter(1));
                A0_fit->setError(func->GetParError(1));


                sprintf(title, "amc_fit_ratio_y%i_m%i_pt%i_rap%i", year -2000, m_idx, pt_idx, rap_idx);

                TH1F *h_ratio = (TH1F *) h_cost1->Clone(title);
                auto to_remove = h_ratio->GetFunction("func");
                h_ratio->GetListOfFunctions()->Remove(to_remove);

                for(int i=1; i <= n_bins; i++){
                    float bin_cent = h_cost1->GetBinCenter(i);
                    float fit_val = func->Eval(bin_cent);
                    h_fit->SetBinContent(i, fit_val);
                    h_fit->SetBinError(i, 0.);
                }
                h_ratio->Divide(h_fit);


                sprintf(title, "M %.0f-%.0f GeV, pt %.0f-%.0f GeV |Y| %.1f-%.1f; cos(#theta_{*})", m_low, m_high, pt_low, pt_high, rap_low, rap_high);
                h_cost1->SetTitle(title);
                sprintf(title, "%sy%i_m%.0f_ptbin%i_rapbin%i_fit.png", plot_dir, year -2000, m_low, pt_idx, rap_idx);



                TCanvas *c = draw_ratio_plot(string(title), h_cost1, h_ratio, "cos(#theta)", "MC/fit", 0.7, 1.3);
                c->Print(title);
                //h_cost2->Scale(1./h_cost2->Integral() / bin_size);
                //h_cost2->Draw("same");
                //h_cost2->SetLineColor(kBlack);


                printf("AFB: %.4f +/- %.4f \n", func->GetParameter(0), func->GetParError(0));
                printf("A0: %.3f +/- %.3f \n", func->GetParameter(1), func->GetParError(1));


                if(write_out){
                    f_out->cd();
                    h_ratio->Write();
                    A0_fit->Write();
                }
            }

        }
    }
}
