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
#include "../../utils/PlotUtils.C"
#include "../../utils/ScaleFactors.C"

int make_amc_gen_cost(TTree *t_gen, TH1F *h_cost_st, 
        float m_low, float m_high, float pt_low, float pt_high, int year){
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


    A0_helpers A0_helper_;
    setup_A0_helper(&A0_helper_, year);


    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        if(m >= m_low && m <= m_high && sig_event){
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            float my_cost = cost_st;
            if(pt >= pt_low && pt <= pt_high){
                if(gen_weight >0) nEvents++;
                else  nEvents--;

                float denom = get_reweighting_denom(A0_helper_, my_cost, m, pt);
                //float denom = 3./8.*(1.+my_cost*my_cost + 0.5 * 0.05* (1. - 3. *my_cost*my_cost));
                float reweight = 1./denom;
                //printf("%.3f \n", reweight * gen_weight);
                h_cost_st->Fill(my_cost, reweight *gen_weight);
                h_cost_st->Fill(-my_cost, reweight *gen_weight);

            }

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}




void test_rw(){

    int year = 2018;
    TFile *f_gen = TFile::Open("../analyze/output_files/DY18_gen_level_april17.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");

    float A0_fit[n_m_bins][n_pt_bins];
    float A0_fit_unc[n_m_bins][n_pt_bins];

    char plot_dir[] = "Misc_plots/A0_fits/";



    int n_bins = 40;

    TH1F *h_cost1 = new TH1F("h_cost1", "", n_bins, -1., 1.);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, -1., 1.);
    float bin_size = 2./n_bins;

    int nEvents = 0;
    //gROOT->SetBatch(1);
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.47);
    gStyle->SetStatY(0.9);


    float m_low = 150.;
    float m_high = 200.;
    //float m_low = 150.;
    //float m_high = 200.;
    //float m_mid = 0.5 * (m_low + m_high);
    printf("Mass range from %.0f to %.0f \n", m_low, m_high);
    h_cost1->Reset();
    float pt_low = 0.;
    float pt_high = 30.;
    char title[100];
    TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
    nEvents = make_amc_gen_cost(t_gen_mu,  h_cost1, m_low, m_high, pt_low, pt_high, year);
    nEvents += make_amc_gen_cost(t_gen_el,  h_cost1,m_low, m_high, pt_low, pt_high, year);

    h_cost1->Scale(1./h_cost1->Integral());

    h_cost1->Draw();

}
