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
#include "../../utils/bins.h"
#include "../../utils/HistUtils.C"
#include "../../utils/PlotUtils.C"

int make_amc_gen_cost(TTree *t_gen, TH2D *h_2d, float m_low, float m_high){
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





    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        //bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 30.;
        if(m >= m_low && m <= m_high && sig_event){
            nEvents++;
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */
            float my_cost = cost_st;
            h_2d->Fill(m, my_cost, gen_weight * 1000.);

        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}

void normalize(TH2D *h){
    for(int i=1; i<=h->GetNbinsX(); i++){
        for(int j=1; j<=h->GetNbinsY(); j++){
            float xw = h->GetXaxis()->GetBinWidth(i);
            float yw = h->GetYaxis()->GetBinWidth(j);
            float content = h->GetBinContent(i,j);
            float err = h->GetBinError(i,j);
            //printf("i,j xw, yw %i %i %.2f %.2f \n", i,j, xw,yw);

            h->SetBinContent(i,j,content/(xw*yw));
            h->SetBinError(i,j,err/(xw*yw));

        }
    }
}




void LQ_rw_denom(){

    bool write_out = true;
    char *out_file = "../analyze/SFs/2016/LQ_rw_test.root";
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_april17.root");

    TFile * f_out;
    if(write_out)
        f_out = TFile::Open(out_file, "RECREATE");
    

    //TFile *f_gen = TFile::Open("../MuMu17_dy_gen.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/A0_fits/";



  //  int n_cost_bins = 20;
    //float cost_bins[] = {-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.0};
    //float cost_bins[] = {-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};
    //float cost_bins[] = {0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};

    int n_LQ_pt_bins = 1;
    //float LQ_pt_bins[] = {0., 20., 60., 100., 10000};
    float LQ_pt_bins[] = {0., 10000};

   // int n_LQ_m_bins = 37;
   // float m_LQ_bins[] = {350., 375., 400., 425., 450., 475., 500., 525., 550., 575., 600., 625., 650., 675., 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 
   //     1200., 1250., 1300., 1350., 1400., 1450., 1500., 1600., 1700., 1800., 1900.,  2000.,  2500.,   3000., };

    TH2D *h_mu = new TH2D("h_mu", "", n_lq_m_bins, lq_m_bins,   n_cost_bins, cost_bins);
    TH2D *h_el = new TH2D("h_el", "", n_lq_m_bins, lq_m_bins,   n_cost_bins, cost_bins);

    TH1D *h_1d = new TH1D("h1", "", n_cost_bins, cost_bins);


    int nEvents = 0;
    float m_low = lq_m_bins[0];
    float m_high = 14000.;

    make_amc_gen_cost(t_gen_mu,  h_mu, m_low,m_high);
    make_amc_gen_cost(t_gen_el, h_el, m_low, m_high);
    
    /*
    TH1D *h_el_cost_pre = h_el->ProjectionY("h_el_cost_pre", 10,10);
    TH1D *h_mu_cost_pre = h_mu->ProjectionY("h_mu_cost_pre", 10,10);
    make_ratio_plot("LQ_cos_theta_pre.png", h_el_cost_pre, "El",h_mu_cost_pre, "Mu", "El/Mu", "cos(#theta_{*})", false, true);
    */

    normalize(h_mu);
    normalize(h_el);


    /*
    TH1D *h_el_cost = h_el->ProjectionY("h_el_cost", 10,10);
    TH1D *h_mu_cost = h_mu->ProjectionY("h_mu_cost", 10,10);
    make_ratio_plot("LQ_cos_theta.png", h_el_cost, "El",h_mu_cost, "Mu", "El/Mu", "cos(#theta_{*})", false, true);

    TH1D *h_el_m = h_el->ProjectionX("h_el_m");
    TH1D *h_mu_m = h_mu->ProjectionX("h_mu_m");
    make_ratio_plot("LQ_m.png", h_el_m, "El",h_mu_m, "Mu", "El/Mu", "M (GeV)", false, true);


    h_el_cost->SetLineColor(kBlue);
    h_mu_cost->SetLineColor(kRed);
    h_el_cost->Draw();
    h_mu_cost->Draw("same");
    */

    if(write_out){
        f_out->cd();
        h_mu->Write();
        h_el->Write();
        f_out->Print();
        f_out->Close();
    }

}
