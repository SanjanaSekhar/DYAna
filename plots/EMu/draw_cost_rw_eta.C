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
#include "../../utils/HistMaker.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"



const int year = 2018;
char *plot_dir = "Misc_plots/emu_cost_reweights/";

void make_emu_cost_rap_cut_hist(TTree *t1, TH1F *h_cost, bool is_data = false, 
        int year=2016, float m_low = 150., float m_high = 999999., float Y_low = 0., float Y_high = 3.0, bool ss = false, bool costrw = false){
    Long64_t size  =  t1->GetEntries();
    TempMaker tm(t1, is_data, year);
    tm.do_emu = true;
    if(!is_data) tm.do_emu_costrw = costrw;

    tm.setup();
    int nEvents=0;

    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        float rap = tm.cm.Rapidity();
        bool pass  = tm.m >= m_low && tm.m <= m_high && tm.met_pt < met_cut && tm.has_no_bjets && abs(rap) > Y_low && abs(rap) < Y_high;
        if(pass){
            nEvents++;
            tm.doCorrections();
            tm.getEvtWeight();

            if(ss) h_cost->Fill(abs(tm.cost), tm.evt_weight);
            else h_cost->Fill(tm.cost, tm.evt_weight);


        }
    }
    printf("Selected %i events \n", nEvents);
    tm.finish();
}



void Fakerate_est_emu_rap_cut(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TH1F *h_cost, int flag1 = FLAG_MUONS, 
        int year=2016, float m_low = 150., float m_high = 10000., float Y_low = 0., float Y_high = 3.0, bool reweight = true){
    FakeRate el_FR, mu_FR;
    //TH2D *FR;
    setup_new_el_fakerate(&el_FR, year);
    setup_new_mu_fakerate(&mu_FR, year);
    //FR.h->Print();
    for (int l=0; l<=2; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        Float_t m, xF, cost, jet1_btag, jet2_btag, gen_weight;
        Float_t jet1_pt, jet2_pt, pu_SF;

        Float_t el_id_SF, el_reco_SF;
        Float_t era1_iso_SF, era1_id_SF;
        Float_t era2_iso_SF, era2_id_SF;
        Float_t evt_fakerate, lep1_fakerate, lep2_fakerate, el1_eta, el1_pt, mu1_eta, mu1_pt;
        TLorentzVector *el = 0;
        TLorentzVector *mu = 0;
        Int_t iso_lep;
        Float_t met_pt, el1_charge, mu1_charge;
        Int_t nJets, no_bjets;
        nJets = 2;
        pu_SF=1;
        t->SetBranchAddress("cost", &cost);
        t->SetBranchAddress("mu1_charge", &mu1_charge);
        t->SetBranchAddress("el1_charge", &el1_charge);
        t->SetBranchAddress("met_pt", &met_pt);
        t->SetBranchAddress("jet2_btag", &jet2_btag);
        t->SetBranchAddress("jet1_btag", &jet1_btag);
        t->SetBranchAddress("jet1_pt", &jet1_pt);
        t->SetBranchAddress("jet2_pt", &jet2_pt);
        //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
        //t1->SetBranchAddress("el_fakerate", &el1_fakerate);
        t->SetBranchAddress("el1_pt", &el1_pt);
        t->SetBranchAddress("mu1_pt", &mu1_pt);
        t->SetBranchAddress("el1_eta", &el1_eta);
        t->SetBranchAddress("mu1_eta", &mu1_eta);
        t->SetBranchAddress("nJets", &nJets);
        t->SetBranchAddress("el", &el);
        t->SetBranchAddress("mu", &mu);
        t->SetBranchAddress("has_nobjets", &no_bjets);

        if(l==0 || l==2 ){
            t->SetBranchAddress("iso_lep", &iso_lep);
        }
        if(l==2){
            t->SetBranchAddress("el_id_SF", &el_id_SF);
            t->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
            t->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t->SetBranchAddress("era2_id_SF", &era2_id_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            if(l==0){
                //iso_lep: 0 = muons, 1 electrons
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                else  lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = lep1_fakerate/(1-lep1_fakerate);
            }
            if(l==1){
                lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                lep2_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                evt_fakerate = -(lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));
            }
            if(l==2){
                Double_t era1_SF = era1_id_SF;
                Double_t era2_SF = era2_id_SF;

                Double_t mu_SF, mu_lumi;
                if(year ==2016){
                    mu_SF = (era1_SF *bcdef_lumi16 + era2_SF * gh_lumi16) / 
                    (bcdef_lumi16 + gh_lumi16);
                    mu_lumi=mu_lumi16;
                }
                if(year==2017){
                    mu_SF = era1_SF;
                    mu_lumi=mu_lumi17;
                }
                if(year==2018){
                    mu_SF = (era1_SF *mu_lumi18_era1 + era2_SF * mu_lumi18_era2) / 
                    (mu_lumi18);
                    mu_lumi=mu_lumi18;
                }

                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF *pu_SF *
                    mu_SF  * 1000. * mu_lumi;
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                else            lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = -(lep1_fakerate * mc_weight)/(1-lep1_fakerate);
            }

            TLorentzVector cm = *el + *mu;
            float m = cm.M();
            float rap = cm.Rapidity();
            bool opp_sign =  ((abs(mu1_charge - el1_charge)) > 0.01);
            bool pass = m>= m_low && m <= m_high && met_pt < met_cut  && no_bjets && opp_sign && abs(rap) > Y_low && abs(rap) < Y_high && 
                ((flag1 == FLAG_MUONS && mu1_pt > 27.) || (flag1 == FLAG_ELECTRONS && el1_pt > 29.));
            if(pass){
                //if(l==3) printf("Evt fr %.2e \n", evt_fakerate);
                h_cost->Fill(cost, evt_fakerate);
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h_cost->Integral());
    }

    if(reweight){
        fakes_costrw_helper h_rw;
        setup_fakes_costrw_helper(&h_rw, year);
        //h_cost->Print("range");
        fakes_cost_reweight(h_cost, h_rw.el_rw);
        //h_cost->Print("range");
    }


    printf("Total fakerate est is %.0f \n", h_cost->Integral());

}


int n_bins = 8;
TH1F *emu_data_cost = new TH1F("emu_data_cost", "Data", n_bins, -1.,1.);
TH1F *emu_bkg_cost = new TH1F("emu_bkg_cost", "DiBoson ttbar and wt", n_bins, -1.,1);
TH1F *emu_dy_cost = new TH1F("emu_dy_cost", "DY", n_bins, -1.,1);
TH1F *emu_QCD_cost = new TH1F("emu_QCD_cost", "QCD", n_bins, -1.,1);


void make_emu_cost_ratio(float m_low, float m_high, float Y_low, float Y_high, bool write_out){

    emu_data_cost->Reset();
    emu_bkg_cost->Reset();
    emu_dy_cost->Reset();
    emu_QCD_cost->Reset();

    bool ss = false;
    make_emu_cost_rap_cut_hist(t_emu_data, emu_data_cost, true,   year, m_low, m_high, Y_low, Y_high, ss);
    make_emu_cost_rap_cut_hist(t_emu_ttbar,  emu_bkg_cost, false,  year, m_low, m_high, Y_low, Y_high, ss);
    make_emu_cost_rap_cut_hist(t_emu_diboson,  emu_bkg_cost, false,  year, m_low, m_high, Y_low, Y_high, ss);
    make_emu_cost_rap_cut_hist(t_emu_wt,  emu_bkg_cost, false,   year, m_low, m_high, Y_low, Y_high, ss);
    make_emu_cost_rap_cut_hist(t_emu_dy,  emu_dy_cost, false,   year, m_low, m_high, Y_low, Y_high, ss);


    bool fakes_reweight = false;
    Fakerate_est_emu_rap_cut(t_emu_WJets, t_emu_QCD, t_emu_WJets_contam,  emu_QCD_cost, FLAG_MUONS, year, m_low, m_high, Y_low, Y_high, fakes_reweight);


    bool scale_error=true;
    float qcd_err = 0.4;
    bool add_err = true;
    if(scale_error){
        setHistError(emu_QCD_cost, qcd_err, add_err);
    }

    symmetrize1d(emu_bkg_cost);
    symmetrize1d(emu_QCD_cost);

    char plot_file[100];
    char h_name[100];
    sprintf(h_name, "emu%i_cost_data_sub", year % 2000);
    TH1F *h_emu_data_sub = (TH1F *) emu_data_cost->Clone(h_name);
    h_emu_data_sub->Add(emu_dy_cost, -1);
    h_emu_data_sub->Add(emu_QCD_cost, -1);

    symmetrize1d(h_emu_data_sub);




    bool logy = false;

    sprintf(h_name, "%i: M %.0f-%0.f GeV, |Y| %.1f - %.1f", year, m_low, m_high, Y_low, Y_high);



    TCanvas *c_emu_plot = make_ratio_plot(string("emu_cost_comparison"), h_emu_data_sub, "EMu Data - Other Backgrounds", emu_bkg_cost, "ttbar + tW + diboson MC Estimate", "ratio", 
            "e#mu cos(#theta)", logy, false, 0.5, 1.5, h_name);

    int rap_bin = 0;
    float rap_low_ = Y_low + 0.01;
    while(rap_bin < n_y_bins){
        if(rap_low_ >= y_bins[rap_bin]) rap_bin ++;
        else break;
    }
    rap_bin--;
    if(write_out) {
        sprintf(plot_file, "%sy%i_emu_m%.0f_to_%.0f_rap%i_cost_rw.png", plot_dir, year - 2000, m_low, m_high, rap_bin);
        c_emu_plot->Print(plot_file);
    }

}


void draw_cost_rw_eta(){


    setTDRStyle();
    init_emu(year);
    init_emu_indv_bkgs(year);
    gROOT->SetBatch(1);



    float m_low, m_high, Y_low, Y_high;



    for(int i=0; i < n_emu_rw_m_bins; i++){

        m_low = emu_rw_m_bins[i];
        m_high =emu_rw_m_bins[i+1];
        Y_low = 0.;
        Y_high = 2.4;
        
        make_emu_cost_ratio(m_low, m_high, Y_low, Y_high, true);
    }


    for(int i=0; i < n_y_bins; i++){

        m_low = 170.;
        m_high =2000.;
        Y_low = y_bins[i];
        Y_high = y_bins[i+1];
        
        make_emu_cost_ratio(m_low, m_high, Y_low, Y_high, true);
    }


}
