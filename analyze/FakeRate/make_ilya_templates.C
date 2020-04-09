
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "TH2D.h"
#include "Math/Functor.h"
#include "../../utils/TemplateMaker_systematics.C"
#include "../../utils/root_files.h"

float gen_ilya_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree* t_QCD_contam, TH2F *h, 
        int year,  Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool ss = false){
    h->Sumw2();
    TH2D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    if(flag1 == FLAG_MUONS) setup_new_mu_fakerate(&FR, year);
    else setup_new_el_fakerate(&FR, year);
    //TH2D *FR;
    h_err = (TH2D *) FR.h->Clone("h_err");
    h_err->Reset();
    for (int l=0; l<=3; l++){
        TTree *t;
        bool is_data, is_one_iso;
        if (l==0){
            t = t_WJets;
            is_data =true;
            is_one_iso = true;
        }
        if (l==1){
            t = t_QCD;
            is_data =true;
            is_one_iso = false;
        }
        if (l==2){
            t = t_WJets_contam;
            is_data =false;
            is_one_iso = true;
        }
        if (l==3){
            t = t_QCD_contam;
            is_data =false;
            is_one_iso = false;
        }
        TempMaker tm(t, is_data, year);
        tm.is_one_iso = is_one_iso;
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;
        tm.setup();

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);

            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.not_cosmic;

            bool opp_sign = false;
            if(flag1 == FLAG_MUONS) opp_sign = ((abs(tm.mu1_charge - tm.mu2_charge)) > 0.01);
            else opp_sign = ((abs(tm.el1_charge - tm.el2_charge)) > 0.01);

            if(!ss) pass = pass && opp_sign;
            if(pass){
                double evt_reweight = 0.;

                if(flag1 == FLAG_MUONS){
                    Double_t mu1_fakerate, mu2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = mu1_fakerate/(1-mu1_fakerate);
                    }

                    if(l==1){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                    }
                    if(l==2){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate )/(1-mu1_fakerate);
                    }
                    if(l==3){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                    }
                    if(l==0 && tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f));
                    if(l==0 && tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f));
                }
                else{
                    Double_t el1_fakerate, el2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = el1_fakerate/(1-el1_fakerate);
                    }
                    if(l==1){
                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                    }
                    if(l==2){

                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate )/(1-el1_fakerate);
                    }
                    if(l==3){

                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                    }
                    if(l==0 && tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f));
                    if(l==0 && tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f));
                }
                double tot_evt_weight = 0.;
                if(is_data) tot_evt_weight = evt_reweight; 
                else{
                    tot_evt_weight = evt_reweight * tm.getEvtWeight();
                }

                if(opp_sign) tot_weight_os += tot_evt_weight;
                else tot_weight_ss += tot_evt_weight;

                Double_t y = abs(tm.cm.Rapidity());

                h->Fill(y, tm.m, tot_evt_weight);


            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());

    }
    printf("Performing fakes cleanup (removing neg. bins) \n");
    cleanup_template(h);

    printf("Total Fakerate weight Weight is %.2f \n", h->Integral());
    if(h->Integral() < 0.){
        h->Scale(0.);
        printf("zeroing Fakes template \n");
    }
    float scaling = 1.;
    if(ss) scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
    return scaling;
}


void make_ilya_templates(){

    int year = 2016;

    printf("Initializing files \n");
    init(year);
    //init_emu(year);

    /*
    if(year == 2016){
        f_mumu_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu16_data_oct30.root");
        f_mumu_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu16_fakes_contam_mlow_nov6.root");

        f_elel_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl16_data_mlow_nov1.root");
        f_elel_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl16_fakes_contam_mlow_nov4.root");
    }

    if(year == 2017){
        f_mumu_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu17_data_mlow_oct21.root");
        f_mumu_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu17_fakes_contam_mlow_nov6.root");

        f_elel_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl17_data_mlow_nov1.root");
        f_elel_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl17_fakes_contam_mlow_nov4.root");
    }
    else if(year == 2018){
        f_mumu_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu18_data_mlow_oct21.root");
        f_mumu_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/MuMu18_fakes_contam_mlow_nov6.root");

        f_elel_data = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl18_data_mlow_nov1.root");
        f_elel_WJets_contam = TFile::Open("root://cmseos.fnal.gov//store/user/oamram/Analysis_ntuples/ElEl18_fakes_contam_mlow_nov4.root");
    }


    t_mumu_WJets = (TTree *)f_mumu_data->Get("T_WJets");
    t_mumu_QCD = (TTree *)f_mumu_data->Get("T_QCD");


    t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
    t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");


    t_elel_WJets = (TTree *)f_elel_data->Get("T_WJets");
    t_elel_QCD = (TTree *)f_elel_data->Get("T_QCD");


    t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
    t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");
    */


    float eta_bins[] = {0., 1.0, 1.25, 1.5, 2.4};
    float mass_bins[] = {150, 171, 200, 250, 320, 510, 700, 1000, 14000}; 
    int n_mbins = 8;
    int n_etabins = 4;

    TH2F *h_elel = new TH2F("h_fakes_elel", "Electron pairs fakes est", n_etabins, eta_bins, n_mbins, mass_bins);
    TH2F *h_mumu = new TH2F("h_fakes_mumu", "Muon pairs fakes est", n_etabins, eta_bins, n_mbins, mass_bins);

    TH2F *h_dummy = new TH2F("h_dummy", "dummy", 100, -100., 100., 100, -100, 100);

    float m_low = 150.;
    float m_high = 14000.;
    bool ss = true;

    bool incl_ss = true;
    bool ss_binning = false;
    bool use_xF = false;
    /*
     gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_dummy, year, m_low, m_high, 
            FLAG_ELECTRONS, incl_ss, ss_binning, use_xF);
     h_dummy->Reset();
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_dummy, year, m_low, m_high, FLAG_MUONS, 
            incl_ss, ss_binning, use_xF);
            */
    float elel_sign_scaling = gen_ilya_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel, year, m_low, m_high, FLAG_ELECTRONS, ss);
    float mumu_sign_scaling = gen_ilya_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu, year, m_low, m_high, FLAG_MUONS, ss);

    char fname[80], mu_plot_name[80], el_plot_name[80];
    sprintf(fname, "lepton_pair_fakes_est_%i.root", year);
    sprintf(mu_plot_name, "mu_fakes_est_%i.png", year);
    sprintf(el_plot_name, "el_fakes_est_%i.png", year);
    printf("Writing out to %s \n", fname);
    TFile *fout = new TFile(fname, "RECREATE");
    fout->cd();
    h_elel->Write();
    h_mumu->Write();
    fout->Close();

    TH1D *h_mu_proj = h_mumu->ProjectionY();
    TH1D *h_el_proj = h_elel->ProjectionY();

    TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
    h_mu_proj->Draw();
    c1->Print(mu_plot_name);

    TCanvas *c2 = new TCanvas("c2", "", 1000, 800);
    h_el_proj->Draw();
    c2->Print(el_plot_name);
}

    
