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
#include "../../utils/LQ_TemplateMaker_systematics.C"





void LQ_make_gen_templates(){

    int year = 2016;
    bool do_ptrw = false;
    float m_LQ = 1000.;
    string fout_name = string("combine/templates/LQm1000_gen_templates16_020222.root");
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_aug4.root");
    gROOT->SetBatch(1);

    //TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");

    TFile * fout = TFile::Open(fout_name.c_str(), "RECREATE");

    char dirname[40];

    string sys = "";
    
    TH3F* h_data = new TH3F("h_data", "Data template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_data->SetDirectory(0);
    TH3F* h_sym = new TH3F("h_sym", "sym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_sym->SetDirectory(0);
    TH3F* h_asym = new TH3F("h_asym", "asym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_asym->SetDirectory(0);
    TH3F* h_pl = new TH3F("h_pl", "pl template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_pl->SetDirectory(0);
    TH3F* h_mn = new TH3F("h_mn", "mn template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_mn->SetDirectory(0);
    TH3F* h_alpha = new TH3F("h_alpha", "alpha template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_alpha->SetDirectory(0);
    TH3F* h_LQpure_u = new TH3F("h_LQpure_u", "LQpure_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_LQpure_u->SetDirectory(0);
    TH3F* h_LQpure_d = new TH3F("h_LQpure_d", "LQpure_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_LQpure_d->SetDirectory(0);
    TH3F* h_LQint_u = new TH3F("h_LQint_u", "LQint_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_LQint_u->SetDirectory(0);
    TH3F* h_LQint_d = new TH3F("h_LQint_d", "LQint_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
    h_LQint_d->SetDirectory(0);

    
        printf("\n \n Start making gen level templates for LQ+DY");


        int nEvents = 0;


//nEvents += make_gen_temps(t_gen_mu, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);
        nEvents += make_gen_temps(t_gen_el, h_sym, h_asym, h_alpha, h_LQpure_u, h_LQpure_d, h_LQint_u, h_LQint_d,  m_LQ, do_ptrw, year, sys);



        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);
       // h_LQpure_u->Scale(0.5);
       // h_LQpure_d->Scale(0.5);
       // h_LQint_u->Scale(0.5);
       // h_LQint_d->Scale(0.5);

        make_pl_mn_templates(h_sym, h_asym, h_pl, h_mn); 
        /*
        float scale_ = nEvents / h_raw->Integral();

        h_raw->Scale(scale_);
        h_pl->Scale(scale_);
        h_mn->Scale(scale_);
        h_alpha->Scale(scale_);
        */
        fout->cd();
        snprintf(dirname, 10, "LQ");
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        //h_uncut->Write();
        //h_raw->Write();

        h_pl->Write();
        h_mn->Write();
        h_alpha->Write();
        h_sym->Write();
        h_asym->Write();
        h_LQpure_u->Write();
        h_LQpure_d->Write();
        h_LQint_u->Write();
        h_LQint_d->Write();

        //h_uncut->Reset();
        //h_raw->Reset();
        h_pl->Reset();
        h_mn->Reset();
        h_alpha->Reset();
        h_asym->Reset();
        h_sym->Reset();
        h_LQpure_u->Reset();
        h_LQpure_d->Reset();
        h_LQint_u->Reset();
        h_LQint_d->Reset();

    


    fout->Close();
    printf("Templates written to %s \n", fout_name.c_str());

}
