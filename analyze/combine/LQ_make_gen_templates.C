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
//#include "../../utils/HistUtils.C"
#include "../../utils/ScaleFactors.C"
#include "../../utils/PlotUtils.C"
//#include "../../utils/LQ_TemplateMaker_systematics.C"
#include "LQ_TemplateUtils.h"





void LQ_make_gen_templates(){

    for(int year=2016;year<=2018;year++){
        bool do_ptrw = false;
        float m_LQ = 1000.;

        char fout_name[200];
        sprintf(fout_name,"combine/templates/LQm%i_gen_templates%i_020222.root",int(m_LQ),year%2000);
	string fout_n = string(fout_name, 200);

        char genfile_name[200];
        sprintf(genfile_name,"../analyze/output_files/DY%i_gen_level_aug4.root",year%2000);
        string genfile_n = string(genfile_name,200);
	TFile *f_gen = TFile::Open(genfile_n.c_str());

        TFile *f_gen_data = TFile::Open("../analyze/root_files/LQ_m1000_test_020122.root");
        gROOT->SetBatch(1);

        //TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
        TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");
        TTree *t_gen_data = (TTree *) f_gen_data->Get("T_lhe");

        //calculate the total gen_weights to scale the data temps 
        float gen_weight, sum_weights;
        t_gen_el->SetBranchAddress("gen_weight", &gen_weight);
        for(int i=0; i<t_gen->GetEntries();i++){
      t_gen->GetEntry(i);
      sum_weights+=gen_weight;
    }

        TFile * fout = TFile::Open(fout_n.c_str(), "RECREATE");

        char dirname[40], title[300];

        string sys = "";
        printf("Starting year %i",year);
        
        sprintf(title, "ee%i_data_obs", year %2000);
        TH3F* h_data = new TH3F(title, "Data template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_data->SetDirectory(0);

        sprintf(title, "ee%i_sym", year %2000);
        TH3F* h_sym = new TH3F(title, "sym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_sym->SetDirectory(0);

        sprintf(title, "ee%i_asym", year %2000);
        TH3F* h_asym = new TH3F(title, "asym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_asym->SetDirectory(0);

        sprintf(title, "ee%i_alpha", year %2000);
        TH3F* h_alpha = new TH3F(title, "alpha template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_alpha->SetDirectory(0);
        sprintf(title, "ee%i_LQpure_u", year %2000);
        TH3F* h_LQpure_u = new TH3F(title, "LQpure_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQpure_u->SetDirectory(0);

        sprintf(title, "ee%i_LQpure_d", year %2000);
        TH3F* h_LQpure_d = new TH3F(title, "LQpure_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQpure_d->SetDirectory(0);

        sprintf(title, "ee%i_LQint_u", year %2000);
        TH3F* h_LQint_u = new TH3F(title, "LQint_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQint_u->SetDirectory(0);

        sprintf(title, "ee%i_LQint_d", year %2000);
        TH3F* h_LQint_d = new TH3F(title, "LQint_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQint_d->SetDirectory(0);

        sprintf(title, "ee%i_raw", year %2000);
        TH3F* h_raw = new TH3F(title, "raw template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_raw->SetDirectory(0);

        TH1F *h1_raw, *h1_data, *h1_pl, *h1_mn, *h1_asym, *h1_sym, *h1_alpha, *h1_LQpure_u, *h1_LQint_u,*h1_LQpure_d, *h1_LQint_d;


        
           // printf("\n \n Start making gen level templates for LQ+DY\n");


        int nEvents = 0;
        int nEvents_data = 0;

    //nEvents += make_gen_temps(t_gen_mu, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);
        nEvents += make_gen_temps(t_gen_el, h_raw, h_sym, h_asym, h_alpha, h_LQpure_u, h_LQpure_d, h_LQint_u, h_LQint_d,  m_LQ, do_ptrw, year, sys);
    	//printf("Finished make_gen_temps, nEvents = %i\n",nEvents);
        nEvents_data += make_gen_data_temps(t_gen_data, h_data, year, sum_weights);

        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);
           // h_LQpure_u->Scale(0.5);
           // h_LQpure_d->Scale(0.5);
           // h_LQint_u->Scale(0.5);
           // h_LQint_d->Scale(0.5);
        printf("Starting convert3d\n");

        h1_data = convert3d(h_data);
        h1_raw = convert3d(h_raw);
        h1_sym = convert3d(h_sym);
        h1_asym = convert3d(h_asym);
        h1_alpha = convert3d(h_alpha);
        h1_LQpure_u = convert3d(h_LQpure_u);
        h1_LQint_u = convert3d(h_LQint_u);
        h1_LQpure_d = convert3d(h_LQpure_d);
        h1_LQint_d = convert3d(h_LQint_d);
        delete h_sym,h_asym,h_alpha,h_LQpure_u,h_LQpure_d,h_LQint_u,h_LQint_d;

        printf("Starting make_pl_mn\n");
            //h1_sym->Print("range");
    	//h1_asym->Print("range");

        // n_y_bins -= 1;
        int n_1d_bins = get_n_1d_bins(n_y_bins, n_cost_bins);
	sprintf(title, "ee%i_fpl", year %2000);
        h1_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_pl->SetDirectory(0);
	sprintf(title, "ee%i_fmn", year %2000);
        h1_mn = new TH1F(title, "Minus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mn->SetDirectory(0);


        make_pl_mn_templates(h1_sym, h1_asym, h1_pl, h1_mn); 
        printf("Finished make_pl_mn\n");
        /*
        float scale_ = nEvents / h1_raw->Integral();

        h1_raw->Scale(scale_);
        h1_pl->Scale(scale_);
        h1_mn->Scale(scale_);
        h1_alpha->Scale(scale_);
        h1_LQpure_u->Scale(scale_);
        h1_LQpure_d->Scale(scale_);
        h1_LQint_u->Scale(scale_);
        h1_LQint_d->Scale(scale_);


        scale_ = nEvents_data / h1_data->Integral();
        h1_data->Scale(scale_);
        */
        fout->cd();
        snprintf(dirname, 10, "LQ");
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

            //h_uncut->Write();
            //h_raw->Write();
        h1_data->Write();
        h1_pl->Write();
        h1_mn->Write();
        h1_alpha->Write();
        h1_sym->Write();
        h1_asym->Write();
        h1_LQpure_u->Write();
        h1_LQpure_d->Write();
        h1_LQint_u->Write();
        h1_LQint_d->Write();

            //h_uncut->Reset();
            //h_raw->Reset();
        h1_data->Reset();
        h1_pl->Reset();
        h1_mn->Reset();
        h1_alpha->Reset();
        h1_asym->Reset();
        h1_sym->Reset();
        h1_LQpure_u->Reset();
        h1_LQpure_d->Reset();
        h1_LQint_u->Reset();
        h1_LQint_d->Reset();

        


        fout->Close();
        printf("Templates written to %s \n", fout_n.c_str());

    }
}
