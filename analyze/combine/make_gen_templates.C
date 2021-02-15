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
#include "../../utils/TemplateMaker_systematics.C"





void make_gen_templates(){

    int year = 2016;
    bool do_ptrw = true;
    string fout_name = string("combine/templates/y16_gen_temps.root");
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_nov13.root");
    gROOT->SetBatch(1);

    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");

    TFile * fout = TFile::Open(fout_name.c_str(), "RECREATE");

    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;

    string sys = "";
    //int i_max = 1;

    int n_bins = 40;
    TH1F *h_uncut = new TH1F("cost_uncut", "", n_bins, -1., 1.);
    TH1F *h_raw = new TH1F("cost_data_obs", "", n_bins, -1., 1.);
    TH1F *h_sym = new TH1F("cost_sym", "", n_bins, -1., 1.);
    TH1F *h_asym = new TH1F("cost_asym", "", n_bins, -1., 1.);
    TH1F *h_pl = new TH1F("cost_pl", "", n_bins, -1., 1.);
    TH1F *h_mn = new TH1F("cost_mn", "", n_bins, -1., 1.);
    TH1F *h_alpha = new TH1F("cost_alpha", "", n_bins, -1., 1.);

    float m_low, m_high;

    

    for(int i=i_start; i<i_max; i++){

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);


        int nEvents = 0;


        nEvents += make_gen_temps(t_gen_mu, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);
        nEvents += make_gen_temps(t_gen_el, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);

        Double_t dnB, dnF;

        Double_t nB = h_raw->IntegralAndError(1,n_bins/2, dnB);
        Double_t nF = h_raw->IntegralAndError(n_bins/2 + 1,n_bins, dnF);
        Double_t n_tot = nB + nF;
        
        Double_t AFB = ((nF - nB))/((nF+nB));
        Double_t dAFB_v2 = sqrt( pow(dnB * 2. * nF / (n_tot*n_tot),2) + pow(dnF * 2. * nB / (n_tot*n_tot),2));

        printf("Counting AFB %.3f +/- %.3f \n", AFB, dAFB_v2);
        printf("F, B : %.3e %.3e \n", nF, nB);
        h_raw->Print("range");




        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);

        make_pl_mn_templates(h_sym, h_asym, h_pl, h_mn); 

        float scale_ = nEvents / h_raw->Integral();

        h_raw->Scale(scale_);
        h_pl->Scale(scale_);
        h_mn->Scale(scale_);
        h_alpha->Scale(scale_);

        fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        h_uncut->Write();
        h_raw->Write();
        h_pl->Write();
        h_mn->Write();
        h_alpha->Write();

        h_uncut->Reset();
        h_raw->Reset();
        h_pl->Reset();
        h_mn->Reset();
        h_alpha->Reset();
        h_asym->Reset();
        h_sym->Reset();


    }


    fout->Close();
    printf("Templates written to %s \n", fout_name.c_str());

}
