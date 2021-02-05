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


TH1F *h_uncut, *h_raw, *h_sym, *h_asym, *h_pl, *h_mn, *h_alpha;


double fit_fcn(double *x, double *par){

    double xx = x[0];

    double AFB = par[0];
    double A0 = par[1];

    double alpha = 3./4./(2-A0);
    double norm = 3./4./(2+alpha);

    double rAlpha = alpha * norm;
    double rPl = (norm + AFB);
    double rMn = (norm - AFB);

    double x_alpha = h_alpha->GetBinContent(h_alpha->FindBin(xx));
    double x_pl = h_pl->GetBinContent(h_pl->FindBin(xx));
    double x_mn = h_mn->GetBinContent(h_mn->FindBin(xx));

    return rAlpha * x_alpha + rPl*x_pl + rMn*x_mn;
}






void make_gen_sys_templates(){

    int year = 2016;
    bool do_ptrw = true;
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_nov13.root");
    gROOT->SetBatch(1);

    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;

    string sys = "";
    //int i_max = 1;

    int n_bins = 40;
    h_uncut = new TH1F("cost_uncut", "", n_bins, -1., 1.);
    h_raw = new TH1F("cost_data_obs", "", n_bins, -1., 1.);
    h_sym = new TH1F("cost_sym", "", n_bins, -1., 1.);
    h_asym = new TH1F("cost_asym", "", n_bins, -1., 1.);
    h_pl = new TH1F("cost_pl", "", n_bins, -1., 1.);
    h_mn = new TH1F("cost_mn", "", n_bins, -1., 1.);
    h_alpha = new TH1F("cost_alpha", "", n_bins, -1., 1.);

    float m_low, m_high;

    

    for(int i=i_start; i<i_max; i++){

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);


        int nEvents = 0;


        nEvents += make_gen_temps(t_gen_mu, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);
        nEvents += make_gen_temps(t_gen_el, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);


        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);

        make_pl_mn_templates(h_sym, h_asym, h_pl, h_mn); 

        float scale_ = nEvents / h_raw->Integral();

        h_raw->Scale(scale_);
        h_pl->Scale(scale_);
        h_mn->Scale(scale_);
        h_alpha->Scale(scale_);

        TF1 my_func("fit_fcn", fit_fcn, -1., 1., 2);
        my_func.SetParLimits(0, -0.7, 0.7);
        my_func.SetParLimits(1, -0.5, 0.5);
        my_func.SetParameter(0, 0.6);
        my_func.SetParameter(1, 0.08);

        h_raw->Fit(&my_func, "L");
        my_func.Print();


        h_uncut->Reset();
        h_raw->Reset();
        h_pl->Reset();
        h_mn->Reset();
        h_alpha->Reset();
        h_asym->Reset();
        h_sym->Reset();


    }


}
