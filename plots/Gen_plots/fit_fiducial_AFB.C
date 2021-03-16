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

void poisson_errs(TH1 *h){
    for(int i=1; i<= h->GetNbinsX(); i++){
        float cont = h->GetBinContent(i);
        h->SetBinError(i, sqrt(cont));
    }
}


double fit_fcn(double *x, double *par){

    double xx = x[0];

    double AFB = par[0];
    double A0 = par[1];

    double alpha = 2.*A0/(2.-A0);
    double norm = 3./4./(2+alpha);

    double rAlpha = alpha * norm;
    double rPl = (norm + AFB);
    double rMn = (norm - AFB);

    double x_alpha = h_alpha->GetBinContent(h_alpha->FindBin(xx));
    double x_pl = h_pl->GetBinContent(h_pl->FindBin(xx));
    double x_mn = h_mn->GetBinContent(h_mn->FindBin(xx));

    return rAlpha * x_alpha + rPl*x_pl + rMn*x_mn;
}


std::tuple<float,float,float,float,float,float> get_afb_fid(float m_low, float m_high, int year, string sys){

    static bool init = false;

    static TFile *f_gen;
    static TTree *t_gen_mu, *t_gen_el;
    int n_bins = 40;
    float bin_size = 2./n_bins;

    static TF1 func, my_func;
    bool do_ptrw = true;

    if(!init){

        f_gen = TFile::Open("../analyze/output_files/DY17_gen_level_nov13.root");
        gROOT->SetBatch(1);

        t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
        t_gen_el = (TTree *) f_gen->Get("T_gen_el");


        char dirname[40];


        h_uncut = new TH1F("cost_uncut", "", n_bins, -1., 1.);
        h_raw = new TH1F("cost_data_obs", "", n_bins, -1., 1.);
        h_sym = new TH1F("cost_sym", "", n_bins, -1., 1.);
        h_asym = new TH1F("cost_asym", "", n_bins, -1., 1.);
        h_pl = new TH1F("cost_pl", "", n_bins, -1., 1.);
        h_mn = new TH1F("cost_mn", "", n_bins, -1., 1.);
        h_alpha = new TH1F("cost_alpha", "", n_bins, -1., 1.);

        func = TF1("func_uncut", "3./8.*(1 + x*x + ([1]/2.)*(1-3*x*x)) + [0]*x", -1., 1.);

        func.SetParameter(0,0.6);
        func.SetParameter(1,0.08);
        func.SetParName(0, "AFB");
        func.SetParName(1, "A0");

        my_func = TF1 ("fit_fcn", fit_fcn, -1., 1., 2);
        my_func.SetParLimits(0, -0.7, 0.7);
        my_func.SetParLimits(1, -0.5, 0.5);
        my_func.SetParameter(0, 0.6);
        my_func.SetParameter(1, 0.08);
        my_func.SetParName(0, "AFB");
        my_func.SetParName(1, "A0");

        init = true;
    }

    int nEvents = 0;


    h_uncut->Reset();
    h_raw->Reset();
    h_pl->Reset();
    h_mn->Reset();
    h_alpha->Reset();
    h_asym->Reset();
    h_sym->Reset();


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

    poisson_errs(h_raw);


    h_uncut->Scale(1./h_uncut->Integral() / bin_size);
    h_uncut->Fit(&func);

    h_raw->Fit(&my_func, "L");




    /*
    TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
    h_uncut->Draw();
    c1->Print("test.png");
    */

    float afb_full = func.GetParameter(0);
    float a0_full = func.GetParameter(1);

    float afb_fid = my_func.GetParameter(0);
    float a0_fid = my_func.GetParameter(1);

    float afb_fid_unc = max(my_func.GetParError(0), func.GetParError(0));
    float a0_fid_unc = max(my_func.GetParError(1), func.GetParError(0));


    return std::make_tuple(afb_full, a0_full, afb_fid, a0_fid, afb_fid_unc, a0_fid_unc);
}









void fit_fiducial_AFB(){

    int year = 2017;

    float m_low, m_high;
    int i_start=1;

    int i_max = n_m_bins;
    //int i_max = 1;

    vector<string>  sys_labels_raw = { "pdf", "RENORM", "FAC", "REFAC","ptcut"};

    vector<string> sys_labels;

    for(auto iter = sys_labels_raw.begin(); iter !=sys_labels_raw.end(); iter++){

          auto cpy = *iter;
          sys_labels.push_back(iter->append("Up"));
          sys_labels.push_back(cpy.append("Down"));
      }
        
    float afb_full, a0_full, afb_fid, a0_fid, afb_fid_unc, a0_fid_unc;
    float afb_shift[8], a0_shift[8], afb_shift_stat_unc[8], afb_shift_sys_unc[8], a0_shift_stat_unc[8], a0_shift_sys_unc[8];
    

    for(int i=i_start; i<i_max; i++){

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);




        string sys = string("");

        afb_shift_sys_unc[i] =0;
        a0_shift_sys_unc[i] =0;

        std::tie(afb_full, a0_full, afb_fid, a0_fid, afb_fid_unc, a0_fid_unc) = get_afb_fid(m_low, m_high, year, sys);

        afb_shift[i] = afb_full - afb_fid;
        a0_shift[i] = a0_full - a0_fid;
        afb_shift_stat_unc[i] = afb_fid_unc;
        a0_shift_stat_unc[i] = a0_fid_unc;


        printf("Nominal: \n");
        printf("AFB shift is: %.3f +/- %.3f \n", afb_shift[i], afb_shift_stat_unc[i]);
        printf("A0 shift is: %.3f +/- %.3f \n", a0_shift[i], a0_shift_stat_unc[i]);

        continue;

        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){

            std::tie(afb_full, a0_full, afb_fid, a0_fid, afb_fid_unc, a0_fid_unc) = get_afb_fid(m_low, m_high, year, *iter);
            float afb_shift_temp = afb_full - afb_fid;
            float a0_shift_temp = a0_full - a0_fid;

            printf("Sys %s: \n", iter->c_str());
            printf("Shift diff is: %.3f, %.3f \n", afb_shift[i] - afb_shift_temp, a0_shift[i] - a0_shift_temp);

            afb_shift_sys_unc[i] += pow(afb_shift[i] - afb_shift_temp, 2); 
            a0_shift_sys_unc[i] += pow(a0_shift[i] - a0_shift_temp, 2); 

        }

        afb_shift_sys_unc[i] = sqrt(afb_shift_sys_unc[i]);
        a0_shift_sys_unc[i] = sqrt(a0_shift_sys_unc[i]);

    }
    for(int i=i_start; i<i_max; i++){

        printf("\n\n Final for bin %i: \n", i);
        printf("AFB shift: %.3f +/- %.3f (stat) +/- %.3f (sys) \n", afb_shift[i], afb_shift_stat_unc[i], afb_shift_sys_unc[i]);
        printf("A0 shift: %.3f +/- %.3f (stat) +/- %.3f (sys) \n", a0_shift[i], a0_shift_stat_unc[i], a0_shift_sys_unc[i]);
    }



    


}
