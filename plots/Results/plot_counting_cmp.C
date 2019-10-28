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
#include "results.h"

void fill_fb_dilu(TTree *t_dy, TH2D *h_for, TH2D *h_back, TH2D *h_corr, TH2D *h_inc, int flag1 = FLAG_MUONS, int year = 2016){
        TempMaker tm(t_dy, false, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;
        tm.is_gen_level = true;
        tm.setup();
        int nEvents=0;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = tm.m > 150. && tm.met_pt < 50.  && tm.has_no_bjets && tm.not_cosmic;

            if(pass){
                nEvents++;
                double mult = tm.cost * tm.cost_st;
                double weight = tm.getEvtWeight();

                Double_t y = abs(tm.cm.Rapidity());

                if(tm.cost > 0.){
                    h_for->Fill(y, tm.m, weight);
                }
                else{
                    h_back->Fill(y, tm.m, weight);
                }

                if(mult > 0.) h_corr->Fill(y, tm.m, weight);
                else h_inc->Fill(y, tm.m, weight);
            }
        }
        printf("selected %i events \n", nEvents);

}


void plot_counting_cmp(){
    int year = 2017;
    int flag = FLAG_MUONS;
    init(year);
    TTree *t_mc;
    if(flag == FLAG_MUONS) t_mc = t_mumu_mc;
    else t_mc = t_elel_mc;


    const int n_eta_bins = 4;
    double eta_bins[] = {0., 1.0, 1.25, 1.5, 2.4};

    TH2D *h_for = new TH2D("h_for", "Foward", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_back = new TH2D("h_back", "Backward", n_eta_bins, eta_bins, n_m_bins, m_bins);

    TH2D *h_corr = new TH2D("h_corr", "Correct sign", n_eta_bins, eta_bins, n_m_bins, m_bins);
    TH2D *h_inc = new TH2D("h_inc", "Incorrect sign", n_eta_bins, eta_bins, n_m_bins, m_bins);

    TGraph *g_mc[n_eta_bins];
    TGraph *g_fitted[n_eta_bins];

    fill_fb_dilu(t_mumu_mc, h_for, h_back, h_corr, h_inc, FLAG_MUONS, year);
    h_for->Print();
    h_back->Print();
    h_corr->Print();
    h_inc->Print();

    Double_t afb_mc[8], afb_fit[8], afb_fit_unc[8];
    for(int i=1; i<= n_eta_bins; i++){
        for(int j=1; j<= n_m_bins; j++){
            double n_for = h_for->GetBinContent(i,j);
            double n_back = h_back->GetBinContent(i,j);
            double n_corr = h_corr->GetBinContent(i,j);
            double n_inc = h_inc->GetBinContent(i,j);

            double e_for = h_for->GetBinError(i,j);
            double e_back = h_back->GetBinError(i,j);
            double e_corr = h_corr->GetBinError(i,j);
            double e_inc = h_inc->GetBinError(i,j);

            double dilu = (n_corr - n_inc)/ (n_corr + n_inc);

            afb_mc[j-1] = (n_for - n_back) / (n_for + n_back);
            afb_fit[j-1] = dilu * y_comb[j-1];
            afb_fit_unc[j-1] = dilu * y_comb_errs[j-1];

        }
        g_mc[i-i] = new TGraphErrors(n_m_bins, m, afb_mc, m_errs, nullptr); 
        g_fitted[i-i] = new TGraphErrors(n_m_bins, m, afb_fit, m_errs, afb_fit_unc); 
    }
    g_mc[0]->Draw();
    g_fitted[0]->Draw("same");
}
