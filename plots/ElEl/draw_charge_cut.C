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
#include "../../utils/root_files.h"
#include "../../utils/PlotUtils.C"




void draw_charge_cut(){
    const int type = FLAG_ELECTRONS;
    int year = 2016;
    char *fin = "../analyze/output_files/2017/ElEl17_dy_oct15.root";


    TFile *f = TFile::Open(fin);

    TTree *T_sig = (TTree *) gDirectory->Get("T_sig");
    TTree *T_ss = (TTree *) gDirectory->Get("T_ss");




    TH1F *sig_m = new TH1F("sig_m", "Opposite sign", 30, 150, 2000);
    TH1F *sig_cut_m = new TH1F("sig_cut_m", "Opposite sign", 30, 150, 2000);

    TH1F *ss_m = new TH1F("ss_m", "Same Sign", 30, 150, 2000);
    TH1F *ss_cut_m = new TH1F("ss_cut_m", "Same Sign", 30, 150, 2000);

    TH1F *dummy= new TH1F("dummy", "Same Sign", 30, 150, 2000);


    bool do_RC = false;
    float m_low = 150.;
    float m_high= 20000.;

    make_m_cost_pt_xf_hist(T_sig, sig_m, dummy, dummy, dummy, dummy, dummy,  false, type,  do_RC, year, m_low, m_high);
    make_m_cost_pt_xf_hist(T_ss, ss_m, dummy, dummy, dummy, dummy, dummy,  false, type,  do_RC, year, m_low, m_high);

    TFile *f_dummy = (TFile*) TFile::Open("dummy.root", "RECREATE");
    TTree *T_sig_cut = T_sig->CopyTree("el1_gc * el2_gc > 0.1");
    TTree *T_ss_cut = T_ss->CopyTree("el1_gc * el2_gc > 0.1");

    make_m_cost_pt_xf_hist(T_sig_cut, sig_cut_m, dummy, dummy, dummy, dummy, dummy,  false, type,  do_RC, year, m_low, m_high);
    make_m_cost_pt_xf_hist(T_ss_cut, ss_cut_m, dummy, dummy, dummy, dummy, dummy,  false, type,  do_RC, year, m_low, m_high);


    sig_m->SetLineColor(kBlue);
    sig_m->SetLineWidth(3);

    ss_m->SetLineColor(kBlue);
    ss_m->SetLineWidth(3);

    sig_cut_m->SetLineColor(kRed);
    sig_cut_m->SetLineWidth(3);

    ss_cut_m->SetLineColor(kRed);
    ss_cut_m->SetLineWidth(3);

    make_ratio_plot("ElEl_sig_charge_cut.pdf", sig_cut_m, "Charge Quality Cut", sig_m, "No Cut", "Cut/No Cut", "os ee M_{ee}", true, true);
    make_ratio_plot("ElEl_ss_charge_cut.pdf", ss_cut_m, "Charge Quality Cut", ss_m, "No Cut", "Cut/No Cut", "ss ee M_{ee}", true, true);

}



