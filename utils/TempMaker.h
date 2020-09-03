#ifndef DYAFB_TEMPMAKERH
#define DYAFB_TEMPMAKERH
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "ScaleFactors.C"
#include "BTagUtils.C"
#include "HistUtils.C"
#include "bins.h"


Float_t bcdef_lumi16 = 3.119 + 4.035 + 4.270 +  2.615 + 5.809;
Float_t gh_lumi16 =  8.754 + 7.655;
Float_t mu_lumi16 = bcdef_lumi16 + gh_lumi16;
Float_t el_lumi16 = 35.9;

Float_t mu_lumi17 = 42.14;
Float_t el_era1_lumi = 4.823 + 9.664;
Float_t el_era2_lumi = 4.252 + 9.278 + 13.540;
Float_t el_lumi17 = 41.52;

Float_t mu_lumi18_era1 = 9.08;
Float_t mu_lumi18_era2 = 52.23;
Float_t mu_lumi18 = mu_lumi18_era1 + mu_lumi18_era2;
Float_t el_lumi18 = 59.74;

float met_cut = 100;


//Average renorm and fac scale reweights for each mass bin (to remove
//normalization component)

//2016
//
double *h_R_up, *h_R_down, *h_F_up, *h_F_down, *h_RF_up, *h_RF_down;

double h_R_up16[8] = { 0.966, 0.964, 0.963, 0.959, 0.955, 0.951, 0.946, 0.939, }; 
double h_R_down16[8] = { 1.023, 1.025, 1.024, 1.026, 1.029, 1.030, 1.034, 1.040, }; 
double h_F_up16[8] = { 1.049, 1.044, 1.037, 1.029, 1.018, 1.001, 0.987, 0.964, }; 
double h_F_down16[8] = { 0.942, 0.947, 0.955, 0.965, 0.977, 0.998, 1.014, 1.041, }; 
double h_RF_up16[8] = { 1.027, 1.020, 1.011, 0.999, 0.982, 0.963, 0.953, 0.943, }; 
double h_RF_down16[8] = { 0.965, 0.973, 0.982, 0.994, 1.006, 1.016, 1.025, 1.033, };


//2017/8
double h_R_up17[8] = { 1.044, 1.039, 1.033, 1.026, 1.019, 1.009, 0.998, 0.984, }; 
double h_R_down17[8] = { 0.947, 0.954, 0.960, 0.969, 0.979, 0.992, 1.004, 1.020, }; 
double h_F_up17[8] = { 0.976, 0.975, 0.975, 0.973, 0.968, 0.961, 0.956, 0.950, }; 
double h_F_down17[8] = { 1.020, 1.020, 1.018, 1.018, 1.021, 1.025, 1.028, 1.033, }; 
double h_RF_up17[8] = { 1.023, 1.018, 1.012, 1.003, 0.990, 0.973, 0.965, 0.955, }; 
double h_RF_down17[8] = { 0.972, 0.979, 0.983, 0.991, 0.999, 1.008, 1.016, 1.025, }; 



class TempMaker{
    public:
        TempMaker(TTree *t, bool isdata = false, int year = 2016);
        ~TempMaker();
        void setup();
        void setup_systematic(const string &s_label);
        void getEvent(int i);
        void doCorrections();
        float getEvtWeight();
        void fixRFNorm(TH2 *h, int mbin, int year);
        void finish();
        float getReweightingDenom();
        float getLQReweightingDenom(int flag);



        int year = 2016;
        bool is_data = false;
        bool is_gen_level = false;
        bool is_one_iso = false;

        bool do_electrons = false;
        bool do_muons = false;
        bool do_emu = false;
        bool do_RC = true;

        int iso_lep = -1;        
        float el_lumi = 0.;

        string sys_label = string("");
        TTree *t_in;

        Float_t m, xF, cost, gen_weight, reweight_a, reweight_s, reweight_alpha, jet1_btag, jet2_btag, cost_st, gen_cost;
        Float_t evt_weight;
        Float_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
        Float_t era2_HLT_SF, era2_iso_SF, era2_id_SF;
        Float_t el_id_SF, el_reco_SF, el_HLT_SF, pu_SF, pu_SF_up, pu_SF_down;
        Float_t jet1_pt, jet2_pt, jet1_eta, jet2_eta;
        Float_t mu1_pt, mu1_eta, mu2_pt, mu2_eta;
        Float_t el1_pt, el1_eta, el2_pt, el2_eta;
        Float_t mu_R_up, mu_R_down, mu_F_up, mu_F_down, 
                 mu_RF_up, mu_RF_down, pdf_up, pdf_down;
        Float_t mu_p_SF, mu_m_SF, mu_p_SF_up, mu_m_SF_up, mu_p_SF_down, mu_m_SF_down, alphaS_up, alphaS_down;
        Float_t cost_pt, met_pt, el1_charge, el2_charge, mu1_charge, mu2_charge;
        Float_t prefire_SF, prefire_SF_up, prefire_SF_down;
        Float_t dummy;
        Float_t jet1_btag_SF = 1.0;
        Float_t jet2_btag_SF = 1.0;

        Float_t evt_pdfweight;
        Float_t pdf_weights[60];
        TLorentzVector *lep_p=0;
        TLorentzVector *lep_m=0;
        TLorentzVector *gen_lep_p=0;
        TLorentzVector *gen_lep_m=0;
        TLorentzVector cm, gen_cm;
        TLorentzVector *mu = 0;
        TLorentzVector *el = 0;
        Float_t pt, gen_m, gen_pt;
        Int_t nJets, pu_NtrueInt, jet1_flavour, jet2_flavour;

        int inc_id1, inc_id2;



        Int_t has_no_bjets = 1;
        bool not_cosmic = true;

        int count = 0;

        //reweight MC to match data dilepton pt distribution
        bool do_ptrw = false;

        //reweight MC backgrounds based on correction derived from data emu
        bool do_emu_costrw = false;

        Long64_t nEntries;


        //systematics flags

        Float_t one = 1.0;
        Float_t *systematic = &one;
        int sys_shift = 0;

        // do_sys: 0 = nominal, 1= var up, -1 = var down
        int do_pdf_sys = 0;
        int do_btag_sys = 0;
        int do_pileup_sys = 0;
        int do_A0_sys = 0;
        int do_ptrw_sys = 0;
        int do_emu_costrw_sys = 0;

        int do_muHLT_barrel_sys = 0;
        int do_muID_barrel_sys = 0;
        int do_muISO_barrel_sys = 0;

        int do_muHLT_endcap_sys = 0;
        int do_muID_endcap_sys = 0;
        int do_muISO_endcap_sys = 0;
        int do_muRC_sys = 0;

        int do_elID_barrel_sys = 0;
        int do_elHLT_barrel_sys = 0;
        int do_elRECO_barrel_sys = 0;
        int do_elID_endcap_sys = 0;
        int do_elHLT_endcap_sys = 0;
        int do_elRECO_endcap_sys = 0;
        int el_SF_pt_range = 0;

        int do_prefire_sys = 0;

        int do_elScale_sys = 0;
        int do_elSmear_sys = 0;
        
        float elp_rescale, elm_rescale;
};


//global systematics stuff so not tied to single instance of class
el_SFs el_SF;
mu_SFs era1_SFs, era2_SFs;
pileup_systematics pu_sys;
ptrw_helper ptrw_SFs; 
emu_costrw_helper emu_costrw;
A0_helpers A0_helper; 
LQ_rw_helper LQ_helper;

#ifndef STAND_ALONE
BTag_readers b_reader;
BTag_effs btag_effs;
#endif

void setup_all_SFs(int year){
    printf("Setting up SF's \n");
#ifndef STAND_ALONE
    setup_btag_SFs(&b_reader, &btag_effs, year);
#endif
    setup_el_SF(&el_SF, year);
    setup_mu_SFs(&era1_SFs, &era2_SFs,  year);
    setup_pileup_systematic(&pu_sys, year); 
    setup_A0_helper(&A0_helper, year);
    setup_LQ_rw_helper(&LQ_helper, year);
    setup_ptrw_helper(&ptrw_SFs, year);
    setup_emu_costrw_helper(&emu_costrw, year);
}

#endif
