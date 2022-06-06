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


Float_t bcdef_lumi16 = 5.819 + 2.617 + 4.276 + 4.066 + 3.135;
Float_t gh_lumi16 =  7.642 + 8.723;
Float_t mu_lumi16 = bcdef_lumi16 + gh_lumi16;
Float_t el_lumi16 = mu_lumi16;

Float_t mu_lumi17 = 41.48;
Float_t el_lumi17 = mu_lumi17;

Float_t mu_lumi18_era1 = 8.95;
Float_t mu_lumi18_era2 = 50.79;
Float_t mu_lumi18 = mu_lumi18_era1 + mu_lumi18_era2;
Float_t el_lumi18 = mu_lumi18;

float met_cut = 100;




class TempMaker{
    public:
        TempMaker(TTree *t, bool isdata = false, int year = 2016);
        ~TempMaker();
        void setup();
        void setup_systematic(const string &s_label);
        void getEvent(int i);
        void doCorrections();
        float getEvtWeight(bool add_btag_SF);
        void fixRFNorm(TH2 *h, int mbin, int year);
        float fixRFNorm(int mbin, int year);
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

        int counter =0;

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
        Float_t mu_prefire_SF, mu_prefire_SF_up, mu_prefire_SF_down;
        Float_t top_ptrw;
        Float_t dummy;
        Float_t jet1_btag_SF = 1.0;
        Float_t jet2_btag_SF = 1.0;

        Float_t evt_pdfweight, nnpdf30_weight;
        Float_t pdf_weights[60];
        TLorentzVector *lep_p=0;
        TLorentzVector *lep_m=0;
        TLorentzVector *gen_lep_p=0;
        TLorentzVector *gen_lep_m=0;
        TLorentzVector cm, gen_cm;
        TLorentzVector *mu = 0;
        TLorentzVector *el = 0;
        Float_t pt, gen_m, gen_pt, gen_rap;
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

        int do_muID_SYS_sys = 0;
        int do_muISO_SYS_sys = 0;

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
        int do_elID_SYS_sys = 0;
        int do_elRECO_SYS_sys = 0;
        int el_SF_pt_range = 0;

        int do_prefire_sys = 0;
        int do_mu_prefire_sys = 0;

        int do_elScale_sys = 0;
        int do_elSmear_sys = 0;
        
        float elp_rescale, elm_rescale;

        int btag_mc_eff_idx = 0;
};


//global systematics stuff so not tied to single instance of class
el_SFs el_SF;
mu_SFs era1_SFs, era2_SFs;
pileup_systematics pu_sys;
ptrw_helper ptrw_SFs; 
emu_costrw_helper emu_costrw;
A0_helpers A0_helper; 
LQ_rw_helper LQ_helper;
RF_pdf_norm_helper RF_pdf_helper;

#ifndef STAND_ALONE
BTag_readers b_reader;
BTag_effs btag_effs;
#endif

void setup_all_SFs(int year){
    //printf("Setting up SF's \n");
#ifndef STAND_ALONE
    bool setup_btag_systematics = true;
    setup_btag_SFs(&b_reader, &btag_effs, year, setup_btag_systematics);
#endif
    setup_el_SF(&el_SF, year);
    setup_mu_SFs(&era1_SFs, &era2_SFs,  year);
    setup_pileup_systematic(&pu_sys, year); 
    setup_A0_helper(&A0_helper, year);
    setup_LQ_rw_helper(&LQ_helper, year);
    setup_ptrw_helper(&ptrw_SFs, year);
    setup_emu_costrw_helper(&emu_costrw, year);
    setup_RF_pdf_norm_helper(&RF_pdf_helper, year);
}

#endif
