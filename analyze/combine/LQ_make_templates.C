
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
//#include"Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
//#include "LQ_make_ss_templates.C"
//#include "make_emu_templates.C"
#include "LQ_TemplateUtils.h"




TH1F *h1_elel_asym, *h1_elel_sym; 
TH1F *h1_mumu_asym, *h1_mumu_sym; 
TH1F *h1_elel_pl, *h1_elel_mn, *h1_elel_alpha, *h1_elel_LQpure_u, *h1_elel_LQint_u,*h1_elel_LQpure_d, *h1_elel_LQint_d,*h1_elel_LQpure_u_vec, *h1_elel_LQint_u_vec,*h1_elel_LQpure_d_vec, *h1_elel_LQint_d_vec, *h1_elel_db, *h1_elel_top,  *h1_elel_tautau, *h1_elel_data, *h1_elel_mc, *h1_elel_qcd, *h1_elel_gam;
TH1F *h1_mumu_pl, *h1_mumu_mn, *h1_mumu_alpha, *h1_mumu_LQpure_u, *h1_mumu_LQint_u,*h1_mumu_LQpure_d, *h1_mumu_LQint_d, *h1_mumu_LQpure_u_vec, *h1_mumu_LQint_u_vec ,*h1_mumu_LQpure_d_vec, *h1_mumu_LQint_d_vec ,*h1_mumu_db, *h1_mumu_top, *h1_mumu_tautau, *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd, *h1_mumu_gam;
int n_1d_bins = n_lq_m_bins*(std::round(std::ceil(n_y_bins/2.) * n_cost_bins + std::floor(n_y_bins/2.) * (n_cost_bins-2)));
//Double_t m_LQ;
bool make_ud = true;

//take m_LQ from command line

void make_data_templates(int year, bool scramble_data = true, bool fake_data = true){

    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }


    char title[100];
    sprintf(title, "ee%i_data_obs", year %2000);

    TH3F* h_elel_data = new TH3F(title, "Data template of (m,x_f, cost_r)",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_elel_data->SetDirectory(0);
    sprintf(title, "mumu%i_data_obs", year %2000);
    TH3F* h_mumu_data = new TH3F(title, "Data template of (m,x_f, cost_r)",
            n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_mumu_data->SetDirectory(0);
    bool ss = false;

    if(!fake_data){
        gen_data_template(t_elel_data, h_elel_data,  year, FLAG_ELECTRONS,  scramble_data, ss, use_xF);
        gen_data_template(t_mumu_data, h_mumu_data,  year,  FLAG_MUONS, scramble_data, ss, use_xF);
    }
    else{
        float Afb = 0.61;
        float A0 = 0.06;
        printf("Making fake data \n");
        bool incl_ss = true;
        bool ss_binning = false;
        float elel_sign_scaling, elel_err, mumu_sign_scaling, mumu_err;
        //printf("making ElEl fakes template \n");
        gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_data, year, 
            FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, "");
        //printf("making MuMu fakes template \n");
        gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_data, year, 
            FLAG_MUONS, incl_ss, ss_binning, use_xF, "");
        //one_mc_template(t_mumu_mc, A0, Afb, h_mumu_data, year, m_low, m_high, FLAG_MUONS, use_xF, "");
        //one_mc_template(t_elel_mc, A0, Afb, h_elel_data, year, m_low, m_high, FLAG_ELECTRONS, use_xF, "");
        TTree *mu_ts[3] = {t_mumu_mc, t_mumu_tautau, t_mumu_gamgam };
        TTree *el_ts[3] = {t_elel_mc, t_elel_tautau, t_elel_gamgam };

        //do processes without emu correction first
        bool emu_costrw = false;
        gen_combined_background_template(3, mu_ts, h_mumu_data, year,  FLAG_MUONS,  ss, use_xF, emu_costrw,  "");
        gen_combined_background_template(3, el_ts, h_elel_data, year,  FLAG_ELECTRONS,  ss, use_xF,  emu_costrw, "");

        TTree *mu_ts_rw[3] = {t_mumu_ttbar, t_mumu_wt, t_mumu_diboson};
        TTree *el_ts_rw[3] = {t_elel_ttbar, t_elel_wt, t_elel_diboson};

        emu_costrw = true;
        gen_combined_background_template(3, mu_ts_rw, h_mumu_data, year,  FLAG_MUONS,  ss, use_xF, emu_costrw,  "");
        gen_combined_background_template(3, el_ts_rw, h_elel_data, year, FLAG_ELECTRONS,  ss, use_xF,  emu_costrw, "");
    }

    h1_elel_data = convert3d(h_elel_data);
    h1_mumu_data = convert3d(h_mumu_data);
    

    printf("Integral of data templates are %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    printf("Made data templates \n");
    delete h_elel_data, h_mumu_data;
}


//changed
void make_qcd_templates(int year, const string &sys_label){

     if(sys_label.empty() || sys_label.find("fakes")){
    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }
    
    char title[100];
    //titles taken care of in conversion

    sprintf(title, "ee%i_qcd%s", year %2000, sys_label.c_str());
    TH3F* h_elel_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_elel_qcd->SetDirectory(0);
    sprintf(title, "mumu%i_qcd%s", year %2000, sys_label.c_str());
    TH3F* h_mumu_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_mumu_qcd->SetDirectory(0);
    bool incl_ss = true;
    bool ss_binning = false;
    float elel_sign_scaling, elel_err, mumu_sign_scaling, mumu_err;
    printf("making ElEl fakes template \n");
    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, 
                FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, sys_label);
    printf("making MuMu fakes template \n");
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, FLAG_MUONS, 
                incl_ss, ss_binning, use_xF, sys_label);
   


    printf("Integral of QCD templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral());

    symmetrize3d(h_mumu_qcd);
    symmetrize3d(h_elel_qcd);

    h1_mumu_qcd = convert3d(h_mumu_qcd);
    h1_elel_qcd = convert3d(h_elel_qcd);

    printf("Made qcd templates \n");
    delete h_elel_qcd, h_mumu_qcd;
    }
}

//changed
void make_mc_templates(int year, Double_t m_LQ, const string &sys_label){
	
    
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos && sys_label.find("emu") == string::npos){
        printf("Doing mu only \n");
        do_mu = true;
        do_el = false;
    }
    if(sys_label.find("el") != string::npos){
        printf("Doing el only \n");
        do_mu = false;
        do_el = true;
    }
    else{
        do_mu = true;
        do_el = true;
    }
    if(sys_label.find("fakes") != string::npos){
        do_el = false;
        do_mu = false;
    }
    bool ss= false;
    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }
    
    char title[100];
    if(do_mu){
        printf("making muon mc templates \n");
        sprintf(title, "mumu%i_sym%s", year %2000, sys_label.c_str());
        auto h_mumu_sym = new TH3F(title, "Symmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
        auto h_mumu_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_alpha->SetDirectory(0);
        sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
        auto h_mumu_asym = new TH3F(title, "Asymmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
        //u LQ
		sprintf(title, "mumu%i_LQpure_u%s", year %2000, sys_label.c_str());
        auto h_mumu_LQpure_u = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQpure_u->SetDirectory(0);
        sprintf(title, "mumu%i_LQint_u%s", year %2000, sys_label.c_str());
        auto h_mumu_LQint_u = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQint_u->SetDirectory(0);
        //d LQ
        sprintf(title, "mumu%i_LQpure_d%s", year %2000, sys_label.c_str());
        auto h_mumu_LQpure_d = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQpure_d->SetDirectory(0);
        sprintf(title, "mumu%i_LQint_d%s", year %2000, sys_label.c_str());
        auto h_mumu_LQint_d = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQint_d->SetDirectory(0);
        //u LQ vec
        sprintf(title, "mumu%i_LQpure_u_vec%s", year %2000, sys_label.c_str());
        auto h_mumu_LQpure_u_vec = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQpure_u_vec->SetDirectory(0);
        sprintf(title, "mumu%i_LQint_u_vec%s", year %2000, sys_label.c_str());
        auto h_mumu_LQint_u_vec = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQint_u_vec->SetDirectory(0);
        //d LQ vec
        sprintf(title, "mumu%i_LQpure_d_vec%s", year %2000, sys_label.c_str());
        auto h_mumu_LQpure_d_vec = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQpure_d_vec->SetDirectory(0);
        sprintf(title, "mumu%i_LQint_d_vec%s", year %2000, sys_label.c_str());
        auto h_mumu_LQint_d_vec = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQint_d_vec->SetDirectory(0);

        sprintf(title, "mumu%i_top%s", year %2000, sys_label.c_str());
        auto h_mumu_top = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_top->SetDirectory(0);
        sprintf(title, "mumu%i_db%s", year %2000, sys_label.c_str());
        auto h_mumu_db = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_db->SetDirectory(0);
        sprintf(title, "mumu%i_tautau%s", year %2000, sys_label.c_str());
        auto h_mumu_tautau = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_tautau->SetDirectory(0);
        sprintf(title, "mumu%i_gam%s", year %2000, sys_label.c_str());
        auto h_mumu_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_gam->SetDirectory(0);

        printf("Making mumu mc \n");
        //gen_mc_template includes m_LQ
        gen_mc_SM_template(t_mumu_mc,  h_mumu_sym, h_mumu_asym, h_mumu_alpha, year, FLAG_MUONS, use_xF, sys_label );
        gen_mc_LQ_template(t_mumu_mc,  h_mumu_LQpure_u, h_mumu_LQint_u,h_mumu_LQpure_d, h_mumu_LQint_d, h_mumu_LQpure_u_vec, h_mumu_LQint_u_vec, h_mumu_LQpure_d_vec, h_mumu_LQint_d_vec, year, m_LQ, FLAG_MUONS, make_ud, use_xF, sys_label );

        TTree *mumu_ts[2] = {t_mumu_ttbar, t_mumu_wt};
        printf("Making mumu back \n");
        bool emu_costrw = true;
        gen_combined_background_template(2, mumu_ts, h_mumu_top, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);


        mumu_ts[0] = t_mumu_diboson;
        gen_combined_background_template(1, mumu_ts, h_mumu_db, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);
        emu_costrw = false;

        mumu_ts[0] = t_mumu_tautau;
        printf("Making mumu tautau \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_tautau, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);

        mumu_ts[0] = t_mumu_gamgam;
        printf("Making mumu gamgam \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_gam, year, FLAG_MUONS,  ss, use_xF, emu_costrw, sys_label);



       // symmetrize3d(h_mumu_gam);
        symmetrize3d(h_mumu_top);
        //symmetrize3d(h_mumu_db);

        h1_mumu_sym = convert3d(h_mumu_sym);
        h1_mumu_asym = convert3d(h_mumu_asym);
        h1_mumu_top = convert3d(h_mumu_top);
        h1_mumu_db = convert3d(h_mumu_db);
        h1_mumu_tautau = convert3d(h_mumu_tautau);
        h1_mumu_gam = convert3d(h_mumu_gam);
        h1_mumu_alpha = convert3d(h_mumu_alpha);
        h1_mumu_LQpure_u = convert3d(h_mumu_LQpure_u);
        h1_mumu_LQint_u = convert3d(h_mumu_LQint_u);
        h1_mumu_LQpure_d = convert3d(h_mumu_LQpure_d);
        h1_mumu_LQint_d = convert3d(h_mumu_LQint_d);
        h1_mumu_LQpure_u_vec = convert3d(h_mumu_LQpure_u_vec);
        h1_mumu_LQint_u_vec = convert3d(h_mumu_LQint_u_vec);
        h1_mumu_LQpure_d_vec = convert3d(h_mumu_LQpure_d_vec);
        h1_mumu_LQint_d_vec = convert3d(h_mumu_LQint_d_vec);
        delete h_mumu_alpha, h_mumu_top, h_mumu_db, h_mumu_tautau, h_mumu_gam;

    }
    if(do_el){
        printf("making electron mc templates \n");
        sprintf(title, "ee%i_sym%s", year %2000, sys_label.c_str());
        auto h_elel_sym = new TH3F(title, "Symmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        sprintf(title, "ee%i_alpha%s", year %2000, sys_label.c_str());
        auto h_elel_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_alpha->SetDirectory(0);
        sprintf(title, "ee%i_asym%s", year %2000, sys_label.c_str());
        auto h_elel_asym = new TH3F(title, "Asymmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);
        // LQ u quark
        sprintf(title, "ee%i_LQpure_u%s", year %2000, sys_label.c_str());
        auto h_elel_LQpure_u = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQpure_u->SetDirectory(0);
        sprintf(title, "ee%i_LQint_u%s", year %2000, sys_label.c_str());
        auto h_elel_LQint_u = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQint_u->SetDirectory(0);
        // LQ d quark
        sprintf(title, "ee%i_LQpure_d%s", year %2000, sys_label.c_str());
        auto h_elel_LQpure_d = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQpure_d->SetDirectory(0);
        sprintf(title, "ee%i_LQint_d%s", year %2000, sys_label.c_str());
        auto h_elel_LQint_d = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQint_d->SetDirectory(0);
        // LQ vec u quark
        sprintf(title, "ee%i_LQpure_u_vec%s", year %2000, sys_label.c_str());
        auto h_elel_LQpure_u_vec = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQpure_u_vec->SetDirectory(0);
        sprintf(title, "ee%i_LQint_u_vec%s", year %2000, sys_label.c_str());
        auto h_elel_LQint_u_vec = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQint_u_vec->SetDirectory(0);
        // LQ vec d quark
        sprintf(title, "ee%i_LQpure_d_vec%s", year %2000, sys_label.c_str());
        auto h_elel_LQpure_d_vec = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQpure_d_vec->SetDirectory(0);
        sprintf(title, "ee%i_LQint_d_vec%s", year %2000, sys_label.c_str());
        auto h_elel_LQint_d_vec = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQint_d_vec->SetDirectory(0);

        sprintf(title, "ee%i_top%s", year %2000, sys_label.c_str());
        auto h_elel_top = new TH3F(title, "Combined background template",
               n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_top->SetDirectory(0);
        sprintf(title, "ee%i_db%s", year %2000, sys_label.c_str());
        auto h_elel_db = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_db->SetDirectory(0);
        sprintf(title, "ee%i_tautau%s", year %2000, sys_label.c_str());
        auto h_elel_tautau = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_tautau->SetDirectory(0);
        sprintf(title, "ee%i_gam%s", year %2000, sys_label.c_str());
        auto h_elel_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_gam->SetDirectory(0);

        printf("starting elel dy \n");
        //gen_mc_template includes m_LQ
        gen_mc_SM_template(t_elel_mc,  h_elel_sym, h_elel_asym, h_elel_alpha, year, FLAG_ELECTRONS, use_xF, sys_label );
        gen_mc_LQ_template(t_elel_mc,  h_elel_LQpure_u, h_elel_LQint_u, h_elel_LQpure_d, h_elel_LQint_d, h_elel_LQpure_u_vec, h_elel_LQint_u_vec, h_elel_LQpure_d_vec, h_elel_LQint_d_vec, year, m_LQ, FLAG_ELECTRONS, make_ud, use_xF, sys_label );

        TTree *elel_ts[2] = {t_elel_ttbar, t_elel_wt};
        bool emu_costrw = true;
        gen_combined_background_template(2, elel_ts, h_elel_top, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);

        elel_ts[0] = t_elel_diboson;
        gen_combined_background_template(1, elel_ts, h_elel_db, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);
        emu_costrw = false;

        elel_ts[0] = t_elel_tautau;
        gen_combined_background_template(1, elel_ts, h_elel_tautau, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);

        elel_ts[0] = t_elel_gamgam;
        printf("Making ee gamgam \n");
        gen_combined_background_template(1, elel_ts, h_elel_gam, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);


       // symmetrize3d(h_elel_gam);
        symmetrize3d(h_elel_top);
        //symmetrize3d(h_elel_db);
        
        h1_elel_sym = convert3d(h_elel_sym);
        h1_elel_asym = convert3d(h_elel_asym);
        h1_elel_top = convert3d(h_elel_top);
        h1_elel_db = convert3d(h_elel_db);
        h1_elel_tautau = convert3d(h_elel_tautau);
        h1_elel_gam = convert3d(h_elel_gam);
        h1_elel_alpha = convert3d(h_elel_alpha);
        h1_elel_LQpure_u = convert3d(h_elel_LQpure_u);
        h1_elel_LQint_u = convert3d(h_elel_LQint_u);
        h1_elel_LQpure_d = convert3d(h_elel_LQpure_d);
        h1_elel_LQint_d = convert3d(h_elel_LQint_d);
        h1_elel_LQpure_u_vec = convert3d(h_elel_LQpure_u_vec);
        h1_elel_LQint_u_vec = convert3d(h_elel_LQint_u_vec);
        h1_elel_LQpure_d_vec = convert3d(h_elel_LQpure_d_vec);
        h1_elel_LQint_d_vec = convert3d(h_elel_LQint_d_vec);
        delete h_elel_alpha, h_elel_top, h_elel_db, h_elel_tautau, h_elel_gam;
    }

}


//changed
void convert_mc_templates(int year, const string &sys_label){
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos && sys_label.find("emu") == string::npos){
         printf("Doing mu only \n");
        do_mu = true;
        do_el = false;
    }
    else if(sys_label.find("el") != string::npos){
        printf("Doing el only \n");
        do_mu = false;
        do_el = true;
    }
    else{
        do_mu = true;
        do_el = true;
    }
    if(sys_label.find("fakes") != string::npos){
        do_el = false;
        do_mu = false;
    }
    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }
        //merge highest rap bin for high mass templates
    //if(m_low > 550.) n_var1_bins -= 1;
    int n_1d_bins = get_n_1d_bins(n_var1_bins, n_cost_bins);
    if(do_mu){

        char title[100];
        sprintf(title, "mumu%i_fpl%s", year%2000, sys_label.c_str());
        h1_mumu_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mumu_pl->SetDirectory(0);
        sprintf(title, "mumu%i_fmn%s", year%2000, sys_label.c_str());
        h1_mumu_mn = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mumu_mn->SetDirectory(0);
        make_pl_mn_templates(h1_mumu_sym, h1_mumu_asym, h1_mumu_pl, h1_mumu_mn);


        h1_mumu_sym->Reset();
        h1_mumu_asym->Reset();


    }

    if(do_el){

        char title[100];
        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_elel_pl->SetDirectory(0);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_elel_mn->SetDirectory(0);
        make_pl_mn_templates(h1_elel_sym, h1_elel_asym, h1_elel_pl, h1_elel_mn);

        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl->SetName(title);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn->SetName(title);

        h1_elel_sym->Reset();
        h1_elel_asym->Reset();

    }
}

void write_out_non_sys_templates(){

    h1_mumu_data->Write();
    h1_elel_data->Write();
  

    h1_mumu_data->Reset();
    h1_elel_data->Reset();
    
}


void write_out_templates(const string &sys_label){

    bool do_mu, do_el;
      bool do_fakes = false;
    if(sys_label.empty()) do_fakes = true;

    if(sys_label.find("mu") != string::npos && sys_label.find("emu") == string::npos){
        do_mu = true;
        do_el = false;
    }
    else if(sys_label.find("el") != string::npos){
        do_mu = false;
        do_el = true;
    }
    else{
        do_mu = true;
        do_el = true;
    }
    if(sys_label.find("fakes") != string::npos){
        do_el = false;
        do_mu = false;
        do_fakes = true;
    }

    if(do_mu){
        h1_mumu_top->Write();
        h1_mumu_db->Write();
        h1_mumu_tautau->Write();
        h1_mumu_gam->Write();
        h1_mumu_alpha->Write();
        h1_mumu_pl->Write();
        h1_mumu_mn->Write();
        h1_mumu_LQpure_u->Write();
        h1_mumu_LQint_u->Write();
        h1_mumu_LQpure_d->Write();
        h1_mumu_LQint_d->Write();
        h1_mumu_LQpure_u_vec->Write();
        h1_mumu_LQint_u_vec->Write();
        h1_mumu_LQpure_d_vec->Write();
        h1_mumu_LQint_d_vec->Write();


        h1_mumu_top->Reset();
        h1_mumu_db->Reset();
        h1_mumu_tautau->Reset();
        h1_mumu_gam->Reset();
        h1_mumu_alpha->Reset();
        h1_mumu_pl->Reset();
        h1_mumu_mn->Reset();
        h1_mumu_LQpure_u->Reset();
        h1_mumu_LQint_u->Reset();
        h1_mumu_LQpure_d->Reset();
        h1_mumu_LQint_d->Reset();
        h1_mumu_LQpure_u_vec->Reset();
        h1_mumu_LQint_u_vec->Reset();
        h1_mumu_LQpure_d_vec->Reset();
        h1_mumu_LQint_d_vec->Reset();
    }

    if(do_el){
        h1_elel_top->Write();
        h1_elel_db->Write();
        h1_elel_tautau->Write();
        h1_elel_gam->Write();
        h1_elel_alpha->Write();
        h1_elel_pl->Write();
        h1_elel_mn->Write();
        h1_elel_LQpure_u->Write();
        h1_elel_LQint_u->Write();
        h1_elel_LQpure_d->Write();
        h1_elel_LQint_d->Write();
        h1_elel_LQpure_u_vec->Write();
        h1_elel_LQint_u_vec->Write();
        h1_elel_LQpure_d_vec->Write();
        h1_elel_LQint_d_vec->Write();


        h1_elel_top->Reset();
        h1_elel_db->Reset();
        h1_elel_tautau->Reset();
        h1_elel_gam->Reset();
        h1_elel_alpha->Reset();
        h1_elel_pl->Reset();
        h1_elel_mn->Reset();
        h1_elel_LQpure_u->Reset();
        h1_elel_LQint_u->Reset();
        h1_elel_LQpure_d->Reset();
        h1_elel_LQint_d->Reset();
        h1_elel_LQpure_u_vec->Reset();
        h1_elel_LQint_u_vec->Reset();
        h1_elel_LQpure_d_vec->Reset();
        h1_elel_LQint_d_vec->Reset();

    }
     if(do_fakes){

        h1_mumu_qcd->Write();
        h1_elel_qcd->Write();
        h1_mumu_qcd->Reset();
        h1_elel_qcd->Reset();
    }


}

void LQ_make_templates(int year = -1, string fout_name_temp = "", int iJob =-1, Double_t m_LQ=1000.){
 
   
   // year =2016;
    if(fout_name_temp == "") fout_name_temp = string("combine/templates/test.root");
    if(year == -1) year = 2017;

    bool scramble_data =false ;
    bool fake_data =true; //use mc instead of data
    use_xF = false;

    

    printf("Initializing files \n");
    init(year);
    //init_emu(year);
    printf("Setting up SFs... ");
    setup_all_SFs(year);
    
    printf("   done \n");

//    vector<string> sys_labels {""};

    //for(int i=1;i<=4;i++){
    
    //m_LQ = 1000.;   
    char templates_name[200];
    sprintf(templates_name,"output_files/LQm%i_nonsys_templates%i.root",int(m_LQ),year%2000);
    string fout_name = string(templates_name,200);
    TFile * fout = TFile::Open(fout_name.c_str(), "RECREATE");

    char dirname[40];
    printf("=========================\n m_LQ = %f, year = %d, fake_data = %d \n=========================\n",m_LQ,year,fake_data );

    
/*
    int i_start=0;
    int i_max = n_lq_m_bins;
    if(iJob >0){ 
        i_start =iJob;
        i_max = iJob +1;
    }
*/

    //for(int i=i_start; i<i_max; i++){
    //for(int i=0; i<1; i++){
        fout->cd();
        snprintf(dirname, 10, "LQ");
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

       // m_low = lq_m_bins[i];
        //m_high = lq_m_bins[i+1];
        printf("\n \n Start making templates ");

        make_data_templates(year, scramble_data, fake_data);
        
        //make_ss_data_templates(year);
        //make_ss_mc_templates(year);
        //make_ss_qcd_templates(year);

        string sys_label = string("");
        make_qcd_templates(year,sys_label);
        make_mc_templates(year, m_LQ, sys_label);
        convert_mc_templates(year, sys_label);

        fout->cd();
        gDirectory->cd(dirname);
        write_out_non_sys_templates();
        //write_out_ss_templates();
        write_out_templates(sys_label);
        //write_groups(year, f_log);
   // }


    fout->Close();
    printf("Templates written to %s \n", fout_name.c_str());
    }

//}

