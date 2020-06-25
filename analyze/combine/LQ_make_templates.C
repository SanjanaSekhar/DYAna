
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
TH1F *h1_elel_pl, *h1_elel_mn, *h1_elel_alpha, *h1_elel_LQpure, *h1_elel_LQint, *h1_elel_back,  *h1_elel_dy_gg, *h1_elel_data, *h1_elel_mc, *h1_elel_qcd, *h1_elel_gam;
TH1F *h1_mumu_pl, *h1_mumu_mn, *h1_mumu_alpha, *h1_mumu_LQpure, *h1_mumu_LQint, *h1_mumu_back,  *h1_mumu_dy_gg, *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd, *h1_mumu_gam;
Double_t m_LQ;

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
        //one_mc_template(t_mumu_mc, A0, Afb, h_mumu_data, year, m_low, m_high, FLAG_MUONS, use_xF, "");
        //one_mc_template(t_elel_mc, A0, Afb, h_elel_data, year, m_low, m_high, FLAG_ELECTRONS, use_xF, "");
        TTree *mu_ts[7] = {t_mumu_mc, t_mumu_nosig, t_mumu_tautau, t_mumu_ttbar, t_mumu_diboson, t_mumu_wt, t_mumu_gamgam};
        TTree *el_ts[7] = {t_elel_mc, t_elel_nosig, t_elel_tautau, t_elel_ttbar, t_elel_diboson, t_elel_wt, t_elel_gamgam};
        //TTree *mu_ts[6] = {t_mumu_nosig, t_mumu_tautau, t_mumu_ttbar, t_mumu_diboson, t_mumu_wt, t_mumu_gamgam};
        //TTree *el_ts[6] = {t_elel_nosig, t_elel_tautau, t_elel_ttbar, t_elel_diboson, t_elel_wt, t_elel_gamgam};
        gen_combined_background_template(7, mu_ts, h_mumu_data, year, FLAG_MUONS,  ss, use_xF,  "");
        gen_combined_background_template(7, el_ts, h_elel_data, year, FLAG_ELECTRONS,  ss, use_xF,  "");
    }

    h1_elel_data = convert3d(h_elel_data);
    h1_mumu_data = convert3d(h_mumu_data);
    

    printf("Integral of data templates are %.2f %.2f \n", h1_elel_data->Integral(), h1_mumu_data->Integral()); 
    printf("Made data templates \n");
    delete h_elel_data, h_mumu_data;
}


//changed
void make_qcd_templates(int year){

    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }
    
    char title[100];
    //titles taken care of in conversion

    sprintf(title, "ee%i_qcd", year %2000);
    TH3F* h_elel_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_elel_qcd->SetDirectory(0);
    sprintf(title, "mumu%i_qcd", year %2000);
    TH3F* h_mumu_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_mumu_qcd->SetDirectory(0);
    bool incl_ss = true;
    bool ss_binning = false;
    float elel_sign_scaling, elel_err, mumu_sign_scaling, mumu_err;
    printf("making ElEl fakes template \n");
    //::tie(elel_sign_scaling, elel_err) = gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, m_low, m_high, 
          //  FLAG_ELECTRONS, incl_ss, ss_binning, use_xF);
    printf("making MuMu fakes template \n");
    incl_ss = false; // muons use os only for their fakes
    //std::tie(mumu_sign_scaling, mumu_err) = gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, m_low, m_high, FLAG_MUONS, 
         //   incl_ss, ss_binning, use_xF);
   


    printf("Integral of QCD templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral());

    symmetrize3d(h_mumu_qcd);
    symmetrize3d(h_elel_qcd);

    h1_mumu_qcd = convert3d(h_mumu_qcd);
    h1_elel_qcd = convert3d(h_elel_qcd);
    //convert_qcd_to_param_hist(h_elel_qcd, f_log, elel_sign_scaling, elel_err, FLAG_ELECTRONS);
    //convert_qcd_to_param_hist(h_mumu_qcd, f_log, mumu_sign_scaling, mumu_err, FLAG_MUONS);

    printf("Made qcd templates \n");
    delete h_elel_qcd, h_mumu_qcd;
}

//changed
void make_mc_templates(int year, const string &sys_label){
	
    
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos){
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
        //LQ pure
		sprintf(title, "mumu%i_LQpure%s", year %2000, sys_label.c_str());
        auto h_mumu_LQpure = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQpure->SetDirectory(0);
        //LQ interference
        sprintf(title, "mumu%i_LQint%s", year %2000, sys_label.c_str());
        auto h_mumu_LQint = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_LQint->SetDirectory(0);

        sprintf(title, "mumu%i_bk%s", year %2000, sys_label.c_str());
        auto h_mumu_back = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_back->SetDirectory(0);
        sprintf(title, "mumu%i_dy_gg%s", year %2000, sys_label.c_str());
        auto h_mumu_dy_gg = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_dy_gg->SetDirectory(0);
        sprintf(title, "mumu%i_gam%s", year %2000, sys_label.c_str());
        auto h_mumu_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_gam->SetDirectory(0);

        printf("Making mumu mc \n");
        //gen_mc_template includes m_LQ
        gen_mc_template(t_mumu_mc,  h_mumu_sym, h_mumu_asym, h_mumu_alpha, h_mumu_LQpure, h_mumu_LQint, year, m_LQ, FLAG_MUONS, use_xF, sys_label );
        TTree *mumu_ts[3] = {t_mumu_ttbar, t_mumu_wt, t_mumu_diboson};
        printf("Making mumu back \n");
        gen_combined_background_template(3, mumu_ts, h_mumu_back, year, FLAG_MUONS,  ss, use_xF,  sys_label);
        mumu_ts[0] = t_mumu_nosig;
        mumu_ts[1] = t_mumu_tautau;
        printf("Making mumu nosig \n");
        gen_combined_background_template(2, mumu_ts, h_mumu_dy_gg, year, FLAG_MUONS,  ss, use_xF,  sys_label);

        mumu_ts[0] = t_mumu_gamgam;
        printf("Making mumu gamgam \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_gam, year, FLAG_MUONS,  ss, use_xF, sys_label);


        symmetrize3d(h_mumu_gam);
        symmetrize3d(h_mumu_back);
        symmetrize3d(h_mumu_dy_gg);

        h1_mumu_sym = convert3d(h_mumu_sym);
        h1_mumu_asym = convert3d(h_mumu_asym);
        h1_mumu_back = convert3d(h_mumu_back);
        h1_mumu_dy_gg = convert3d(h_mumu_dy_gg);
        h1_mumu_gam = convert3d(h_mumu_gam);
        h1_mumu_alpha = convert3d(h_mumu_alpha);
        h1_mumu_LQpure = convert3d(h_mumu_LQpure);
        h1_mumu_LQint = convert3d(h_mumu_LQint);
        delete h_mumu_alpha, h_mumu_back, h_mumu_dy_gg, h_mumu_gam;

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

        sprintf(title, "ee%i_LQpure%s", year %2000, sys_label.c_str());
        auto h_elel_LQpure = new TH3F(title, "LQpure template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQpure->SetDirectory(0);

        sprintf(title, "ee%i_LQint%s", year %2000, sys_label.c_str());
        auto h_elel_LQint = new TH3F(title, "LQint template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_LQint->SetDirectory(0);

        sprintf(title, "ee%i_bk%s", year %2000, sys_label.c_str());
        auto h_elel_back = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_back->SetDirectory(0);
        sprintf(title, "ee%i_dy_gg%s", year %2000, sys_label.c_str());
        auto h_elel_dy_gg = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_dy_gg->SetDirectory(0);
        sprintf(title, "ee%i_gam%s", year %2000, sys_label.c_str());
        auto h_elel_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_gam->SetDirectory(0);

        printf("starting elel dy \n");
        //gen_mc_template includes m_LQ
        //int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F *h_asym, TH3F *h_alpha, TH3F *h_LQpure, TH3F *h_LQint,
        //int year, Double_t m_LQ, int flag1 = FLAG_MUONS, bool use_xF = false,
        //const string &sys_label = "" )
        gen_mc_template(t_elel_mc, h_elel_sym, h_elel_asym, h_elel_alpha, h_elel_LQpure, h_elel_LQint, year, m_LQ, FLAG_ELECTRONS,  use_xF, sys_label);
        TTree *elel_ts[3] = {t_elel_ttbar, t_elel_wt, t_elel_diboson};
        gen_combined_background_template(3, elel_ts, h_elel_back, year, FLAG_ELECTRONS, ss, use_xF, sys_label);
        elel_ts[0] = t_elel_nosig;
        elel_ts[1] = t_elel_tautau;
        gen_combined_background_template(2, elel_ts, h_elel_dy_gg, year, FLAG_ELECTRONS, ss, use_xF, sys_label);

        elel_ts[0] = t_elel_gamgam;
        printf("Making ee gamgam \n");
        gen_combined_background_template(1, elel_ts, h_elel_gam, year,FLAG_ELECTRONS, ss, use_xF, sys_label);


        symmetrize3d(h_elel_gam);
        symmetrize3d(h_elel_back);
        symmetrize3d(h_elel_dy_gg);
        
        h1_elel_sym = convert3d(h_elel_sym);
        h1_elel_asym = convert3d(h_elel_asym);
        h1_elel_back = convert3d(h_elel_back);
        h1_elel_dy_gg = convert3d(h_elel_dy_gg);
        h1_elel_gam = convert3d(h_elel_gam);
        h1_elel_alpha = convert3d(h_elel_alpha);
        h1_elel_LQpure = convert3d(h_elel_LQpure);
        h1_elel_LQint = convert3d(h_elel_LQint);
        delete h_elel_alpha, h_elel_back, h_elel_dy_gg, h_elel_gam;
    }

}


//changed
void convert_mc_templates(int year, const string &sys_label){
    bool do_mu, do_el;
    if(sys_label.find("mu") != string::npos){
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
    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }
    if(do_mu){

        char title[100];
        sprintf(title, "mumu%i_fpl%s", year%2000, sys_label.c_str());
        h1_mumu_pl = new TH1F(title, "Plus template of DY", n_var1_bins * n_cost_bins * n_lq_m_bins, 0, n_var1_bins * n_cost_bins * n_lq_m_bins);
        h1_mumu_pl->SetDirectory(0);
        sprintf(title, "mumu%i_fmn%s", year%2000, sys_label.c_str());
        h1_mumu_mn = new TH1F(title, "Plus template of DY", n_var1_bins * n_cost_bins * n_lq_m_bins, 0, n_var1_bins * n_cost_bins * n_lq_m_bins);
        h1_mumu_mn->SetDirectory(0);
        make_pl_mn_templates(h1_mumu_sym, h1_mumu_asym, h1_mumu_pl, h1_mumu_mn);


        h1_mumu_sym->Reset();
        h1_mumu_asym->Reset();


    }

    if(do_el){

        char title[100];
        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl = new TH1F(title, "Plus template of DY", n_var1_bins * n_cost_bins *n_lq_m_bins, 0, n_var1_bins * n_cost_bins * n_lq_m_bins);
        h1_elel_pl->SetDirectory(0);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn = new TH1F(title, "Plus template of DY", n_var1_bins * n_cost_bins * n_lq_m_bins, 0, n_var1_bins * n_cost_bins * n_lq_m_bins);
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
    h1_mumu_qcd->Write();
    h1_elel_qcd->Write();

    h1_mumu_data->Reset();
    h1_elel_data->Reset();
    h1_mumu_qcd->Reset();
    h1_elel_qcd->Reset();
}


void write_out_templates(const string &sys_label){

    bool do_mu, do_el;

    if(sys_label.find("mu") != string::npos){
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


    if(do_mu){
        h1_mumu_back->Write();
        h1_mumu_dy_gg->Write();
        h1_mumu_gam->Write();
        h1_mumu_alpha->Write();
        h1_mumu_pl->Write();
        h1_mumu_mn->Write();
        h1_mumu_LQpure->Write();
        h1_mumu_LQint->Write();


        h1_mumu_back->Reset();
        h1_mumu_dy_gg->Reset();
        h1_mumu_gam->Reset();
        h1_mumu_alpha->Reset();
        h1_mumu_pl->Reset();
        h1_mumu_mn->Reset();
        h1_mumu_LQpure->Reset();
        h1_mumu_LQint->Reset();
    }

    if(do_el){
        h1_elel_back->Write();
        h1_elel_dy_gg->Write();
        h1_elel_gam->Write();
        h1_elel_alpha->Write();
        h1_elel_pl->Write();
        h1_elel_mn->Write();
        h1_elel_LQpure->Write();
        h1_elel_LQint->Write();


        h1_elel_back->Reset();
        h1_elel_dy_gg->Reset();
        h1_elel_gam->Reset();
        h1_elel_alpha->Reset();
        h1_elel_pl->Reset();
        h1_elel_mn->Reset();
        h1_elel_LQpure->Reset();
        h1_elel_LQint->Reset();

    }


}

void LQ_make_templates(int year = 2016, int nJobs = 6, int iJob =-1){

    
    for(int year=2016; year<=2016; year++){
    

    bool scramble_data =false ;
    bool fake_data =true; //use mc instead of data
    use_xF = false;

    

    printf("Initializing files \n");
    init(year);
    //init_emu(year);
    printf("Setting up SFs... ");
    setup_all_SFs(year);
    
    printf("   done \n");

    vector<string> sys_labels {""};

    for(int i=1;i<=1;i++){
    
    Double_t m_LQ = 1000.*i;   
    char templates_name[100];
    sprintf(templates_name,"combine/templates/LQm%i_normed_templates%i.root",int(m_LQ),year%2000);
    const TString fout_name(templates_name);
    TFile * fout = TFile::Open(fout_name, "RECREATE");

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
        make_qcd_templates(year);
        //make_ss_data_templates(year);
        //make_ss_mc_templates(year);
        //make_ss_qcd_templates(year);

        string sys_label = string("");
        make_mc_templates(year, sys_label);
        convert_mc_templates(year, sys_label);

        fout->cd();
        gDirectory->cd(dirname);
        write_out_non_sys_templates();
        //write_out_ss_templates();
        write_out_templates(sys_label);
        //write_groups(year, f_log);
   // }


    fout->Close();
    printf("Templates written to %s \n", fout_name.Data());
    }
}
}
