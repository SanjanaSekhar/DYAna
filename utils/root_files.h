#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

TFile *f_elel_mc, *f_elel_back, *f_elel_data, *f_elel_QCD, *f_elel_WJets, *f_elel_WJets_contam, *f_elel_QCD_contam;
TTree *t_elel_mc, *t_elel_back, *t_elel_data, *t_elel_QCD, *t_elel_WJets, *t_elel_WJets_contam, *t_elel_QCD_contam, *t_elel_nosig;

TFile *f_mumu_mc, *f_mumu_back, *f_mumu_data, *f_mumu_QCD, *f_mumu_WJets, *f_mumu_WJets_contam, *f_mumu_QCD_contam;
TTree *t_mumu_mc, *t_mumu_back, *t_mumu_data, *t_mumu_QCD, *t_mumu_WJets, *t_mumu_WJets_contam, *t_mumu_QCD_contam, *t_mumu_nosig;

TFile *f_emu_dy, *f_emu_back, *f_emu_data, *f_emu_QCD, *f_emu_WJets, *f_emu_WJets_contam;
TTree *t_emu_dy, *t_emu_back, *t_emu_data, *t_emu_QCD, *t_emu_WJets, *t_emu_WJets_contam;

TTree *t_elel_ss_dy, *t_elel_ss_back, *t_elel_ss_data; 
TTree *t_mumu_ss_dy, *t_mumu_ss_back, *t_mumu_ss_data ;

TFile *f_emu_ss_data, *f_emu_ss_ttbar, *f_emu_ss_dy, *f_emu_ss_diboson;
TTree *t_emu_ss_data, *t_emu_ss_ttbar, *t_emu_ss_dy, *t_emu_ss_diboson;

TFile *f_mumu_ttbar, *f_mumu_wt, *f_mumu_diboson;
TTree *t_mumu_ttbar, *t_mumu_wt, *t_mumu_diboson;

TFile *f_emu_ttbar, *f_emu_wt, *f_emu_diboson;
TTree *t_emu_ttbar, *t_emu_wt, *t_emu_diboson;

TFile *f_elel_ttbar, *f_elel_wt, *f_elel_diboson;
TTree *t_elel_ttbar, *t_elel_wt, *t_elel_diboson;

TFile *f_mumu_gamgam, *f_elel_gamgam;
TTree *t_mumu_gamgam, *t_elel_gamgam;

int n_xf_bins = 4;
Float_t xf_bins[] = {0., 0.04, 0.07, 0.10, 1.0};
int n_cost_bins = 10;
Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
int n_cost_ss_bins = n_cost_bins/2;
Float_t cost_ss_bins[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0};
int n_m_bins = 8;
Double_t m_bins[] = {150, 171, 200,  250, 320, 510, 700, 1000, 14000};

Double_t amc_alpha[n_m_bins] =        {0.056, 0.056, 0.047, 0.055, 0.042, 0.030, 0.018, 0.012};
Double_t amc_alpha_unc[n_m_bins] =    {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.007, 0.007};


void init(int year){
    //MC templates
    printf("init year %i  \n", year);
    if(year == 2016){
        f_elel_data = TFile::Open("../analyze/output_files/2016/ElEl16_data_nov1.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_dy_nov1.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");

        f_elel_back = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_comb_back_nov4.root");
        t_elel_back = (TTree *) f_elel_back ->Get("T_sig");
        t_elel_ss_back = (TTree *)f_elel_back->Get("T_ss");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2016/ElEl16_photInd_nov4.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");


        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2016/ElEl16_fakes_contam_nov4.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2016/MuMu16_data_oct30.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");

        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_dy_nov4.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");

        f_mumu_back = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_comb_back_oct31.root");
        t_mumu_back = (TTree *) f_mumu_back ->Get("T_sig");
        t_mumu_ss_back = (TTree *)f_mumu_back->Get("T_ss");

        f_mumu_gamgam = TFile::Open("../analyze/output_files/2016/MuMu16_photind_oct31.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");


        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2016/MuMu16_fakes_contam_oct31.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
        return;
    }
    if (year == 2017){

        f_elel_data = TFile::Open("../analyze/output_files/2017/ElEl17_data_nov1.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_dy_nov1.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");

        f_elel_back = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_comb_back_nov4.root");
        t_elel_back = (TTree *) f_elel_back ->Get("T_sig");
        t_elel_ss_back = (TTree *)f_elel_back->Get("T_ss");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2017/ElEl17_photInd_nov4.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");


        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2017/ElEl17_fakes_contam_nov4.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2017/MuMu17_data_oct21.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");

        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_dy_nov4.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");

        f_mumu_back = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_comb_back_oct21.root");
        t_mumu_back = (TTree *) f_mumu_back ->Get("T_sig");
        t_mumu_ss_back = (TTree *)f_mumu_back->Get("T_ss");

        f_mumu_gamgam = TFile::Open("../analyze/output_files/2017/MuMu17_photInd_oct31.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");


        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2017/MuMu17_fakes_contam_oct21.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
        return;

    }
    if (year ==2018){

        f_elel_data = TFile::Open("../analyze/output_files/2018/ElEl18_data_nov1.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_dy_nov1.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");

        f_elel_back = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_comb_back_nov4.root");
        t_elel_back = (TTree *) f_elel_back ->Get("T_sig");
        t_elel_ss_back = (TTree *)f_elel_back->Get("T_ss");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2018/ElEl18_photInd_nov4.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");


        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2018/ElEl18_fakes_contam_nov4.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2018/MuMu18_data_oct21.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");

        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_dy_nov4.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");

        f_mumu_back = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_comb_back_oct21.root");
        t_mumu_back = (TTree *) f_mumu_back ->Get("T_sig");
        t_mumu_ss_back = (TTree *)f_mumu_back->Get("T_ss");

        f_mumu_gamgam = TFile::Open("../analyze/output_files/2018/MuMu18_photInd_oct21.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");


        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2018/MuMu18_fakes_contam_oct28.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
    }


    

}
void init_emu(int year){


    if(year == 2016){
        f_emu_data = TFile::Open("../analyze/output_files/2016/EMu16_data_sep15.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");

                             
        f_emu_back = TFile::Open("../analyze/output_files/2016/EMu16_comb_back_sep15.root");
        t_emu_back = (TTree *)f_emu_back->Get("T_sig");

        f_emu_dy = TFile::Open("../analyze/output_files/2016/EMu16_dy_sep15.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2016/EMu16_fakes_contam_sep15.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }
    if (year == 2017){
        f_emu_data = TFile::Open("../analyze/output_files/2017/EMu17_data_sep15.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");

                             
        f_emu_back = TFile::Open("../analyze/output_files/2017/EMu17_comb_back_sep15.root");
        t_emu_back = (TTree *)f_emu_back->Get("T_sig");

        f_emu_dy = TFile::Open("../analyze/output_files/2017/EMu17_dy_sep15.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2017/EMu17_fakes_contam_sep15.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }
    if (year == 2018){
        f_emu_data = TFile::Open("../analyze/output_files/2018/EMu18_data_sep15.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");

        f_emu_back = TFile::Open("../analyze/output_files/2018/EMu18_comb_back_sep15.root");
        t_emu_back = (TTree *)f_emu_back->Get("T_sig");

        f_emu_dy = TFile::Open("../analyze/output_files/2018/EMu18_dy_sep15.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2018/EMu18_fakes_contam_sep15.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }



}




void init_emu_indv_bkgs(int year){
    if(year == 2016){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_TTbar_sep15.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_WT_sep15.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_diboson_sep15.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
    if(year == 2017){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_TTbar_sep15.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_WT_sep15.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_diboson_sep15.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
    if(year == 2018){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_TTbar_sep15.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_WT_sep15.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_diboson_sep15.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
}


void init_indv_bkgs(int year){
    if(year == 2016){
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_ttbar_nov4.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_WT_oct29.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_diboson_oct29.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_ttbar_nov4.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_wt_nov4.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_diboson_nov4.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
    }
    if(year == 2017){
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_ttbar_nov4.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_WT_oct21.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_diboson_oct21.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_ttbar_nov4.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_wt_nov4.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_diboson_nov4.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
    }
    if(year == 2018){
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_ttbar_nov4.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_WT_oct21.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_diboson_oct21.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_ttbar_nov4.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_wt_nov4.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_diboson_nov4.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
    }
}




