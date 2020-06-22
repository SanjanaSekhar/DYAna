#ifndef ROOT_FILES
#define ROOT_FILES

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "bins.h"

TFile *f_elel_mc, *f_elel_data, *f_elel_QCD, *f_elel_WJets, *f_elel_WJets_contam, *f_elel_QCD_contam;
TTree *t_elel_mc, *t_elel_data, *t_elel_QCD, *t_elel_WJets, *t_elel_WJets_contam, *t_elel_QCD_contam, *t_elel_nosig, *t_elel_tautau;

TFile *f_mumu_mc, *f_mumu_data, *f_mumu_QCD, *f_mumu_WJets, *f_mumu_WJets_contam, *f_mumu_QCD_contam;
TTree *t_mumu_mc, *t_mumu_data, *t_mumu_QCD, *t_mumu_WJets, *t_mumu_WJets_contam, *t_mumu_QCD_contam, *t_mumu_nosig, *t_mumu_tautau;

TFile *f_emu_dy, *f_emu_data, *f_emu_QCD, *f_emu_WJets, *f_emu_WJets_contam;
TTree *t_emu_dy, *t_emu_data, *t_emu_QCD, *t_emu_WJets, *t_emu_WJets_contam;

TTree *t_elel_ss_dy,  *t_elel_ss_data; 
TTree *t_mumu_ss_dy,  *t_mumu_ss_data ;

TFile *f_emu_ss_data, *f_emu_ss_ttbar, *f_emu_ss_dy, *f_emu_ss_diboson;
TTree *t_emu_ss_data, *t_emu_ss_ttbar, *t_emu_ss_dy, *t_emu_ss_diboson;

TFile *f_mumu_ttbar, *f_mumu_wt, *f_mumu_diboson;
TTree *t_mumu_ttbar, *t_mumu_wt, *t_mumu_diboson;
TTree *t_mumu_ss_ttbar, *t_mumu_ss_wt, *t_mumu_ss_diboson;

TFile *f_emu_ttbar, *f_emu_wt, *f_emu_diboson;
TTree *t_emu_ttbar, *t_emu_wt, *t_emu_diboson;

TFile *f_elel_ttbar, *f_elel_wt, *f_elel_diboson;
TTree *t_elel_ttbar, *t_elel_wt, *t_elel_diboson;
TTree *t_elel_ss_ttbar, *t_elel_ss_wt, *t_elel_ss_diboson;

TFile *f_mumu_gamgam, *f_elel_gamgam;
TTree *t_mumu_gamgam, *t_elel_gamgam;

void init_indv_bkgs(int year){
    if(year == 2016){
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_ttbar_april9.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");
        t_mumu_ss_ttbar = (TTree *) f_mumu_ttbar ->Get("T_ss");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_wt_april9.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");
        t_mumu_ss_wt = (TTree *) f_mumu_wt ->Get("T_ss");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_diboson_june11.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");
        t_mumu_ss_diboson = (TTree *) f_mumu_diboson ->Get("T_ss");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_ttbar_april9.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");
        t_elel_ss_ttbar = (TTree *) f_elel_ttbar ->Get("T_ss");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_wt_april9.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");
        t_elel_ss_wt = (TTree *) f_elel_wt ->Get("T_ss");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_diboson_june11.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
        t_elel_ss_diboson = (TTree *) f_elel_diboson ->Get("T_ss");
    }
    if(year == 2017){ //---------------------------------------------------------------------------------------
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_ttbar_april9.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");
        t_mumu_ss_ttbar = (TTree *) f_mumu_ttbar ->Get("T_ss");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_wt_april9.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");
        t_mumu_ss_wt = (TTree *) f_mumu_wt ->Get("T_ss");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_diboson_june11.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");
        t_mumu_ss_diboson = (TTree *) f_mumu_diboson ->Get("T_ss");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_ttbar_april9.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");
        t_elel_ss_ttbar = (TTree *) f_elel_ttbar ->Get("T_ss");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_wt_april9.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");
        t_elel_ss_wt = (TTree *) f_elel_wt ->Get("T_ss");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_diboson_june11.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
        t_elel_ss_diboson = (TTree *) f_elel_diboson ->Get("T_ss");
    }
    if(year == 2018){ //---------------------------------------------------------------------------------------
        f_mumu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_ttbar_april9.root");
        t_mumu_ttbar = (TTree *) f_mumu_ttbar ->Get("T_sig");
        t_mumu_ss_ttbar = (TTree *) f_mumu_ttbar ->Get("T_ss");

        f_mumu_wt = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_wt_april9.root");
        t_mumu_wt = (TTree *) f_mumu_wt ->Get("T_sig");
        t_mumu_ss_wt = (TTree *) f_mumu_wt ->Get("T_ss");

        f_mumu_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_diboson_june11.root");
        t_mumu_diboson = (TTree *) f_mumu_diboson ->Get("T_sig");
        t_mumu_ss_diboson = (TTree *) f_mumu_diboson ->Get("T_ss");


        // ---------------------------------------------------------------------------------------------

        f_elel_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_ttbar_april9.root");
        t_elel_ttbar = (TTree *) f_elel_ttbar ->Get("T_sig");
        t_elel_ss_ttbar = (TTree *) f_elel_ttbar ->Get("T_ss");

        f_elel_wt = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_wt_april9.root");
        t_elel_wt = (TTree *) f_elel_wt ->Get("T_sig");
        t_elel_ss_wt = (TTree *) f_elel_wt ->Get("T_ss");

        f_elel_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_diboson_june11.root");
        t_elel_diboson = (TTree *) f_elel_diboson ->Get("T_sig");
        t_elel_ss_diboson = (TTree *) f_elel_diboson ->Get("T_ss");
    }
}


void init_mc(int year){
    init_indv_bkgs(year);

    if(year == 2016){
        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2016/ElEl16_dy_april17.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_tautau= (TTree *) f_elel_mc->Get("T_tautau");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");


        f_elel_gamgam = TFile::Open("../analyze/output_files/2016/ElEl16_phot_ind_april9.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");

//--------------------------------------------------
        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2016/MuMu16_dy_april17.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_tautau = (TTree *) f_mumu_mc->Get("T_tautau");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");

        f_mumu_gamgam = TFile::Open("../analyze/output_files/2016/MuMu16_phot_ind_april9.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");
    }
    else if(year == 2017){

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2017/ElEl17_dy_april17.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_tautau= (TTree *) f_elel_mc->Get("T_tautau");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2017/ElEl17_phot_ind_april9.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");

    //------------------------------------------------------------------------------
    
        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2017/MuMu17_dy_april17.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_tautau = (TTree *) f_mumu_mc->Get("T_tautau");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");

        f_mumu_gamgam = TFile::Open("../analyze/output_files/2017/MuMu17_phot_ind_april9.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");


    }
    else if(year == 2018){

        f_elel_mc = (TFile*) TFile::Open("../analyze/output_files/2018/ElEl18_dy_april17.root");
        t_elel_mc = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig = (TTree *) f_elel_mc->Get("T_DY_back");
        t_elel_tautau= (TTree *) f_elel_mc->Get("T_tautau");
        t_elel_ss_dy = (TTree *)f_elel_mc->Get("T_ss");

        f_elel_gamgam = TFile::Open("../analyze/output_files/2018/ElEl18_phot_ind_april9.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");

        //-----------------------------------------------------------------------------------------------

        f_mumu_mc = (TFile*) TFile::Open("../analyze/output_files/2018/MuMu18_dy_april17.root");
        t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig = (TTree *) f_mumu_mc->Get("T_DY_back");
        t_mumu_tautau = (TTree *) f_mumu_mc->Get("T_tautau");
        t_mumu_ss_dy = (TTree *)f_mumu_mc->Get("T_ss");


        f_mumu_gamgam = TFile::Open("../analyze/output_files/2018/MuMu18_phot_ind_april9.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");

    }
}
    



void init(int year){
    //MC templates
    printf("init year %i  \n", year);
    init_mc(year);
    if(year == 2016){
        f_elel_data = TFile::Open("../analyze/output_files/2016/ElEl16_data_april2.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");



        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2016/ElEl16_fakes_contam_june11.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2016/MuMu16_data_april2.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");



        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2016/MuMu16_fakes_contam_june11.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
        return;
    }
    if (year == 2017){

        f_elel_data = TFile::Open("../analyze/output_files/2017/ElEl17_data_april2.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");



        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2017/ElEl17_fakes_contam_june11.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2017/MuMu17_data_april2.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");

        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2017/MuMu17_fakes_contam_june11.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
        return;

    }
    if (year ==2018){

        f_elel_data = TFile::Open("../analyze/output_files/2018/ElEl18_data_april2.root");
        t_elel_data = (TTree *)f_elel_data->Get("T_sig"); 
        t_elel_ss_data = (TTree *)f_elel_data->Get("T_ss");
        t_elel_WJets = (TTree *) f_elel_data->Get("T_WJets");
        t_elel_QCD = (TTree *) f_elel_data->Get("T_QCD");



        f_elel_WJets_contam = TFile::Open("../analyze/output_files/2018/ElEl18_fakes_contam_june11.root");
        t_elel_WJets_contam = (TTree *)f_elel_WJets_contam->Get("T_WJets");
        t_elel_QCD_contam = (TTree *)f_elel_WJets_contam->Get("T_QCD");



        //-------------------------------------------------------------------------------

        f_mumu_data = TFile::Open("../analyze/output_files/2018/MuMu18_data_april2.root");
        t_mumu_data = (TTree *)f_mumu_data->Get("T_sig"); 
        t_mumu_ss_data = (TTree *)f_mumu_data->Get("T_ss");
        t_mumu_WJets = (TTree *) f_mumu_data->Get("T_WJets");
        t_mumu_QCD = (TTree *) f_mumu_data->Get("T_QCD");


        f_mumu_WJets_contam = TFile::Open("../analyze/output_files/2018/MuMu18_fakes_contam_june11.root");
        t_mumu_WJets_contam = (TTree *)f_mumu_WJets_contam->Get("T_WJets");
        t_mumu_QCD_contam = (TTree *)f_mumu_WJets_contam->Get("T_QCD");
    }


}


void init_emu(int year){


    if(year == 2016){
        f_emu_data = TFile::Open("../analyze/output_files/2016/EMu16_data_april2.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");


        f_emu_dy = TFile::Open("../analyze/output_files/2016/EMu16_dy_april2.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2016/EMu16_fakes_contam_june11.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }
    if (year == 2017){
        f_emu_data = TFile::Open("../analyze/output_files/2017/EMu17_data_april2.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");

                             
        f_emu_dy = TFile::Open("../analyze/output_files/2017/EMu17_dy_april2.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2017/EMu17_fakes_contam_june11.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }
    if (year == 2018){
        f_emu_data = TFile::Open("../analyze/output_files/2018/EMu18_data_april2.root");
        t_emu_data = (TTree *)f_emu_data->Get("T_sig");
        t_emu_WJets = (TTree *)f_emu_data->Get("T_WJets");
        t_emu_QCD = (TTree *)f_emu_data->Get("T_QCD");


        f_emu_dy = TFile::Open("../analyze/output_files/2018/EMu18_dy_april2.root");
        t_emu_dy = (TTree *)f_emu_dy->Get("T_sig");

        f_emu_WJets_contam = TFile::Open("../analyze/output_files/2018/EMu18_fakes_contam_june11.root");
        t_emu_WJets_contam = (TTree *)f_emu_WJets_contam->Get("T_WJets");
    }



}




void init_emu_indv_bkgs(int year){
    if(year == 2016){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_ttbar_april2.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_wt_april2.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2016/EMu16_diboson_june11.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
    if(year == 2017){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_ttbar_april2.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_wt_april2.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2017/EMu17_diboson_june11.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
    if(year == 2018){
        f_emu_ttbar = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_ttbar_april2.root");
        t_emu_ttbar = (TTree *) f_emu_ttbar ->Get("T_sig");

        f_emu_wt = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_wt_april2.root");
        t_emu_wt = (TTree *) f_emu_wt ->Get("T_sig");

        f_emu_diboson = (TFile*) TFile::Open("../analyze/output_files/2018/EMu18_diboson_june11.root");
        t_emu_diboson = (TTree *) f_emu_diboson ->Get("T_sig");
    }
}





#endif
