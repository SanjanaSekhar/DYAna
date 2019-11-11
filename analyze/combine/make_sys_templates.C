#include "make_templates.C"
TTree *t_elel_mc_raw, *t_elel_nosig_raw, *t_elel_back_raw;
TTree *t_mumu_mc_raw, *t_mumu_nosig_raw, *t_mumu_back_raw;

void init_condor_files(int year){
    printf("init year %i  \n", year);
    if(year == 2016){

        f_elel_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/ElEl16_dy_sep12.root");
        t_elel_mc_raw = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig_raw = (TTree *) f_elel_mc->Get("T_DY_back");

        f_elel_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/ElEl16_comb_back_sep14.root");
        t_elel_back_raw = (TTree *) f_elel_back ->Get("T_sig");

        f_elel_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/ElEl16_photInd_sep11.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");


        //-------------------------------------------------------------------------------

        f_mumu_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/MuMu16_dy_sep12.root");
        t_mumu_mc_raw = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig_raw = (TTree *) f_mumu_mc->Get("T_DY_back");

        f_mumu_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/MuMu16_comb_back_sep14.root");
        t_mumu_back_raw = (TTree *) f_mumu_back ->Get("T_sig");

        f_mumu_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2016/MuMu16_photInd_sep11.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");

        return;
    }
    if (year == 2017){


        f_elel_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/ElEl17_dy_sep10.root");
        t_elel_mc_raw = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig_raw = (TTree *) f_elel_mc->Get("T_DY_back");

        f_elel_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/ElEl17_comb_back_sep10.root");
        t_elel_back_raw = (TTree *) f_elel_back ->Get("T_sig");

        f_elel_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/ElEl17_photInd_sep14.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");

        //-------------------------------------------------------------------------------


        f_mumu_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/MuMu17_dy_sep10.root");
        t_mumu_mc_raw = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig_raw = (TTree *) f_mumu_mc->Get("T_DY_back");

        f_mumu_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/MuMu17_comb_back_sep10.root");
        t_mumu_back_raw = (TTree *) f_mumu_back ->Get("T_sig");

        f_mumu_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2017/MuMu17_photInd_sep14.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");

        return;

    }
    if (year ==2018){

        f_elel_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/ElEl18_dy_sep13.root");
        t_elel_mc_raw = (TTree *)f_elel_mc->Get("T_sig");
        t_elel_nosig_raw = (TTree *) f_elel_mc->Get("T_DY_back");

        f_elel_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/ElEl18_comb_back_sep14.root");
        t_elel_back_raw = (TTree *) f_elel_back ->Get("T_sig");

        f_elel_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/ElEl18_photInd_sep11.root");
        t_elel_gamgam = (TTree *)f_elel_gamgam->Get("T_sig");



        //-------------------------------------------------------------------------------


        f_mumu_mc = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/MuMu18_dy_sep13.root");
        t_mumu_mc_raw = (TTree *)f_mumu_mc->Get("T_sig");
        t_mumu_nosig_raw = (TTree *) f_mumu_mc->Get("T_DY_back");

        f_mumu_back = (TFile*) TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/MuMu18_comb_back_sep14.root");
        t_mumu_back_raw = (TTree *) f_mumu_back ->Get("T_sig");

        f_mumu_gamgam = TFile::Open("root://131.225.204.161:1094//store/user/oamram/Condor_inputs/DY_files/2018/MuMu18_photInd_sep11.root");
        t_mumu_gamgam = (TTree *)f_mumu_gamgam->Get("T_sig");

    }
}



void make_sys_templates(int nJobs = 1, int iJob =0, int year = 2016, int type=0){

    //const TString pdf_fout_name("combine/templates/sep19_2016_pdf.root");
    const TString pdf_fout_name("output_files/sep19_sys.root");
    TFile *pdf_fout = TFile::Open(pdf_fout_name, "RECREATE");

    init_condor_files(year);





    printf("Setting up SFs... ");
    setup_all_SFs(year);
    printf("   done \n");


    
    vector<string> sys_labels;
    if (type ==0){
        for(int i =1; i<= 60; i++){
            char name1[20], name2[20];
            sprintf(name1, "_pdf%iUp", i);
            sprintf(name2, "_pdf%iDown", i);
            sys_labels.push_back(string(name1));
            sys_labels.push_back(string(name2));

        }
    }
    else{
      vector<string> sys_labels_raw =  { "_RENORM", "_FAC", "_alphaDen", "_alphaS"  };
      
      
      
      vector<string>  sys_labels_uncorr = 
        {"_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear", "_muRC", "_muHLT", 
          "_muID", "_muISO",  "_elHLT", "_elID", "_elRECO", "_Pu", "_BTAG" };

      string yr_string; 
      if(year == 2016) yr_string = string("16");
      if(year == 2017) yr_string = string("17");
      if(year == 2018) yr_string = string("18");

      for(auto iter = sys_labels_uncorr.begin(); iter !=sys_labels_uncorr.end(); iter++){
          //add year label for uncorrelated systematics
          sys_labels_raw.push_back(iter->append(yr_string));
      }
      for(auto iter = sys_labels_raw.begin(); iter !=sys_labels_raw.end(); iter++){

          auto cpy = *iter;
          sys_labels.push_back(iter->append("Up"));
          sys_labels.push_back(cpy.append("Down"));
      }
    }
    /*
      cout << "sys labels: ";
      for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
          cout << *iter;
      }
      cout<< endl;
      exit(1);
      */


    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;


    for(int i=i_start; i<i_max; i++){

        pdf_fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);
        w = new RooWorkspace("w", "w");

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);

        char cut_str[100];
        printf("Cutting templates \n");
        sprintf(cut_str, "(m > %.0f) && (m < %0.f)", m_low - 15., m_high+15.);
       
        t_elel_mc = t_elel_mc_raw->CopyTree(cut_str);
        t_elel_nosig = t_elel_nosig_raw->CopyTree(cut_str);
        t_elel_back = t_elel_back_raw->CopyTree(cut_str);

        t_mumu_mc = t_mumu_mc_raw->CopyTree(cut_str);
        t_mumu_nosig = t_mumu_nosig_raw->CopyTree(cut_str);
        t_mumu_back = t_mumu_back_raw->CopyTree(cut_str);

        int i_sys =0;
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            if(i_sys % nJobs == iJob){
                printf("Making MC templates for sys %s \n", (*iter).c_str());

                Double_t alpha_denom = amc_alpha[i];


                if(iter->find("alphaDenUp") != string::npos) alpha_denom = amc_alpha[i] + amc_alpha_unc[i];
                if(iter->find("alphaDenDown") != string::npos) alpha_denom = amc_alpha[i] - amc_alpha_unc[i];

                printf("alpha_denom %.2f \n", alpha_denom);

                make_mc_templates(year, alpha_denom, *iter);
                convert_mc_templates(year, *iter);
            }
            i_sys++;
        }
        pdf_fout->cd();
        gDirectory->cd(dirname);
        w->Write();
        cleanup_mc_templates();
    }


    pdf_fout->Close();
    printf("Templates written to %s \n", pdf_fout_name.Data());

}

