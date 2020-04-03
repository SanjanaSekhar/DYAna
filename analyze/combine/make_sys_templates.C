#include "make_templates.C"



void make_sys_templates(int nJobs = 1, int iJob =0, int year = 2016, int type=0){
    //type0 is pdfs, type1 other sys
    if(nJobs == 0){
        printf("Invalid setting of 0 total jobs! Going to change it to be 1 job (ie this process runs over all systematics) \n");
        nJobs = 1;
    }

    //const TString pdf_fout_name("combine/templates/sep19_2016_pdf.root");
    const TString pdf_fout_name("output_files/sys_test.root");
    TFile *pdf_fout = TFile::Open(pdf_fout_name, "RECREATE");

    init_mc(year);





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

      vector<string> sys_labels_raw =  { "_RENORM", "_FAC", "_REFAC", "_alphaDen", "_alphaS"  };
      //vector<string> sys_labels_raw =  {};
      
      
      
      vector<string>  sys_labels_uncorr = 
        {"_METJER", "_METJEC", "_METHEM", "_prefire", "_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear", "_muRC", "_Pu", "_BTAG",
            "_muHLTBAR", "_muIDBAR", "_muISOBAR",  "_muHLTEND", "_muIDEND", "_muISOEND",  
            "_elHLTBAR", "_elIDBAR", "_elRECOBAR", "_elHLTEND", "_elIDEND", "_elRECOEND" };

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
      cout << "sys labels: ";
      for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
          cout << *iter << " " ;
      }
      cout<< endl;


    char dirname[40];

    int i_start=0;
    int i_max = n_m_bins;


    printf("loop start \n");
    for(int i=i_start; i<i_max; i++){

        pdf_fout->cd();
        snprintf(dirname, 10, "w%i", i);
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        m_low = m_bins[i];
        m_high = m_bins[i+1];
        printf("\n \n Start making templates for mass bin %.0f-%.0f \n", m_low, m_high);
       

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

                pdf_fout->cd();
                gDirectory->cd(dirname);
                write_out_templates(*iter);
            }
            i_sys++;
        }

    }


    pdf_fout->Close();
    printf("Templates written to %s \n", pdf_fout_name.Data());

}

