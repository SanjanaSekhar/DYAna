#include "LQ_make_templates.C"



void LQ_make_sys_templates(int nJobs = 1, int iJob =0, int year = 2016, int type=1, Double_t m_LQ=0.){
    //type0 is pdfs, type1 other sys
  scramble_data = false;
  //type =0;
    if(nJobs == 0){
        printf("Invalid setting of 0 total jobs! Going to change it to be 1 job (ie this process runs over all systematics) \n");
        nJobs = 1;
    }
  // for(year=2017;year<=2018;year++){
   // for(int i=1;i<=3;i+=2){
  //  m_LQ = 1000.;
   // year=2018;
    printf("=========================\n m_LQ = %f, year = %d, SYS type = %i \n=========================\n",m_LQ,year,type );
    char templates_name[100];
    sprintf(templates_name,"output_files/LQm%i_sys%i_templates%i.root",int(m_LQ),type,year%2000);
    
    const TString pdf_fout_name(templates_name);
    //const TString pdf_fout_name("output_files/sys_test.root");
    TFile *pdf_fout = TFile::Open(pdf_fout_name, "RECREATE");
    //year =2017;
    init(year);

    printf("Setting up SFs... ");
    setup_all_SFs(year);
    printf("   done \n");
    
   vector<string> sys_labels, sys_labels_raw;
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

      
      vector<string>  sys_labels_uncorr = 
        { "_BTAGCOR", "_BTAGUNCOR", "_BTAGLIGHT",  "_METJER", "_METJEC", "_METHEM", "_prefire", "_elScaleSyst", "_elScaleStat","_elScaleGain", "_elSmear", "_muRC", "_Pu", 
            "_muHLTBAR", "_muIDBAR", "_muISOBAR",  "_muHLTEND", "_muIDEND", "_muISOEND",  "_muIDSYS", "_muISOSYS",  
            "_elHLTBARPTHIGH", "_elIDBARPTHIGH", "_elRECOBARPTHIGH", "_elHLTENDPTHIGH", "_elIDENDPTHIGH", "_elRECOENDPTHIGH",
            "_elHLTBARPTLOW", "_elIDBARPTLOW", "_elRECOBARPTLOW", "_elHLTENDPTLOW", "_elIDENDPTLOW", "_elRECOENDPTLOW",
            "_ptrw1b", "_ptrw2b", "_ptrw3b", "_ptrw4b", "_ptrw5b", "_ptrw6b", "_ptrw7b",
            "_emucostrw1b", "_emucostrw2b", "_emucostrw3b", "_emucostrw4b",
            "_elfakesrw1b", "_elfakesrw2b", "_elfakesrw3b", "_elfakesrw4b",
            "_mufakesrw1b", "_mufakesrw2b", "_mufakesrw3b", "_mufakesrw4b",
            "_RENORM", "_FAC", "_REFAC","_alphaS",
        };

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

   // int i_start=0;
    //int i_max = n_m_bins;


    //printf("loop start \n");
    //for(int i=i_start; i<i_max; i++){

        pdf_fout->cd();
        snprintf(dirname, 10, "LQ");
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

        //m_low = m_bins[i];
        //m_high = m_bins[i+1];
        printf("\n \n Start making templates for LQ_sys");
       

        int i_sys =0;
        for(auto iter = sys_labels.begin(); iter !=sys_labels.end(); iter++){
            if(i_sys % nJobs == iJob){
                printf("Making MC templates for sys %s \n", (*iter).c_str());

                /*
                Double_t alpha_denom = amc_alpha[i];
                if(iter->find("alphaDenUp") != string::npos) alpha_denom = amc_alpha[i] + amc_alpha_unc[i];
                if(iter->find("alphaDenDown") != string::npos) alpha_denom = amc_alpha[i] - amc_alpha_unc[i];
                printf("alpha_denom %.2f \n", alpha_denom);
                */

                make_mc_templates(year, m_LQ,  *iter);
                 make_qcd_templates(year,  *iter);
                convert_mc_templates(year, *iter);

                pdf_fout->cd();
                gDirectory->cd(dirname);
                write_out_templates(*iter);
            }
            i_sys++;
        }

    //}


    pdf_fout->Close();
    printf("Templates written to %s \n", pdf_fout_name.Data());
  }
//}

