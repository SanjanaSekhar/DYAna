

void CopyDir(TDirectory *source, TDirectory *savdir) {
   //copy all objects and subdirs of directory source as a subdir of the current directory
   //loop on all entries of this directory
   source->ls();
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(key->GetName());
         savdir->cd();
         TTree *newT = T->CloneTree(-1,"fast");
         newT->Write();
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         savdir->cd();
         obj->Write();
         delete obj;
     }
  }
  savdir->cd();
}


void LQ_merge_workspaces(){

  //  for(double i=1;i<=3;i+=0.5){
        int mLQ = 1000;
        for(int year=16;year<=18;year++){

            char *ending="071921";
            char f1_s[180];
            sprintf(f1_s,"root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/templ%i_nonsys_%s_m%i/file_0.root", year, ending, mLQ);
            char fout_s[180];
            sprintf(fout_s,"combine/mtemps/LQm%i_merge_templates%i_%s.root", mLQ, year,ending);
            //const TString f1_s("root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/templ16_nonsys_jul27LQ");
            //const TString fout_s("combine/templates/LQm%i_merge_templates%i.root");
            TFile *f1 = TFile::Open(f1_s, "READ");
            TFile *fout = TFile::Open(fout_s, "RECREATE");
            char dirname[40];


            char sys_base[180];
            sprintf(sys_base,"root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/templ%i_sys_%s_m%i", year, ending, mLQ);
            char pdf_base[180];
            sprintf(pdf_base,"root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/templ%i_pdf_%s_m%i", year, ending, mLQ);
            
            
            int num_sys_files = 25;
            int num_pdf_files = 25;

            std::vector<string> fs;
            for(int i=0; i<num_sys_files; i++){
                char fname[180];
                sprintf(fname, "%s/file_%i.root", sys_base, i);
                fs.push_back(string(fname));
            }
            for(int i=0; i<num_pdf_files; i++){
                char fname[180];
                sprintf(fname, "%s/file_%i.root", pdf_base, i);
                fs.push_back(string(fname));
            }

                snprintf(dirname, 10, "LQ");
                fout->mkdir(dirname);
                fout->cd(dirname);
                TDirectory *savedir = gDirectory;

                f1->cd(dirname);
                CopyDir(gDirectory, savedir);
                savedir->cd();
                for(auto f_name = fs.begin(); f_name !=fs.end(); f_name++){
                    TFile *f2 = TFile::Open(f_name->c_str(), "READ");
                    printf("Opening file: %s \n\n\n", f_name->c_str());
                    f2->cd();
                    f2->cd(dirname);
                    TDirectory *source_dir = gDirectory;
                    CopyDir(source_dir, savedir);
                    f2->Close();
                }

            
            fout->Print();
            fout->Write();
            fout->Close();
        }
    }
//}

