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


void LQ_merge_workspaces2(){


    const TString f1_s("combine/templates/LQm1000_ud_templates16.root");
    const TString fout_s("combine/templates/LQm1000_ud_merge_templates16.root");
    TFile *f1 = TFile::Open(f1_s, "READ");
    TFile *fout = TFile::Open(fout_s, "RECREATE");
    char dirname[40];

/*
    char *sys_base = "root://cmseos.fnal.gov//store/user/oamram/Condor_outputs/templ18_sys_may28";
    char *pdf_base = "root://cmseos.fnal.gov//store/user/oamram/Condor_outputs/templ18_pdf_may28";
    int num_sys_files = 15;
    int num_pdf_files = 15;
*/
    std::vector<string> fs;
   
        char fname[180];
        sprintf(fname, "combine/templates/LQm1000_ud_sys0_templates16.root");
        fs.push_back(string(fname));
    
    
        char fname2[180];
        sprintf(fname2, "combine/templates/LQm1000_ud_sys1_templates16.root");
        fs.push_back(string(fname2));
    
   
    int i_start=0;
    int i_max = 8;



   

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
