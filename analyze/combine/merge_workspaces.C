

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


void merge_workspaces(){


    const TString f1_s("combine/templates/june29_2016.root");
    const TString fout_s("combine/templates/june29_merge_2016.root");
    TFile *f1 = TFile::Open(f1_s, "READ");
    TFile *fout = TFile::Open(fout_s, "RECREATE");
    char dirname[40];


    char *sys_base = "root://cmseos.fnal.gov//store/user/oamram/Condor_outputs/templ16_sys_june30";
    char *pdf_base = "root://cmseos.fnal.gov//store/user/oamram/Condor_outputs/templ16_pdf_june30";
    int num_sys_files = 15;
    int num_pdf_files = 15;

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
   
    int i_start=0;
    int i_max = 8;



    for(int i=i_start; i<i_max; i++){

        snprintf(dirname, 10, "w%i", i);
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

    }
    fout->Print();
    fout->Write();
    fout->Close();
}

