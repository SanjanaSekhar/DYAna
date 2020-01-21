#include "../../utils/NTupleReader.C"



void Add_Pileup(int nJobs =1, int iJob = 0, string fin = "", int year =-1)
{
    if(fin == "") fin = string("EOS_files/2017/DY_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu17_pileup.root", false);
    if (year == -1) year = 2017;
    nt.year = year;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = false;

    nt.fout->cd();

    TH1D *total_pileup = new TH1D ("tot_pileup", "", 100, 0., 100.);


    while(nt.getNextFile()){
        nt.fin->cd("EventCounter");
        TDirectory *subdir = gDirectory;
        TH1D *mc_pileup = (TH1D *)subdir->Get("pileup");
        total_pileup->Add(mc_pileup);

        total_pileup->Print();



    }

    nt.fout->cd();
    total_pileup->Write();
    nt.finish();

    return;
}

