#include "../../utils/ScaleFactors.C"

void make_pu_SFs(){

    int year = 2018;
    char *f_mc_ = "SFs/2018/DY18_pileup_jan21.root";
    char *f_out_ = "SFs/2018/pu_SF.root";

    pileup_systematics pu_sys;
    setup_pileup_systematic_old(&pu_sys, year);

    TFile *f_mc = TFile::Open(f_mc_, "r");
    f_mc->cd();
    TH1D *h_mc = (TH1D *) gDirectory->Get("tot_pileup");
    h_mc->Scale(1./h_mc->Integral());

    TFile *fout =  TFile::Open(f_out_, "RECREATE");

    fout->cd();

    TH1D *h_ratio_nom = new TH1D("pu_ratio", "", 100, 0., 100.);
    TH1D *h_ratio_up = new TH1D("pu_ratio_up", "", 100, 0., 100.);
    TH1D *h_ratio_down = new TH1D("pu_ratio_down", "", 100, 0., 100.);

    h_ratio_nom->Divide(pu_sys.data_pileup_nom, h_mc);
    h_ratio_up->Divide(pu_sys.data_pileup_up, h_mc);
    h_ratio_down->Divide(pu_sys.data_pileup_down, h_mc);

    //h_ratio_nom->Print("all");
    //h_ratio_up->Print("all");
    //h_ratio_down->Print("all");
    fout->cd();

    h_ratio_nom->Write();
    h_ratio_up->Write();
    h_ratio_down->Write();
    fout->Close();

}

