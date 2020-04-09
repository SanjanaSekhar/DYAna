#include "../../utils/ScaleFactors.C"

void make_pu_SFs(){

    //char *f_mc_ = "2016/pileup_profile_Summer16.root";
    //char *f_data_ = "2016/PileupData_GoldenJSON_Full2016.root";

    char *f_mc_ = "2017/mcPileup2017.root";
    char *f_data_ = "2017/PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root";
    //
    //char *f_mc_ = "2018/mcPileup2018.root";
    //char *f_data_ = "2018/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root";

    char *f_out_ = "2017/pu_SF.root";

    TFile *f_data = TFile::Open(f_data_, "r");
    f_data->cd();
    TH1D *h_nom = (TH1D *) f_data->Get("pileup");
    TH1D *h_up = (TH1D *) f_data->Get("pileup_plus");
    TH1D *h_down = (TH1D *) f_data->Get("pileup_minus");
    
    h_nom->Scale(1./h_nom->Integral());
    h_up->Scale(1./h_up->Integral());
    h_down->Scale(1./h_down->Integral());


    TFile *f_mc = TFile::Open(f_mc_, "r");
    f_mc->cd();
    TH1D *h_mc = (TH1D *) gDirectory->Get("pu_mc");
    h_mc->Scale(1./h_mc->Integral());

    TFile *fout =  TFile::Open(f_out_, "RECREATE");

    fout->cd();

    TH1D *h_ratio_nom = (TH1D *) h_nom->Clone("pu_ratio");
    TH1D *h_ratio_up = (TH1D *) h_up->Clone("pu_ratio_up");
    TH1D *h_ratio_down = (TH1D *) h_down->Clone("pu_ratio_down");

    h_ratio_nom->Divide(h_mc);
    h_ratio_up->Divide(h_mc);
    h_ratio_down->Divide(h_mc);

    //h_ratio_nom->Print("all");
    //h_ratio_up->Print("all");
    //h_ratio_down->Print("all");
    fout->cd();

    h_ratio_nom->Write();
    h_ratio_up->Write();
    h_ratio_down->Write();
    fout->Close();

}

