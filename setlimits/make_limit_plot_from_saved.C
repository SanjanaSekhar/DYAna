#include "../plots/tdrstyle.C"
#include "../plots/CMS_lumi.C"
#include "find_kl_limit.C"

void make_limit_plot_from_saved(){

    bool inc_meas_limit = true;
    bool smooth_exp_bands = true;
    float smooth_const = 2.;
    char *fout = "../plots/PAS_plots/limit_obs.pdf";
    bool prelim = true;
    TFile *f_in = TFile::Open("limit_obs.root");

    setTDRStyle();

    TGraphErrors *one_sig_raw = (TGraphErrors *) f_in->Get("one_sigma_band");
    TGraphErrors *two_sig_raw = (TGraphErrors *) f_in->Get("two_sigma_band");
    TGraph *exp = (TGraph *) f_in->Get("exp_limit");
    TGraph *meas = (TGraph *) f_in->Get("obs_limit");

    Int_t n_pts = one_sig_raw->GetN();

    TGraphErrors *one_sig, *two_sig;

    if(smooth_exp_bands){
        one_sig = (TGraphErrors *) one_sig_raw->Clone("one_sigma_band_smoothed");
        two_sig = (TGraphErrors *) two_sig_raw->Clone("two_sigma_band_smoothed");
        for(int i=0; i<n_pts; i++){
            Double_t prev, next;
            Double_t std_dev = one_sig_raw->GetErrorY(i);
            if( i >0) prev = one_sig_raw->GetErrorY(i-1);
            else prev = std_dev;
            if( i <n_pts - 1) next = one_sig_raw->GetErrorY(i+1);
            else next = std_dev;

            if(std_dev > next){
                Double_t avg = (prev + next) / 2.;
                one_sig->SetPointError(i,0., avg);
                two_sig->SetPointError(i,0., 2.*avg);
            }

            //Double_t avg = (smooth_const * prev + std_dev + smooth_const *next) / (2*smooth_const + 1);

        }
    }
    else{
        one_sig = one_sig_raw;
        two_sig = two_sig_raw;
    }

    one_sig->Print();

    TLine *kl_one = new TLine(1400, 1, 4700, 1);
    kl_one->SetLineStyle(9);
    kl_one->SetLineWidth(2);
    kl_one->SetLineColor(kBlue);




    one_sig->SetFillColor(kGreen+1);
    two_sig->SetFillColor(kOrange);

    //exp->Print();
    //one_sig->Print();
    //two_sig->Print();


    exp->SetLineColor(1);
    exp->SetLineStyle(2);
    exp->SetLineWidth(4);
    auto legbb = new TLegend(0.36,0.62,0.70,0.9);
	legbb->SetFillStyle(0);
	legbb->SetFillColor(0);    
	legbb->SetBorderSize(0);  
	legbb->SetTextSize(0.040);
	legbb->SetTextFont(42);
	legbb->AddEntry(exp,"Expected","l");
	legbb->AddEntry(one_sig,"#pm 1 std. deviation","f");
	legbb->AddEntry(two_sig,"#pm 2 std. deviation","f");
    legbb->AddEntry(meas,"Observed","l");
    legbb->AddEntry(kl_one,"#kappa_{L} = 1","l");



    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    c1->SetBottomMargin(0.15);
    two_sig->SetTitle("Z' Limits");
    two_sig->GetYaxis()->SetTitle("#kappa_{L}");
    two_sig->GetYaxis()->CenterTitle();
    two_sig->GetYaxis()->SetTitleOffset(1.0);
    two_sig->GetXaxis()->SetTitle("Z' mass (GeV)");
    two_sig->GetXaxis()->SetTitleOffset(1.1);

    two_sig->Draw("AL3SAME");
    one_sig->Draw("3SAME");
    exp->Draw("LSAME");
    kl_one->Draw("LSAME");

    meas->SetLineColor(1);
    meas->SetLineStyle(1);
    meas->SetLineWidth(4);
    meas->Draw("LSAME");

    legbb->Draw("same");

    printf("idx, mass, expected limit, expected std dev obs limit \n");
    for(int idx=0; idx <n_pts; idx++){
        printf("%i %.2f %.2f %.2f %.2f \n", idx, meas->GetPointX(idx), exp->GetPointY(idx), one_sig->GetErrorY(idx), meas->GetPointY(idx));
    }


    
    if(prelim) writeExtraText = true;
    else writeExtraText = false;
    int iPeriod = -1; 
    CMS_lumi( c1, iPeriod, 11);


    //c1->Print("../plots/Paper_plots/limit_obs.png");
    c1->Print(fout);

}
