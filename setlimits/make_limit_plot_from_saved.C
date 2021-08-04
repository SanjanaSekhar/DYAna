#include "../plots/tdrstyle.C"
#include "../plots/CMS_lumi.C"
#include "find_kl_limit.C"

void make_limit_plot_from_saved(){

    bool inc_meas_limit = true;
    bool smooth_exp_bands = true;
    float smooth_const = 2.;
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





    one_sig->SetFillColor(kGreen+1);
    two_sig->SetFillColor(kOrange);

    //exp->Print();
    //one_sig->Print();
    //two_sig->Print();


    exp->SetLineColor(1);
    exp->SetLineStyle(2);
    exp->SetLineWidth(4);
    auto legbb = new TLegend(0.2,0.6,0.6,0.9);
	legbb->SetFillStyle(1001);
	legbb->SetFillColor(0);    
	legbb->SetBorderSize(0);  
	legbb->SetTextSize(0.040);
	legbb->SetTextFont(42);
	legbb->AddEntry(exp,"Expected","l");
	legbb->AddEntry(one_sig,"#pm 1 std. deviation","f");
	legbb->AddEntry(two_sig,"#pm 2 std. deviation","f");



    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    two_sig->SetTitle("Z' Limits");
    two_sig->GetYaxis()->SetTitle("#kappa_{L}(g_{Z'}/g_{Z})");
    two_sig->GetYaxis()->CenterTitle();
    two_sig->GetYaxis()->SetTitleOffset(1.0);
    two_sig->GetXaxis()->SetTitle("Z' mass (GeV)");

    two_sig->Draw("AL3SAME");
    one_sig->Draw("3SAME");
    exp->Draw("LSAME");

    if(inc_meas_limit){
        meas->SetLineColor(1);
        meas->SetLineStyle(1);
        meas->SetLineWidth(4);
        legbb->AddEntry(meas,"Observed","l");
        meas->Draw("LSAME");
    }

    printf("idx, mass, expected limit, expected std dev obs limit \n");
    for(int idx=0; idx <n_pts; idx++){
        printf("%i %.2f %.2f %.2f %.2f \n", idx, meas->GetPointX(idx), exp->GetPointY(idx), one_sig->GetErrorY(idx), meas->GetPointY(idx));
    }


    legbb->Draw();
    
    writeExtraText = false;
    extraText = "Preliminary";
    int iPeriod = -1; 
    CMS_lumi( c1, iPeriod, 0 );


    c1->Print("limit_obs.png");
    c1->Print("limit_obs.pdf");

}
