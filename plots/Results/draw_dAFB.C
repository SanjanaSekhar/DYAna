
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "Math/Functor.h"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/Colors.h"
#include "results.h"

void draw_dAFB(){
    setTDRStyle();
    char *fout = "Paper_plots/delta_AFB.pdf";
    bool prelim = false;



    float chi2_comb(0.);
    Double_t weighted_avg = 0.;
    Double_t weighted_avg_unc = 0.;
    Double_t sum_weights = 0.;

    /*
    for(int i=0; i<n_m_bins; i++){
        chi2_comb += pow((y_diff[i])/y_diff_errs[i],2);

        weighted_avg += y_diff[i]/(y_diff_errs[i]*y_diff_errs[i]);
        sum_weights += 1./(y_diff_errs[i] * y_diff_errs[i]);
    }
    printf("Chisq is %.2f, p-value is %.3f \n", chi2_comb, ROOT::Math::chisquared_cdf_c(chi2_comb, n_m_bins));

    weighted_avg /= sum_weights;
    weighted_avg_unc = pow(1./sum_weights, 0.5);
    */
    weighted_avg = -0.026;
    weighted_avg_unc = -0.0106;

    Double_t sigma = weighted_avg / weighted_avg_unc;
    printf("Weighted avg is %.3f +/- %.3f. Sigma is %.2f \n", weighted_avg, weighted_avg_unc, sigma);

    TLine *l_zero = new TLine(150, 0, 1350, 0);
    l_zero->SetLineStyle(9);
    l_zero->SetLineWidth(4);
    l_zero->SetLineColor(diboson_c);


    Double_t l_avg_vals[n_m_ext_bins];
    Double_t l_avg_unc_vals[n_m_ext_bins];
    
    for(int i=0; i< n_m_ext_bins; i++){
        l_avg_vals[i] = weighted_avg;
        l_avg_unc_vals[i] = weighted_avg_unc;
    }




    TGraphErrors *l_avg = new TGraphErrors(n_m_ext_bins, m_ext, l_avg_vals, nullptr, nullptr);
    TGraphErrors *l_avg_unc = new TGraphErrors(n_m_ext_bins, m_ext, l_avg_vals, nullptr, l_avg_unc_vals);

    TGraphErrors *g_comb = new TGraphErrors(n_m_bins, m, y_diff, m_err, y_diff_errs);

    g_comb->SetMarkerColor(kBlack);
    g_comb->SetLineColor(kBlack);
    g_comb->SetMarkerStyle(kFullSquare);
    g_comb->SetLineWidth(3);

    g_comb->SetMarkerSize(1.5);

    //TLine *l_avg = new TLine(150, weighted_avg, 1200, weighted_avg);
    //l_avg->SetLineStyle(9);
    l_avg->SetLineWidth(4);
    l_avg->SetLineColor(navy_c);
    l_avg->SetFillColor(light_blue);
    l_avg_unc->SetFillColor(light_blue);

    float yTOffset = 0.68;
    float LS = 0.05;
    float yTS = 0.1;
    float xTS = 0.075;




    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 1000, 800);
    c_m->SetRightMargin(0.05);
    c_m->SetBottomMargin(0.15);

    l_avg_unc->GetYaxis()->SetRangeUser(-0.3, 0.5);
    l_avg_unc->GetXaxis()->SetLimits(100., 1400.);

    l_avg_unc->GetYaxis()->SetTitle("#DeltaA_{FB}");
    l_avg_unc->GetYaxis()->SetNdivisions(505);
    l_avg_unc->GetYaxis()->SetTitleSize(yTS);
    l_avg_unc->GetYaxis()->SetTitleOffset(yTOffset);
    l_avg_unc->GetYaxis()->SetLabelSize(LS);
    l_avg_unc->GetYaxis()->CenterTitle();

    l_avg_unc->GetXaxis()->SetTitle("m (GeV)");
    l_avg_unc->GetXaxis()->SetTitleSize(xTS);
    l_avg_unc->GetXaxis()->SetTitleOffset(0.8);
    l_avg_unc->GetXaxis()->SetLabelSize(LS);
    l_avg_unc->GetXaxis()->SetTickLength(0.04);


    gStyle->SetEndErrorSize(0);


    l_avg_unc->Draw("A3");
    l_avg->Draw("LSAME");
    l_zero->Draw("LSAME");
    g_comb->Draw("PSAME");






    gStyle->SetLegendBorderSize(0);
    
    float x_size = 0.32;
    float y_size = 0.2;
    float leg_text_size = 0.045;



    TLegend *leg1 = new TLegend(x_size, y_size);

    float x_start_m = 0.27;
    float y_start_m = 0.60;


    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    //leg1->SetHeader("#Delta A_{FB} = A_{FB,#mu#mu} #minus A_{FB,ee}", "C");



    leg1->AddEntry(g_comb, " Mass binned measurements ", "lep");
    leg1->AddEntry(l_avg, " Inclusive measurement", "lf");
    leg1->AddEntry(l_zero, " #Delta A_{FB} = 0", "l");
    //leg1->AddEntry(l_avg_unc, " Uncertainty on inclusive measurement ", "f");

    leg1->SetTextSize(leg_text_size);
    leg1->Draw();

    TLatex latext; 
    latext.SetNDC();
    latext.SetTextColor(kBlack);
    latext.SetTextAlign(22); //centered
    latext.SetTextFont(42);
    latext.SetTextSize(0.06);    



    latext.DrawLatex(0.70, 0.88, "#Delta A_{FB} = A_{FB,#mu#mu} #minus A_{FB,ee}");




    c_m->cd();
    int iPeriod = -1; 
    if(prelim) writeExtraText = true;
    else writeExtraText = false;
    draw_CMS = true;
    CMS_lumi(c_m, iPeriod, 11 );
    c_m->Update();

    c_m->Print(fout);
    //c_m->Print("Paper_plots/delta_AFB.png");
    
}


