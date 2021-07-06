
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
#include "results.h"

void draw_AFB_mbins(){
    setTDRStyle();

    bool draw_powheg = false;


    Double_t ratio[n_m_bins], ratio_errs[n_m_bins];
    for(int i=0; i<n_m_bins; i++){
        ratio[i] = y_comb[i]/y_amc[i];
        ratio_errs[i] = y_comb_errs[i]/y_amc[i];
        printf("%f ", ratio_errs[i]);
    }
    float chi2_sep(0.), chi2_comb(0.);
    for(int i=0; i<n_m_bins; i++){
        chi2_comb += pow((y_powheg[i] - y_comb[i])/y_comb_errs[i],2);
        chi2_sep += pow((y_powheg[i] - y_mumu[i])/y_mumu_errs[i],2);
        chi2_sep += pow((y_powheg[i] - y_elel[i])/y_elel_errs[i],2);
    }
    printf("Combined chisq is %.2f, p-value is %.3f \n", chi2_comb, ROOT::Math::chisquared_cdf_c(chi2_comb, n_m_bins));
    printf("sep chisq is %.2f, p-value is %.3f \n", chi2_sep, ROOT::Math::chisquared_cdf_c(chi2_comb, 12));



    TGraph *g_sm_pow = new TGraphErrors(n_m_bins, m, y_powheg);
    TGraph *g_sm_amc = new TGraphErrors(n_m_bins, m, y_amc);
    TGraphErrors *g_comb = new TGraphErrors(n_m_bins, m, y_comb, m_err, y_comb_errs);
    TGraphErrors *g_mumu = new TGraphErrors(n_m_bins, m, y_mumu, m_err, y_mumu_errs);
    TGraphErrors *g_elel = new TGraphErrors(n_m_bins, m, y_elel, m_err, y_elel_errs);
    TGraphErrors *g_ratio = new TGraphErrors(n_m_bins, m, ratio, m_err, ratio_errs);

    g_sm_amc->SetMarkerColor(kCyan-6);
    g_sm_amc->SetLineColor(kCyan-6);
    g_sm_amc->SetLineWidth(4);

    g_sm_pow->SetMarkerColor(kBlue);
    g_sm_pow->SetLineColor(kBlue);
    g_sm_pow->SetLineWidth(4);


    g_comb->SetMarkerColor(kBlack);
    g_elel->SetMarkerColor(kGreen);
    g_mumu->SetMarkerColor(kRed-7);
    
    g_comb->SetLineColor(kBlack);
    g_elel->SetLineColor(kGreen);
    g_mumu->SetLineColor(kRed-7);

    g_comb->SetMarkerStyle(kFullSquare);
    g_mumu->SetMarkerStyle(kFullSquare);
    g_elel->SetMarkerStyle(kFullSquare);

    g_comb->SetLineWidth(2);




    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 1000, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0.012);
    pad1->Draw();
    pad1->cd();
    g_sm_amc->GetYaxis()->SetRangeUser(0.2, 0.75);
    g_sm_amc->GetXaxis()->SetLimits(100., 1400.);
    g_sm_amc->Draw("ALP");
    if(draw_powheg) g_sm_pow->Draw("ALP same");
    g_elel->Draw("PE same");
    g_mumu->Draw("PE same");
    g_comb->Draw("PE same");

    int title_size = 30;
    int ratio_title_size = 20;


    g_sm_amc->GetYaxis()->SetTitle("Forward Backward Asymmetry");
    g_sm_amc->GetYaxis()->SetNdivisions(505);
    g_sm_amc->GetYaxis()->SetTitleSize(title_size);
    g_sm_amc->GetYaxis()->SetTitleFont(43);
    g_sm_amc->GetYaxis()->SetTitleOffset(1.2);
    g_sm_amc->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_sm_amc->GetYaxis()->SetLabelSize(15);

    gStyle->SetLegendBorderSize(0);
    
    float x_size = 0.4;
    float y_size = 0.3;
    float leg_text_size = 0.04;



    TLegend *leg1 = new TLegend(x_size, y_size);

    float x_start_m = 0.2;
    float y_start_m = 0.3;


    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    leg1->AddEntry(g_sm_amc, "Standard Model A_{FB} from aMC@NLO", "l");
    if(draw_powheg) leg1->AddEntry(g_sm_pow, "Standard Model A_{FB} from POWHEG", "l");
    leg1->AddEntry(g_mumu, "#mu#mu Measurement", "p");
    leg1->AddEntry(g_elel, "ee Measurement", "p");
    leg1->AddEntry(g_comb, "Combined Measurement", "p");

    leg1->SetTextSize(leg_text_size);
    leg1->Draw();

    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    g_ratio->Draw("APE");
    g_ratio->GetYaxis()->SetRangeUser(0.8, 1.20);
    g_ratio->GetXaxis()->SetLimits(100., 1400.);
    g_ratio->Draw("APE");


    g_ratio->GetYaxis()->SetTitle("Comb./aMC@NLO");
    g_ratio->GetYaxis()->SetNdivisions(505);
    g_ratio->GetYaxis()->SetTitleSize(ratio_title_size);
    g_ratio->GetYaxis()->SetTitleFont(43);
    //g_ratio->GetYaxis()->SetTitleOffset(1.2);
    g_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetYaxis()->SetLabelSize(15);
    // X axis g_ratio plot settings
    g_ratio->GetXaxis()->SetTitle("M (GeV)");
    g_ratio->GetXaxis()->SetTitleSize(title_size);
    g_ratio->GetXaxis()->SetTitleFont(43);
    g_ratio->GetXaxis()->SetTitleOffset(3.);
    g_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetXaxis()->SetLabelSize(20);
    int iPeriod = -1; 
    writeExtraText = true;
    draw_CMS = true;
    CMS_lumi(pad1, iPeriod, 11 );
    c_m->Update();

    c_m->Print("Paper_plots/AFB_mbins_blind.pdf");
    
}


