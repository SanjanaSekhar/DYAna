
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


    Double_t ratio[n_m_bins], ratio_errs[n_m_bins];
    for(int i=0; i<n_m_bins; i++){
        ratio[i] = y_comb[i]/y_powheg[i];
        ratio_errs[i] = y_comb_errs[i]/y_powheg[i];
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



    //TGraphErrors *g_sm = new TGraphErrors(n_m_bins, m, y_powheg, x_err, y_powheg_errs);
    TGraph *g_sm = new TGraphErrors(n_m_bins, m, y_powheg);
    TGraphErrors *g_comb = new TGraphErrors(n_m_bins, m, y_comb, m_err, y_comb_errs);
    TGraphErrors *g_mumu = new TGraphErrors(n_m_bins, m, y_mumu, m_err, y_mumu_errs);
    TGraphErrors *g_elel = new TGraphErrors(n_m_bins, m, y_elel, m_err, y_elel_errs);
    TGraphErrors *g_ratio = new TGraphErrors(n_m_bins, m, ratio, m_err, ratio_errs);

    g_sm->SetMarkerColor(kBlue);
    g_sm->SetLineColor(kBlue);
    g_sm->SetLineWidth(4);


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
    g_sm->Draw("ALP");
    g_sm->GetYaxis()->SetRangeUser(0.2, 0.75);
    g_sm->GetXaxis()->SetLimits(100., 1400.);
    g_sm->Draw("ALP");
    g_mumu->Draw("PE same");
    g_elel->Draw("PE same");
    g_comb->Draw("PE same");


    g_sm->GetYaxis()->SetTitle("Forward Backward Asymmetry");
    g_sm->GetYaxis()->SetNdivisions(505);
    g_sm->GetYaxis()->SetTitleSize(20);
    g_sm->GetYaxis()->SetTitleFont(43);
    //g_sm->GetYaxis()->SetTitleOffset(1.2);
    g_sm->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_sm->GetYaxis()->SetLabelSize(15);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(g_sm, "Standard Model A_{FB} from POWHEG", "l");
    leg1->AddEntry(g_mumu, "#mu#mu Measurement", "p");
    leg1->AddEntry(g_elel, "ee Measurement", "p");
    leg1->AddEntry(g_comb, "Combined Measurement", "p");
    leg1->Draw();

    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    g_ratio->Draw("APE");
    g_ratio->GetYaxis()->SetRangeUser(0.8, 1.20);
    g_ratio->GetXaxis()->SetLimits(100., 1400.);
    g_ratio->Draw("APE");


    g_ratio->GetYaxis()->SetTitle("Comb./POWHEG");
    g_ratio->GetYaxis()->SetNdivisions(505);
    g_ratio->GetYaxis()->SetTitleSize(20);
    g_ratio->GetYaxis()->SetTitleFont(43);
    //g_ratio->GetYaxis()->SetTitleOffset(1.2);
    g_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetYaxis()->SetLabelSize(15);
    // X axis g_ratio plot settings
    g_ratio->GetXaxis()->SetTitle("M_{ll} (GeV)");
    g_ratio->GetXaxis()->SetTitleSize(20);
    g_ratio->GetXaxis()->SetTitleFont(43);
    g_ratio->GetXaxis()->SetTitleOffset(3.);
    g_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g_ratio->GetXaxis()->SetLabelSize(20);
    int iPeriod = -1; 
    writeExtraText = true;
    draw_CMS = true;
    CMS_lumi(pad1, iPeriod, 11 );
    c_m->Update();
    
}


