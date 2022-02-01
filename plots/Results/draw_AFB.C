
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

void draw_AFB(){
    setTDRStyle();
    char *fout  = "Paper_plots/AFB_mbins_unblind.pdf";
    bool prelim = false;

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



    TGraph *g_sm_pow = new TGraphErrors(n_m_ext_bins, m_ext, y_powheg);
    TGraph *g_sm_amc = new TGraphErrors(n_m_ext_bins, m_ext, y_amc);
    TGraph *g_sm_amc_unc = new TGraphErrors(n_m_bins, m_ext, y_amc, nullptr, y_amc_errs);

    TGraphErrors *g_comb = new TGraphErrors(n_m_bins, m, y_comb, m_err, y_comb_errs);
    TGraphErrors *g_mumu = new TGraphErrors(n_m_bins, m, y_mumu, m_err, y_mumu_errs);
    TGraphErrors *g_elel = new TGraphErrors(n_m_bins, m, y_elel, m_err, y_elel_errs);
    TGraphErrors *g_ratio = new TGraphErrors(n_m_bins, m, ratio, m_err, ratio_errs);

    int fill_col = kGreen -10;
    g_sm_amc->SetMarkerColor(diboson_c);
    g_sm_amc_unc->SetFillColor(fill_col);
    g_sm_amc_unc->SetLineColor(fill_col);
    g_sm_amc->SetLineColor(diboson_c);
    g_sm_amc->SetFillColor(fill_col);
    g_sm_amc->SetLineWidth(4);

    g_sm_pow->SetMarkerColor(kBlue);
    g_sm_pow->SetLineColor(kBlue);
    g_sm_pow->SetLineWidth(4);


    g_comb->SetMarkerColor(kBlack);
    g_elel->SetMarkerColor(ttbar_c);
    g_mumu->SetMarkerColor(DY_c);
    
    g_comb->SetLineColor(kBlack);
    g_elel->SetLineColor(ttbar_c);
    g_mumu->SetLineColor(DY_c);

    g_comb->SetMarkerStyle(kFullSquare);
    g_mumu->SetMarkerStyle(kOpenTriangleUp);
    g_elel->SetMarkerStyle(kOpenTriangleUp);

    g_comb->SetLineWidth(2);


    float yTOffset = 0.6;
    float LS = 0.07;

    float TS = 0.08;
    float yTS = 0.12;
    float rTS = 0.08 * 0.7/0.3;
    float ryTS = 0.07 * 0.7/0.3;
    float rLS = LS * 0.7/0.3;
    float rTOffset = 1. * 0.3 / 0.7 - 0.05;



    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 1000, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0.010);
    pad1->SetRightMargin(0.042);
    pad1->Draw();
    pad1->cd();

    g_sm_amc_unc->GetYaxis()->SetRangeUser(0.45, 0.9);
    g_sm_amc_unc->GetXaxis()->SetLimits(100., 1400.);

    g_sm_amc_unc->GetYaxis()->SetTitle("A_{FB}");
    g_sm_amc_unc->GetYaxis()->SetNdivisions(505);
    g_sm_amc_unc->GetYaxis()->SetTitleSize(yTS);
    g_sm_amc_unc->GetYaxis()->SetTitleOffset(yTOffset);
    g_sm_amc_unc->GetYaxis()->SetLabelSize(LS);
    g_sm_amc_unc->GetYaxis()->CenterTitle();


    gStyle->SetEndErrorSize(0);

    g_sm_amc_unc->Draw("A3");
    g_sm_amc->Draw("L same");
    //if(draw_powheg) g_sm_pow->Draw("L same");
    g_elel->Draw("PE same");
    g_mumu->Draw("PE same");
    g_comb->Draw("PE same");



    gStyle->SetLegendBorderSize(0);
    
    float x_size = 0.4;
    float y_size = 0.3;
    float leg_text_size = 0.045;



    TLegend *leg1 = new TLegend(x_size, y_size);

    float x_start_m = 0.4;
    float y_start_m = 0.58;


    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    leg1->AddEntry(g_sm_amc, " Standard model A_{FB} from aMC@NLO", "fl");
    //leg1->AddEntry(g_sm_amc_unc, " Uncertainty on aMC@NLO", "f");
    if(draw_powheg) leg1->AddEntry(g_sm_pow, " Standard model A_{FB} from POWHEG", "l");
    leg1->AddEntry(g_mumu, " #mu#mu measurement", "pel");
    leg1->AddEntry(g_elel, " ee measurement", "pel");
    leg1->AddEntry(g_comb, " Combined measurement", "pel");

    leg1->SetTextSize(leg_text_size);
    leg1->Draw();

    c_m->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    pad2->SetTopMargin(0.08);
    pad2->SetBottomMargin(0.5);
    pad2->SetRightMargin(0.04);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    g_ratio->Draw("APE");
    g_ratio->GetYaxis()->SetRangeUser(0.8, 1.20);
    g_ratio->GetXaxis()->SetLimits(100., 1400.);
    g_ratio->Draw("APE");


    g_ratio->GetYaxis()->SetTitle("Comb./SM");
    g_ratio->GetYaxis()->SetNdivisions(503);
    g_ratio->GetYaxis()->SetTitleSize(ryTS);
    g_ratio->GetYaxis()->SetTitleOffset(rTOffset);
    g_ratio->GetYaxis()->SetLabelSize(rLS);
    // X axis g_ratio plot settings
    g_ratio->GetXaxis()->SetTitle("m (GeV)");
    g_ratio->GetXaxis()->SetTitleSize(rTS);
    g_ratio->GetXaxis()->SetTitleOffset(1.);
    g_ratio->GetXaxis()->SetLabelSize(rLS);
    g_ratio->GetXaxis()->SetTickLength(0.04);
    int iPeriod = -1; 
    if(prelim) writeExtraText = true;
    else writeExtraText = false;
    draw_CMS = true;
    CMS_lumi(pad1, iPeriod, 11 );
    c_m->Update();

    c_m->Print(fout);
    
}


