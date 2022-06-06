
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

void draw_A0(){
    setTDRStyle();

    bool draw_powheg = false;


    float chi2_sep(0.), chi2_comb(0.);
    for(int i=0; i<n_m_bins; i++){
        chi2_comb += pow((y_powheg[i] - y_comb[i])/y_comb_errs[i],2);
        chi2_sep += pow((y_powheg[i] - y_mumu[i])/y_mumu_errs[i],2);
        chi2_sep += pow((y_powheg[i] - y_elel[i])/y_elel_errs[i],2);
    }
    printf("Combined chisq is %.2f, p-value is %.3f \n", chi2_comb, ROOT::Math::chisquared_cdf_c(chi2_comb, n_m_bins));
    printf("sep chisq is %.2f, p-value is %.3f \n", chi2_sep, ROOT::Math::chisquared_cdf_c(chi2_comb, 12));




    TGraph *g_sm_pow = new TGraphErrors(n_m_bins, m, A0_powheg);

    TGraph *g_sm_amc = new TGraphErrors(n_m_ext_bins, m_ext, A0_amc);
    TGraph *g_sm_amc_unc = new TGraphErrors(n_m_ext_bins, m_ext, A0_amc, nullptr, A0_amc_errs);
    TGraphErrors *g_comb = new TGraphErrors(n_m_bins, m, A0_comb, m_err, A0_comb_errs);
    TGraphErrors *g_mumu = new TGraphErrors(n_m_bins, m, A0_mumu, m_err, A0_mumu_errs);
    TGraphErrors *g_elel = new TGraphErrors(n_m_bins, m, A0_elel, m_err, A0_elel_errs);

    g_sm_amc->Print();



    int fill_col = kGreen -10;
    g_sm_amc->SetMarkerColor(diboson_c);
    g_sm_amc_unc->SetFillColor(fill_col);
    g_sm_amc_unc->SetLineColor(fill_col);
    g_sm_amc->SetLineColor(diboson_c);
    g_sm_amc->SetFillColor(fill_col);
    g_sm_amc->SetLineWidth(3);


    g_sm_pow->SetMarkerColor(kBlue);
    g_sm_pow->SetLineColor(kBlue);
    g_sm_pow->SetLineWidth(2);


    g_comb->SetMarkerColor(kBlack);

    g_elel->SetMarkerColor(navy_c);
    g_mumu->SetMarkerColor(DY_c);
    
    g_comb->SetLineColor(kBlack);
    g_elel->SetLineColor(navy_c);
    g_mumu->SetLineColor(DY_c);

    g_comb->SetMarkerStyle(kFullSquare);
    g_mumu->SetMarkerStyle(kOpenTriangleUp);
    g_elel->SetMarkerStyle(kOpenTriangleUp);

    g_comb->SetMarkerSize(1.5);
    g_mumu->SetMarkerSize(1.5);
    g_elel->SetMarkerSize(1.5);

    g_comb->SetLineWidth(2);
    g_mumu->SetLineWidth(2);
    g_elel->SetLineWidth(2);

    g_comb->SetLineWidth(2);


    float yTOffset = 0.6;
    float LS = 0.05;

    float xTS = 0.085;
    float yTS = 0.11;


    TCanvas *c_m = new TCanvas("c_m", "Histograms", 200, 10, 1000, 800);
    //TPad *pad1 = new TPad("pad1", "pad1", 0.,0.,1.,1.);
    //pad1->SetTopMargin(0.07);
    c_m->SetBottomMargin(0.2);
    c_m->SetRightMargin(0.05);

    g_sm_amc_unc->GetYaxis()->SetRangeUser(-0.55, 0.3);
    g_sm_amc_unc->GetXaxis()->SetLimits(100., 1400.);

    g_sm_amc_unc->GetYaxis()->SetTitle("A_{0}");
    g_sm_amc_unc->GetYaxis()->SetNdivisions(505);
    g_sm_amc_unc->GetYaxis()->SetTitleSize(yTS);
    g_sm_amc_unc->GetYaxis()->SetTitleOffset(yTOffset);
    g_sm_amc_unc->GetYaxis()->SetLabelSize(LS);
    g_sm_amc_unc->GetYaxis()->CenterTitle();

    g_sm_amc_unc->GetXaxis()->SetTitle("m (GeV)");
    g_sm_amc_unc->GetXaxis()->SetTitleSize(xTS);
    g_sm_amc_unc->GetXaxis()->SetTitleOffset(0.8);
    g_sm_amc_unc->GetXaxis()->SetLabelSize(LS);



    g_sm_amc_unc->Draw("A3");
    g_sm_amc->Draw("LSAME");
    if(draw_powheg) g_sm_pow->Draw("ALP same");
    g_elel->Draw("PE same");
    g_mumu->Draw("PE same");
    g_comb->Draw("PE same");



    gStyle->SetLegendBorderSize(0);
    
    float x_size = 0.3;
    float y_size = 0.22;
    float leg_text_size = 0.035;



    gStyle->SetEndErrorSize(0);

    TLegend *leg1 = new TLegend(x_size, y_size);

    float x_start_m = 0.19;
    float y_start_m = 0.317;


    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    leg1->AddEntry(g_sm_amc, " Standard Model A_{0} from aMC@NLO", "lf");
    //leg1->AddEntry(g_sm_amc_unc, " Uncertainty on aMC@NLO", "f");
    if(draw_powheg) leg1->AddEntry(g_sm_pow, " Standard Model A_{0} from POWHEG", "l");
    leg1->AddEntry(g_mumu, " #mu#mu Measurement", "pel");
    leg1->AddEntry(g_elel, " ee Measurement", "pel");
    leg1->AddEntry(g_comb, " Combined Measurement", "pel");

    leg1->SetTextSize(leg_text_size);
    leg1->Draw();


    int iPeriod = -1; 
    writeExtraText = false;
    draw_CMS = true;
    CMS_lumi(c_m, iPeriod, 11 );
    c_m->Update();

    c_m->Print("Paper_plots/A0_mbins_unblind.pdf");
    
}


