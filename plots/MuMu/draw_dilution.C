

#define STAND_ALONE

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
#include "../../utils/HistMaker.C"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"






void draw_dilution(){
    //setTDRStyle();
    gStyle->SetOptStat(0);

    const int nBins_xf = 10;
    Float_t xf_bins[] = {0., 0.02, 0.04, 0.07, 1.0};

    const int nBins_y = 10;
    Float_t y_bins[] = {0., 1., 1.25, 1.5,  2.4};
    TH1F *h_y = new TH1F("h_y", "Y distribution: Mass 150-171 GeV", nBins_y, 0., 2.4);
    TH1F *h_xF = new TH1F("h_xf", "xF distribution: Mass 150-171 GeV", nBins_xf, 0., 0.5);
    TH1F *h_Nc = new TH1F("h_Nc", "Number Correct; |xF|", nBins_y, 0., 2.4);
    TH1F *h_Ni = new TH1F("h_Ni", "Number Incorrect", nBins_y, 0., 2.4);
    h_Nc->Sumw2();
    h_Ni->Sumw2();


    //read event data
    int year = 2018;
    init_mc(year);

    float m_low = 700.;
    float m_high = 1000.;


    TempMaker tm(t_mumu_mc, false, year);

    tm.do_muons = true;
    tm.is_gen_level = true;
    tm.setup();

    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < 50.  && tm.has_no_bjets && tm.not_cosmic;
        if(pass){

            tm.doCorrections();
            tm.getEvtWeight();
            Double_t ratio = tm.cost_st/tm.cost;
            double rap = std::abs(tm.cm.Rapidity());
            h_y->Fill(rap, tm.evt_weight);
            h_xF->Fill(tm.xF, tm.evt_weight);
            if(ratio > 0) h_Nc->Fill(rap, tm.evt_weight);
            if(ratio < 0) h_Ni->Fill(rap, tm.evt_weight);


        }
    }
    h_y->Scale(1./h_y->Integral());
    h_xF->Scale(1./h_xF->Integral());

    TCanvas *c0 = new TCanvas("c0", "canvas", 200,10, 900,700);
    h_y->SetFillColor(kBlue);
    h_y->SetMinimum(0);
    h_y->GetXaxis()->SetTitle("|y|");
    h_y->Draw("hist");
    c0->Update();

    TCanvas *c1 = new TCanvas("c1", "canvas", 200,10, 900,700);
    h_xF->SetFillColor(kBlue);
    h_xF->SetMinimum(0);
    h_xF->GetXaxis()->SetTitle("|xF|");
    h_xF->Draw("hist");
    c1->Update();


    TCanvas *c2 = new TCanvas("c2", "canva", 100,100, 700,700);
    h_Nc ->SetLineColor(kBlue);
    h_Nc -> SetTitle("Guessing Lepton Pair Direction as Incident Quark Direction");
    h_Nc->SetStats(kFALSE);
    h_Nc ->SetLineWidth(2);
    h_Ni ->SetLineColor(kRed);
    h_Ni ->SetLineWidth(2);
    h_Nc->SetMinimum(0);
    h_Nc ->Draw("hist");
    h_Ni ->Draw("hist same");

    TLegend *leg1 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg1->AddEntry(h_Nc, "c_{r} Correct Sign", "f");
    leg1->AddEntry(h_Ni, "c_{r} Incorrect Sign", "f");
    leg1->Draw();
    c2->Update();

    Double_t dilu[nBins_y], dilu_error[nBins_y], xF_center[nBins_y]; 

    for (int i=1; i <= nBins_y; i++){
        double Nc = h_Nc->GetBinContent(i);
        double Ni = h_Ni->GetBinContent(i);

        double Nc_e = h_Nc->GetBinError(i);
        double Ni_e = h_Ni->GetBinError(i);
        
        xF_center[i] = h_Nc->GetBinCenter(i);
        dilu[i-1] = (Nc - Ni)/(Nc + Ni);

        dilu_error[i-1] = sqrt(pow(2*Nc_e * Ni / pow(Nc + Ni,2),2) + pow(2*Ni_e * Nc / pow(Nc + Ni,2),2));
        printf("Num events %.1f, xf %.2f Dilu %.2f error %1.3e \n", Nc+Ni, xF_center[i],
                dilu[i-1], dilu_error[i-1]);
    }

    TGraphErrors *g_dillu = new TGraphErrors(nBins_xf, xF_center, dilu, 0, dilu_error);

    TCanvas *c3 = new TCanvas("c3", "canvas", 200,10, 900,700);
    g_dillu->Draw("A C P");
    g_dillu->SetMarkerStyle(20);
    g_dillu->SetTitle("Dilution Effect: Mass 700-1000 GeV");
    g_dillu->SetMinimum(0);
    g_dillu->GetXaxis()->SetTitle("|y|");
    g_dillu->GetYaxis()->SetTitle("Dilution Factor");

    c3->Update();


    /*
    lumiTextSize     = 0.2;
    lumiTextOffset   = 0.2;
    cmsTextSize      = 0.35;
    cmsTextOffset    = 0.2;  // only used in outOfFrame version
    writeExtraText = true;
    extraText = "Simulation";
    lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    int iPeriod = 0; 
    CMS_lumi( c3, iPeriod, 11 );
    CMS_lumi( c2, iPeriod, 33 );
    */




}
