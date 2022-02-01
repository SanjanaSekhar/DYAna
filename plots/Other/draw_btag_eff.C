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
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/HistMaker.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"



void draw_btag_eff(){
    int year = 2016;
    char *plot_dir = "Paper_plots/mc_btag_effs";
    int flag1;


    setTDRStyle();

    //TFile *f_dy = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_dy.root", "READ");
    TFile *f_dy = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_diboson.root", "READ");
    TFile *f_ttbar = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_ttbar.root", "READ");

    Double_t Eta_bins[] = {0, 0.9, 1.2, 2.1, 2.4}; 
    Int_t nEta_bins = 4;
    Double_t Pt_bins[] = {0,10,20,30,40,50,60,80,100,120,160,200,250,300,400,500};
    Int_t nPt_bins = 15;



    auto labels = {"udsg_eff", "c_eff", "b_eff"};

    for (auto & label : labels){


        TH2D *h_dy = (TH2D *)f_dy->Get(label);
        TH2D *h_ttbar = (TH2D *)f_ttbar->Get(label);

        h_dy->Print("range");
        h_ttbar->Print("range");

        TH1D *h_dy_eta1 = h_dy->ProjectionX("dy_eta1", 1,1);
        TH1D *h_ttbar_eta1 = h_ttbar->ProjectionX("ttbar_eta1", 1,1);


        TH1D *h_dy_eta2 = h_dy->ProjectionX("dy_eta2", 2,2);
        TH1D *h_ttbar_eta2 = h_ttbar->ProjectionX("ttbar_eta2", 2,2);

        TH1D *h_dy_eta3 = h_dy->ProjectionX("dy_eta3", 3,3);
        TH1D *h_ttbar_eta3 = h_ttbar->ProjectionX("ttbar_eta3", 3,3);

        TH1D *h_dy_eta4 = h_dy->ProjectionX("dy_eta4", 4,4);
        TH1D *h_ttbar_eta4 = h_ttbar->ProjectionX("ttbar_eta4", 4,4);

        h_dy_eta1->SetLineColor(kGreen + 2);
        h_ttbar_eta1->SetLineColor(kGreen -7);

        h_dy_eta2->SetLineColor(kRed );
        h_ttbar_eta2->SetLineColor(kRed -7);

        h_dy_eta3->SetLineColor(kBlue );
        h_ttbar_eta3->SetLineColor(kBlue -7);

        h_dy_eta4->SetLineColor(kMagenta + 2);
        h_ttbar_eta4->SetLineColor(kMagenta -8);



        h_dy_eta1->SetLineWidth(2);
        h_dy_eta2->SetLineWidth(2);
        h_dy_eta3->SetLineWidth(2);
        h_dy_eta4->SetLineWidth(2);

        h_ttbar_eta1->SetLineWidth(2);
        h_ttbar_eta2->SetLineWidth(2);
        h_ttbar_eta3->SetLineWidth(2);
        h_ttbar_eta4->SetLineWidth(2);

        TCanvas *c = new TCanvas("c1", "", 900, 900);
        c->SetLogy();


        h_dy_eta4->SetMaximum(1);
        h_dy_eta4->GetXaxis()->SetTitle("p_{T} (GeV)");
        h_dy_eta4->GetYaxis()->SetTitle("Efficiency");

        h_dy_eta4->Draw("hist e");
        h_dy_eta3->Draw("hist e same");
        h_dy_eta2->Draw("hist e same");
        h_dy_eta1->Draw("hist e same");

        h_ttbar_eta4->Draw("hist e same");
        h_ttbar_eta3->Draw("hist e same");
        h_ttbar_eta2->Draw("hist e same");
        h_ttbar_eta1->Draw("hist e same");



        TLegend *leg1 = new TLegend(0.3, 0.3);
        /*
        leg1->AddEntry(h_dy_eta1, "DY MC 0.0 < |y| < 0.9");
        leg1->AddEntry(h_ttbar_eta1, "ttbar MC 0.0 < |y| < 0.9");
        leg1->AddEntry(h_dy_eta2, "DY MC 0.9 < |y| < 1.2");
        leg1->AddEntry(h_ttbar_eta2, "ttbar MC 0.9 < |y| < 1.2");
        leg1->AddEntry(h_dy_eta3, "DY MC 1.2 < |y| < 2.1");
        leg1->AddEntry(h_ttbar_eta3, "ttbar MC 1.2 < |y| < 2.1");
        leg1->AddEntry(h_dy_eta4, "DY MC 2.1 < |y| < 2.4");
        leg1->AddEntry(h_ttbar_eta4, "ttbar MC 2.1 < |y| < 2.4");
        */
        leg1->AddEntry(h_dy_eta1, "Diboson MC 0.0 < |y| < 0.9");
        leg1->AddEntry(h_ttbar_eta1, "ttbar MC 0.0 < |y| < 0.9");
        leg1->AddEntry(h_dy_eta2, "Diboson MC 0.9 < |y| < 1.2");
        leg1->AddEntry(h_ttbar_eta2, "ttbar MC 0.9 < |y| < 1.2");
        leg1->AddEntry(h_dy_eta3, "Diboson MC 1.2 < |y| < 2.1");
        leg1->AddEntry(h_ttbar_eta3, "ttbar MC 1.2 < |y| < 2.1");
        leg1->AddEntry(h_dy_eta4, "Diboson MC 2.1 < |y| < 2.4");
        leg1->AddEntry(h_ttbar_eta4, "ttbar MC 2.1 < |y| < 2.4");

        leg1->SetBorderSize(0);

        leg1->Draw();
        c->Update();

        TLatex latext; 
        latext.SetNDC();
        latext.SetTextColor(kBlack);
        latext.SetTextAlign(22); //centered
        latext.SetTextFont(42);
        latext.SetTextSize(0.04);    

        float l = c->GetLeftMargin();
        float t = c->GetTopMargin();
        float r = c->GetRightMargin();
        float b = c->GetBottomMargin();

        float w = 1-l-r;

        //printf("l %.3f,r %.3f,t %.3f,b %.3f,W %.3f,H %.3f \n",l,r,t,b,W,H);


        latext.DrawLatex(l + w/2, 1 - 0.5*t, label);
        c->Update();

        char plt_name[100];
        sprintf(plt_name, "%s/mc%i_db_cmp_btag_%s_eff.png", plot_dir, year - 2000, label);
        c->Print(plt_name);


    }




}



