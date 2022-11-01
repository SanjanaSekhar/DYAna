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
#include "TRatioPlot.h"
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
#include "../../utils/root_files.h"

const int type = FLAG_MUONS;



void draw_pileup(){
    setTDRStyle();
    for(int year = 2016; year<=2018; year++){

    init(year);
    init_indv_bkgs(year);

    char mu_fname1[100], mu_fname2[100];
    sprintf(mu_fname1,"AN_plots/Mu%i_pileup_before.png",year-2000);
    sprintf(mu_fname2,"AN_plots/Mu%i_pileup_after.png",year-2000);
    TFile *f7;

    if(year == 2016) f7 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_69200.root");
    else if(year == 2017) f7 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_69200.root");
    else if(year == 2018) f7 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_69200.root");

    TH1D *data_pu = (TH1D *) f7->Get("pileup")->Clone();
    data_pu->Scale(1./data_pu->Integral());
    data_pu->SetDirectory(0);

    TH1F *mc_pu_before = new TH1F("mc_pu_before", "MC signal", 100, 0, 100);
    TH1F *mc_nosig_pu_before = new TH1F("mc_nosig_pu_before", "MC signal", 100, 0, 100);
    TH1F *ttbar_pu_before = new TH1F("ttbar_pu_before", "MC signal", 100, 0, 100);
    TH1F *diboson_pu_before = new TH1F("diboson_pu_before", "MC signal", 100, 0, 100);
    TH1F *wt_pu_before = new TH1F("wt_pu_before", "MC signal", 100, 0, 100);

    TH1F *mc_pu_after = new TH1F("mc_pu_after", "MC signal", 100, 0, 100);
    TH1F *mc_nosig_pu_after = new TH1F("mc_nosig_pu_after", "MC signal", 100, 0, 100);
    TH1F *ttbar_pu_after = new TH1F("ttbar_pu_after", "MC signal", 100, 0, 100);
    TH1F *diboson_pu_after = new TH1F("diboson_pu_after", "MC signal", 100, 0, 100);
    TH1F *wt_pu_after = new TH1F("wt_pu_after", "MC signal", 100, 0, 100);

    mc_pu_before->SetFillColor(kRed+1);
    mc_pu_before->SetMarkerColor(kRed+1);
    mc_nosig_pu_before->SetFillColor(kMagenta);
    ttbar_pu_before->SetFillColor(kBlue);
    ttbar_pu_before->SetMarkerStyle(21);
    ttbar_pu_before->SetMarkerColor(kBlue);
    diboson_pu_before->SetFillColor(kGreen+3);
    wt_pu_before->SetFillColor(kOrange+7); 

    mc_pu_after->SetFillColor(kRed+1);
    mc_pu_after->SetMarkerColor(kRed+1);
    mc_nosig_pu_after->SetFillColor(kMagenta);
    ttbar_pu_after->SetFillColor(kBlue);
    ttbar_pu_after->SetMarkerStyle(21);
    ttbar_pu_after->SetMarkerColor(kBlue);
    diboson_pu_after->SetFillColor(kGreen+3);
    wt_pu_after->SetFillColor(kOrange+7); 







    make_pileup_hist(t_mumu_mc, mc_pu_before, mc_pu_after, false, year, type);
    //make_pileup_hist(t_mumu_nosig, mc_nosig_pu_before, mc_nosig_pu_after, false, year, type);
    make_pileup_hist(t_mumu_ttbar, ttbar_pu_before, ttbar_pu_after, false, year, type);
    make_pileup_hist(t_mumu_wt, wt_pu_before, wt_pu_after, false,year, type);
    make_pileup_hist(t_mumu_diboson, diboson_pu_before, diboson_pu_after, false, year, type);









    Double_t pu_before_tot = ttbar_pu_before->Integral() + wt_pu_before->Integral() + diboson_pu_before->Integral() + 
                             mc_pu_before->Integral();

    printf("pu_before_tot: %.3e \n", pu_before_tot);
    printf("%.2e %.2e %.2e %.2e  \n", ttbar_pu_before->Integral(), wt_pu_before->Integral(), diboson_pu_before->Integral(),  mc_pu_before->Integral());
    ttbar_pu_before->Scale(1./pu_before_tot);
    wt_pu_before->Scale(1./pu_before_tot);
    diboson_pu_before->Scale(1./pu_before_tot);
    //mc_nosig_pu_before->Scale(1./pu_before_tot);
    mc_pu_before->Scale(1./pu_before_tot);

    Double_t pu_after_tot = ttbar_pu_after->Integral() + wt_pu_after->Integral() + diboson_pu_after->Integral() + 
                             mc_pu_after->Integral();
    printf("pu_after_tot: %.3e \n", pu_after_tot);

    ttbar_pu_after->Scale(1./pu_after_tot);
    wt_pu_after->Scale(1./pu_after_tot);
    diboson_pu_after->Scale(1./pu_after_tot);
    //mc_nosig_pu_after->Scale(1./pu_after_tot);
    mc_pu_after->Scale(1./pu_after_tot);

    printf(" %.3e %.3e \n", mc_pu_after->Integral(), ttbar_pu_after->Integral());

    THStack *pu_before_stack = new THStack("pu_before_stack", "DiMuon Pileup; # Primary Vertices");
    pu_before_stack->Add(ttbar_pu_before);
    pu_before_stack->Add(wt_pu_before);
    pu_before_stack->Add(diboson_pu_before);
    //pu_before_stack->Add(mc_nosig_pu_before);
    pu_before_stack->Add(mc_pu_before);

    THStack *pu_after_stack = new THStack("pu_after_stack", "Dimuon Pileup; # Primary Vertices");
    pu_after_stack->Add(ttbar_pu_after);
    pu_after_stack->Add(wt_pu_after);
    pu_after_stack->Add(diboson_pu_after);
    //pu_after_stack->Add(mc_nosig_pu_after);
    pu_after_stack->Add(mc_pu_after);


    TCanvas *c_pu_before = new TCanvas("c_pu_before", "Histograms", 200, 10, 900, 700);
    c_pu_before->cd();
    TPad *before_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    before_pad1->SetBottomMargin(0);
    before_pad1->Draw();
    before_pad1->cd();
    pu_before_stack->Draw("hist");
    data_pu->SetMarkerStyle(kFullCircle);
    data_pu->SetMarkerColor(1);
    //pu_before_stack->SetMinimum(1);
    //pu_before_stack->SetMaximum(100000);
    //data_pu->SetMinimum(1);
    //data_pu->SetMaximum(100000);
    data_pu->Draw("P same");
    before_pad1->SetLogy();
    before_pad1->Update();
    TLegend *leg3 = new TLegend(0.65, 0.65, 0.9, 0.8);
    leg3->AddEntry(data_pu, "data", "p");
    leg3->AddEntry(mc_pu_before, "DY (q#bar{q}, qg #bar{q}g)", "f");
    //leg3->AddEntry(mc_nosig_pu_before, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg3->AddEntry(diboson_pu_before, "WW + WZ + ZZ", "f");
    leg3->AddEntry(wt_pu_before, "tW + #bar{t}W", "f");
    leg3->AddEntry(ttbar_pu_before, "t#bar{t}", "f");
    leg3->Draw();
    c_pu_before->cd();

    TPad *before_pad2 = new TPad("before_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    before_pad2->SetBottomMargin(0.2);
    before_pad2->SetGridy();
    before_pad2->Draw();
    before_pad2->cd();
    TList *before_stackHists = pu_before_stack->GetHists();
    TH1* before_mc_sum = (TH1*) before_stackHists->At(0)->Clone();
    before_mc_sum->Reset();

    for (int i=0;i<before_stackHists->GetSize();++i) {
      before_mc_sum->Add((TH1*)before_stackHists->At(i));
    }
    auto before_ratio = (TH1F *) data_pu->Clone("h_before_ratio");
    before_ratio->SetMinimum(0.7);
    before_ratio->SetMaximum(1.3);
    before_ratio->Sumw2();
    before_ratio->SetStats(0);
    before_ratio->Divide(before_mc_sum);
    before_ratio->SetMarkerStyle(21);
    before_ratio->Draw("ep");
    c_pu_before->cd();

    before_ratio->SetTitle("");
    // Y axis before_ratio plot settings
   before_ratio->GetYaxis()->SetTitle("Data/MC");
   before_ratio->GetYaxis()->SetNdivisions(505);
   before_ratio->GetYaxis()->SetTitleSize(20);
   before_ratio->GetYaxis()->SetTitleFont(43);
   before_ratio->GetYaxis()->SetTitleOffset(1.2);
   before_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   before_ratio->GetYaxis()->SetLabelSize(15);
   // X axis before_ratio plot settings
   before_ratio->GetXaxis()->SetTitle("Number of vertices before reweighting");
   before_ratio->GetXaxis()->SetTitleSize(20);
   before_ratio->GetXaxis()->SetTitleFont(43);
   before_ratio->GetXaxis()->SetTitleOffset(3.);
   before_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   before_ratio->GetXaxis()->SetLabelSize(20);

    int iPeriod = 4; 
    writeExtraText = false;
    CMS_lumi(before_pad1, year, 33 );
    c_pu_before->Update();
    c_pu_before->Print(mu_fname1);

    TCanvas *c_pu_after = new TCanvas("c_pu_after", "Histograms", 200, 10, 900, 700);
    c_pu_after->cd();
    TPad *after_pad1 = new TPad("pad1", "pad1", 0.,0.3,0.98,1.);
    after_pad1->SetBottomMargin(0);
    after_pad1->Draw();
    after_pad1->cd();
    pu_after_stack->Draw("hist");
    data_pu->SetMarkerStyle(kFullCircle);
    data_pu->SetMarkerColor(1);
    //pu_after_stack->SetMinimum(1);
    //pu_after_stack->SetMaximum(100000);
    //data_pu->SetMinimum(1);
    //data_pu->SetMaximum(100000);
    data_pu->Draw("P  same");
    after_pad1->SetLogy();
    after_pad1->Update();
    TLegend *leg4 = new TLegend(0.65, 0.65, 0.9, 0.8);
    leg4->AddEntry(data_pu, "data", "p");
    leg4->AddEntry(mc_pu_after, "DY (q#bar{q}, qg #bar{q}g)", "f");
    //leg4->AddEntry(mc_nosig_pu_after, "DY no asymmety(gg, qq, #bar{q}#bar{q})", "f");
    leg4->AddEntry(diboson_pu_after, "WW + WZ + ZZ", "f");
    leg4->AddEntry(wt_pu_after, "tW + #bar{t}W", "f");
    leg4->AddEntry(ttbar_pu_after, "t#bar{t}", "f");
    leg4->Draw();


    writeExtraText = false;
    CMS_lumi(after_pad1, year, 33 );

    c_pu_after->cd();


    TPad *after_pad2 = new TPad("after_pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    after_pad2->SetBottomMargin(0.2);
    after_pad2->SetGridy();
    after_pad2->Draw();
    after_pad2->cd();
    TList *after_stackHists = pu_after_stack->GetHists();
    TH1* after_mc_sum = (TH1*)after_stackHists->At(0)->Clone();
    after_mc_sum->Reset();

    for (int i=0;i<after_stackHists->GetSize();++i) {
      after_mc_sum->Add((TH1*)after_stackHists->At(i));
    }
    auto after_ratio = (TH1F *) data_pu->Clone("h_after_ratio");
    after_ratio->SetMinimum(0.7);
    after_ratio->SetMaximum(1.3);
    after_ratio->Sumw2();
    after_ratio->SetStats(0);
    after_ratio->Divide(after_mc_sum);
    after_ratio->SetMarkerStyle(21);
    after_ratio->Draw("ep");
    c_pu_after->cd();

    after_ratio->SetTitle("Number of Vertices after reweighting");
    // Y axis after_ratio plot settings
   after_ratio->GetYaxis()->SetTitle("Data/MC");
   after_ratio->GetYaxis()->SetNdivisions(505);
   after_ratio->GetYaxis()->SetTitleSize(20);
   after_ratio->GetYaxis()->SetTitleFont(43);
   after_ratio->GetYaxis()->SetTitleOffset(1.2);
   after_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   after_ratio->GetYaxis()->SetLabelSize(15);
   // X axis after_ratio plot settings
   after_ratio->GetXaxis()->SetTitle("Number of vertices after reweighting");
   after_ratio->GetXaxis()->SetTitleSize(20);
   after_ratio->GetXaxis()->SetTitleFont(43);
   after_ratio->GetXaxis()->SetTitleOffset(3.);
   after_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   after_ratio->GetXaxis()->SetLabelSize(20);
	c_pu_after->Update();
	c_pu_after->Print(mu_fname2);

	}
}



