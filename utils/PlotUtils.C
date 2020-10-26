
#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"

float computeChi2(TH1 *h){
    // only use on ratio plots, with expected value of 1
    float sum = 0.;
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        float err = h->GetBinError(i);
        if (val > 0. && err > 0.){
            sum += std::pow((val-1)/err,2);
        }
    }
    return sum;
}

void unzero_bins(TH1 *h){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = std::max(h->GetBinContent(i), 1e-8);
        h->SetBinContent(i,val);
    }
}

void binwidth_normalize(TH1 *h){
    for (int i=1; i <= h->GetNbinsX(); i++){
        float content = h->GetBinContent(i);
        float error = h->GetBinError(i);
        float width = h->GetBinWidth(i);
        h->SetBinContent(i, content/width);
        h->SetBinError(i, error/width);
    }
}


void binwidth_normalize(THStack *h_stack){
    for(auto h: *h_stack->GetHists()){
        binwidth_normalize( (TH1 *) h);
    }
}




        





TCanvas *draw_ratio_plot(std::string title, TH1F *h, TH1F *ratio, char axis_label[80], char ratio_label[80], float ratio_min = 0.01, float ratio_max = 2.){
    TCanvas *c = new TCanvas(title.c_str(), "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad((title+"p1").c_str(), "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    h->SetMinimum(1e-5);
    h->Draw();


    c->cd();
    TPad *pad2 = new TPad((title+"p2").c_str(), "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    ratio->SetMinimum(ratio_min);
    ratio->SetMaximum(ratio_max);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    return c;
}



TCanvas* make_ratio_plot(std::string title, TH1* h1, char h1_label[80], TH1* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false, bool write_out = true, float ratio_min = 0.5, float ratio_max = 1.5){
    //ratio is done as h1/h2

    unzero_bins(h1);
    unzero_bins(h2);

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);

    h1->SetLineWidth(3);
    h2->SetLineWidth(3);

    TCanvas *c = new TCanvas(title.c_str(), "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad((title+"p1").c_str(), "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    float hmax = 1.2 * std::max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(hmax);
    h1->SetMinimum(1e-5);
    h1->Draw("hist E");
    gStyle->SetEndErrorSize(4);
    h2->Draw("hist E same");


    gStyle->SetLegendBorderSize(0);
    //TLegend *leg1 = new TLegend(0.2, 0.2);
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad((title+"p2").c_str(), "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto ratio = (TH1F *) h1->Clone("h_ratio");
    ratio->SetMinimum(ratio_min);
    ratio->SetMaximum(ratio_max);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->Divide(h2);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    //int iPeriod = 4; 
    //CMS_lumi(pad1, iPeriod, 33 );
    if(write_out) c->Print(title.c_str());
    return c;
}


std::tuple<TCanvas*, TPad*> make_stack_ratio_plot(TH1F *h_data,  THStack *h_stack, TLegend *leg, TString label, TString xlabel, 
        float hmax =-1., bool logy = true, bool logx= false){

    TCanvas *c = new TCanvas("c_" + label, "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1" + label, "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    if(logx) pad1->SetLogx();
    h_stack->Draw("hist");
    if(hmax <= 0. ) hmax = 1.2 * std::max(h_stack->GetMaximum(), h_data->GetMaximum());
    h_stack->SetMaximum(hmax);
    h_stack->SetMinimum(0.1);
    gStyle->SetEndErrorSize(4);
    h_data->SetMarkerStyle(kFullCircle);
    h_data->SetMarkerColor(1);
    h_data->DrawCopy("P E same");


    h_stack->GetYaxis()->SetTitleSize(30);
    h_stack->GetYaxis()->SetTitleFont(43);
    h_stack->GetYaxis()->SetTitleOffset(1.2);


    leg->Draw();

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    if(logx) pad2->SetLogx();


    TList *stackHists = h_stack->GetHists();
    TH1* sum = (TH1*)stackHists->At(0)->Clone();
    sum->Reset();

    for (int i=0;i<stackHists->GetSize();++i) {
      sum->Add((TH1*)stackHists->At(i));
    }
    auto h_ratio = (TH1F *) h_data->Clone("h_ratio" + label);
    float center = 1.0;
    bool do_diff = false;
    if(do_diff){
        center = 0.0;
        h_ratio->Add(sum, -1.);
    }
    h_ratio->SetMinimum(center - 0.5);
    h_ratio->SetMaximum(center + 0.5);
    h_ratio->Sumw2();
    h_ratio->SetStats(0);
    h_ratio->Divide(sum);
    h_ratio->SetMarkerStyle(21);
    h_ratio->Draw("ep");
    TLine *l1 = new TLine(150,1,2000,1);
    l1->SetLineStyle(2);
    l1->Draw();
    c->cd();

    h_ratio->SetTitle("");
    // Y axis m_ratio plot settings
   h_ratio->GetYaxis()->SetTitle("Obs/Exp");
   h_ratio->GetYaxis()->SetNdivisions(505);
   h_ratio->GetYaxis()->SetTitleSize(20);
   h_ratio->GetYaxis()->SetTitleFont(43);
   h_ratio->GetYaxis()->SetTitleOffset(1.2);
   h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetYaxis()->SetLabelSize(15);
   // X axis m_ratio plot settings
   h_ratio->GetXaxis()->SetTitle(xlabel);
   h_ratio->GetXaxis()->SetTitleSize(20);
   h_ratio->GetXaxis()->SetTitleFont(43);
   h_ratio->GetXaxis()->SetTitleOffset(3.);
   h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h_ratio->GetXaxis()->SetLabelSize(20);

   float chi2 = computeChi2(h_ratio);
   int n_bins = h_ratio->GetNbinsX();

   printf("Made ratio plot for label %s chi2/dof = %.1f/%i \n", label.Data(), chi2, n_bins);
   return std::make_pair(c, pad1);
}

