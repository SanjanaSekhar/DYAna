#include "fit_amc_gen_cost.C"


void draw_pt_cmp(){

    //gROOT->SetBatch(1);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);

    TFile *f16_gen = TFile::Open("../analyze/output_files/DY16_gen_level_april17.root");
    TTree *t16_gen_mu = (TTree *) f16_gen->Get("T_gen_mu");
    TTree *t16_gen_el = (TTree *) f16_gen->Get("T_gen_el");

    TFile *f17_gen = TFile::Open("../analyze/output_files/DY17_gen_level_april17.root");
    TTree *t17_gen_mu = (TTree *) f17_gen->Get("T_gen_mu");
    TTree *t17_gen_el = (TTree *) f17_gen->Get("T_gen_el");

    TFile *f18_gen = TFile::Open("../analyze/output_files/DY18_gen_level_april17.root");
    TTree *t18_gen_mu = (TTree *) f18_gen->Get("T_gen_mu");
    TTree *t18_gen_el = (TTree *) f18_gen->Get("T_gen_el");

    TH1F *h16_cost = new TH1F("h16_cost", "", 40, -1., 1.);
    TH1F *h17_cost = new TH1F("h17_cost", "", 40,-1., 1.);
    TH1F *h18_cost = new TH1F("h18_cost", "", 40, -1., 1.);

    TH1F *h16_pt = new TH1F("h16_pt", "", 40, 0., 400.);
    TH1F *h17_pt = new TH1F("h17_pt", "", 40, 0., 400.);
    TH1F *h18_pt = new TH1F("h18_pt", "", 40, 0., 400.);

    TH1F *h_dummy = new TH1F("h_dummy", "", 100, 0, 100);
    float m_low = 150.;
    float m_high = 10000;
    float pt_low = 0.;
    float pt_high = 10000.;


    make_amc_gen_cost(t16_gen_mu,  h16_cost, h_dummy, h16_pt, h_dummy, m_low, m_high, pt_low, pt_high);
    make_amc_gen_cost(t16_gen_el,  h16_cost, h_dummy, h16_pt, h_dummy, m_low, m_high, pt_low, pt_high);

    make_amc_gen_cost(t17_gen_mu,  h17_cost, h_dummy, h17_pt, h_dummy, m_low, m_high, pt_low, pt_high);
    make_amc_gen_cost(t17_gen_el,  h17_cost, h_dummy, h17_pt, h_dummy, m_low, m_high, pt_low, pt_high);

    make_amc_gen_cost(t18_gen_mu,  h18_cost, h_dummy, h18_pt, h_dummy, m_low, m_high, pt_low, pt_high);
    make_amc_gen_cost(t18_gen_el,  h18_cost, h_dummy, h18_pt, h_dummy, m_low, m_high, pt_low, pt_high);

    h16_pt->Scale(1./h16_pt->Integral());
    h17_pt->Scale(1./h17_pt->Integral());
    h18_pt->Scale(1./h18_pt->Integral());

    h16_cost->Scale(1./h16_cost->Integral());
    h17_cost->Scale(1./h17_cost->Integral());
    h18_cost->Scale(1./h18_cost->Integral());


    h16_pt->SetLineColor(kRed);
    h17_pt->SetLineColor(kBlue);
    h18_pt->SetLineColor(kGreen);
    h16_pt->SetLineWidth(3);
    h17_pt->SetLineWidth(3);
    h18_pt->SetLineWidth(3);

    h16_cost->SetLineColor(kRed);
    h17_cost->SetLineColor(kBlue);
    h18_cost->SetLineColor(kGreen);

    h16_cost->SetLineWidth(3);
    h17_cost->SetLineWidth(3);
    h18_cost->SetLineWidth(3);

    TH1F *h17_pt_ratio = (TH1F * ) h17_pt->Clone("h17_pt_ratio");
    h17_pt_ratio->Divide(h16_pt);

    TH1F *h18_pt_ratio = (TH1F * )h18_pt->Clone("h18_pt_ratio");
    h18_pt_ratio->Divide(h16_pt);

    TH1F *h17_cost_ratio = (TH1F * )h17_cost->Clone("h17_cost_ratio");
    h17_cost_ratio->Divide(h16_cost);

    TH1F *h18_cost_ratio = (TH1F * )h18_cost->Clone("h18_cost_ratio");
    h18_cost_ratio->Divide(h16_cost);

    TCanvas *c1 = new TCanvas("c1", "", 200, 10, 900, 700);
    TPad *cpad1 = new TPad("cp1", "pad1", 0.,0.3,0.98,1.);
    cpad1->SetBottomMargin(0);
    cpad1->Draw();
    cpad1->cd();
    cpad1->SetLogy();
    h16_pt->Draw("hist");
    h17_pt->Draw("hist same");
    h18_pt->Draw("hist same");

    TLegend *leg1 = new TLegend(0.2, 0.2);
    leg1->AddEntry(h16_pt, "2016");
    leg1->AddEntry(h17_pt, "2017");
    leg1->AddEntry(h18_pt, "2018");

    leg1->Draw();

    

    c1->cd();
    TPad *cpad2 = new TPad("cp2", "pad2", 0.,0,.98,0.3);
    cpad2->SetBottomMargin(0.2);
    cpad2->SetGridy();
    cpad2->Draw();
    cpad2->cd();
    h17_pt_ratio->Draw("ep");
    h18_pt_ratio->Draw("ep same");

    h17_pt_ratio->GetYaxis()->SetTitle("Ratio wrt 2016");
    h17_pt_ratio->GetXaxis()->SetTitle("pT (GeV)");
    h17_pt_ratio->GetYaxis()->SetTitleSize(20);
    h17_pt_ratio->GetYaxis()->SetLabelFont(43);
    h17_pt_ratio->GetYaxis()->SetLabelSize(20);
    h17_pt_ratio->GetYaxis()->SetTitleFont(43);
    h17_pt_ratio->GetXaxis()->SetTitleSize(20);
    h17_pt_ratio->GetXaxis()->SetLabelFont(43);
    h17_pt_ratio->GetXaxis()->SetLabelSize(20);
    h17_pt_ratio->GetXaxis()->SetTitleFont(43);
    h17_pt_ratio->GetXaxis()->SetTitleOffset(3.);
    h17_pt_ratio->GetYaxis()->SetTitleOffset(1.5);


    TCanvas *c2 = new TCanvas("c2", "", 200, 10, 900, 700);
    TPad *cpad3 = new TPad("cp3", "pad1", 0.,0.3,0.98,1.);
    cpad3->SetBottomMargin(0);
    cpad3->Draw();
    cpad3->cd();
    h16_cost->Draw("hist");
    h17_cost->Draw("hist same");
    h18_cost->Draw("hist same");

    leg1->Draw();

    

    c2->cd();
    TPad *cpad4 = new TPad("cp4", "pad2", 0.,0,.98,0.3);
    cpad4->SetBottomMargin(0.2);
    cpad4->SetGridy();
    cpad4->Draw();
    cpad4->cd();
    h17_cost_ratio->Draw("ep");
    h18_cost_ratio->Draw("ep same");

    h17_cost_ratio->GetYaxis()->SetTitle("Ratio wrt 2016");
    h17_cost_ratio->GetXaxis()->SetTitle("cos(#theta)");
    h17_cost_ratio->GetYaxis()->SetTitleSize(20);
    h17_cost_ratio->GetYaxis()->SetTitleFont(43);
    h17_cost_ratio->GetXaxis()->SetTitleSize(20);
    h17_cost_ratio->GetXaxis()->SetTitleFont(43);
    h17_cost_ratio->GetYaxis()->SetLabelSize(20);
    h17_cost_ratio->GetYaxis()->SetLabelFont(43);
    h17_cost_ratio->GetXaxis()->SetLabelSize(20);
    h17_cost_ratio->GetXaxis()->SetLabelFont(43);
    h17_cost_ratio->GetXaxis()->SetTitleOffset(3.);
    h17_cost_ratio->GetYaxis()->SetTitleOffset(1.5);
}

