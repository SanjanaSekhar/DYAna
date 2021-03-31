#include "fit_amc_gen_cost.C"


void draw_year_cmp(){

    char *plot_dir = "Misc_plots/gen_level_year_cmp";
    char *label = "default";

    bool do_ptrw = false;
    bool do_nnpdf_unrw = false;


    //gROOT->SetBatch(1);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);

    TFile *f16_gen = TFile::Open("../analyze/output_files/DY16_gen_level_nov13.root");
    TTree *t16_gen_mu = (TTree *) f16_gen->Get("T_gen_mu");
    TTree *t16_gen_el = (TTree *) f16_gen->Get("T_gen_el");

    TFile *f17_gen = TFile::Open("../analyze/output_files/DY17_gen_level_nov13.root");
    TTree *t17_gen_mu = (TTree *) f17_gen->Get("T_gen_mu");
    TTree *t17_gen_el = (TTree *) f17_gen->Get("T_gen_el");

    TFile *f18_gen = TFile::Open("../analyze/output_files/DY18_gen_level_nov13.root");
    TTree *t18_gen_mu = (TTree *) f18_gen->Get("T_gen_mu");
    TTree *t18_gen_el = (TTree *) f18_gen->Get("T_gen_el");


    int n_m_bins = 61;
    float m_bin_size = 30.;
	float m_bin_low = 170.;
	float m_bin_high = m_bin_low + n_m_bins*m_bin_size;
    TH1F *h16_m = new TH1F("m1", "Data Dimuon Mass Distribution", n_m_bins, m_bin_low, m_bin_high);
    TH1F *h17_m = new TH1F("m2", "Data Dimuon Mass Distribution", n_m_bins, m_bin_low, m_bin_high);
    TH1F *h18_m = new TH1F("m3", "Data Dimuon Mass Distribution", n_m_bins, m_bin_low, m_bin_high);

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
    float rap_low = -100.;
    float rap_high = 100.;


    make_amc_gen_cost(t16_gen_mu,  h16_m, h16_cost, h_dummy, h16_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);
    make_amc_gen_cost(t16_gen_el,  h16_m, h16_cost, h_dummy, h16_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);

    make_amc_gen_cost(t17_gen_mu,  h17_m, h17_cost, h_dummy, h17_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);
    make_amc_gen_cost(t17_gen_el,  h17_m, h17_cost, h_dummy, h17_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);

    make_amc_gen_cost(t18_gen_mu,  h18_m, h18_cost, h_dummy, h18_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);
    make_amc_gen_cost(t18_gen_el,  h18_m, h18_cost, h_dummy, h18_pt, h_dummy, m_low, m_high, pt_low, pt_high, rap_low, rap_high, do_ptrw, do_nnpdf_unrw);


    h16_m->Scale(1./h16_m->Integral());
    h17_m->Scale(1./h17_m->Integral());
    h18_m->Scale(1./h18_m->Integral());

    h16_pt->Scale(1./h16_pt->Integral());
    h17_pt->Scale(1./h17_pt->Integral());
    h18_pt->Scale(1./h18_pt->Integral());

    h16_cost->Scale(1./h16_cost->Integral());
    h17_cost->Scale(1./h17_cost->Integral());
    h18_cost->Scale(1./h18_cost->Integral());


    h16_m->SetLineColor(kRed);
    h17_m->SetLineColor(kBlue);
    h18_m->SetLineColor(kGreen);
    h16_m->SetLineWidth(3);
    h17_m->SetLineWidth(3);
    h18_m->SetLineWidth(3);


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

    TH1F *h17_m_ratio = (TH1F * ) h17_m->Clone("h17_m_ratio");
    h17_m_ratio->Divide(h16_m);

    TH1F *h18_m_ratio = (TH1F * )h18_m->Clone("h18_m_ratio");
    h18_m_ratio->Divide(h16_m);

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

    sprintf(plt_file, "%s%s_pt_cmp.png", plot_dir, label);
    c1->Print(plt_file);

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

    sprintf(plt_file, "%s%s_cost_cmp.png", plot_dir, label);
    c2->Print(plt_file);


    TCanvas *cm = new TCanvas("cm", "", 200, 10, 900, 700);
    TPad *cpadm = new TPad("cpm", "pad1", 0.,0.3,0.98,1.);
    cpadm->SetBottomMargin(0);
    cpadm->Draw();
    cpadm->cd();
    h16_m->Draw("hist");
    h17_m->Draw("hist same");
    h18_m->Draw("hist same");

    leg1->Draw();

    

    cm->cd();
    TPad *cpad4m = new TPad("cp4m", "pad2", 0.,0,.98,0.3);
    cpad4m->SetBottomMargin(0.2);
    cpad4m->SetGridy();
    cpad4m->Draw();
    cpad4m->cd();
    h17_m_ratio->Draw("ep");
    h18_m_ratio->Draw("ep same");

    h17_m_ratio->GetYaxis()->SetTitle("Ratio wrt 2016");
    h17_m_ratio->GetXaxis()->SetTitle("cos(#theta)");
    h17_m_ratio->GetYaxis()->SetTitleSize(20);
    h17_m_ratio->GetYaxis()->SetTitleFont(43);
    h17_m_ratio->GetXaxis()->SetTitleSize(20);
    h17_m_ratio->GetXaxis()->SetTitleFont(43);
    h17_m_ratio->GetYaxis()->SetLabelSize(20);
    h17_m_ratio->GetYaxis()->SetLabelFont(43);
    h17_m_ratio->GetXaxis()->SetLabelSize(20);
    h17_m_ratio->GetXaxis()->SetLabelFont(43);
    h17_m_ratio->GetXaxis()->SetTitleOffset(3.);
    h17_m_ratio->GetYaxis()->SetTitleOffset(1.5);

    sprintf(plt_file, "%s%s_m_cmp.png", plot_dir, label);
    cm->Print(plt_file);
}

