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
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"

const int type = FLAG_MUONS;
int year = 2016;
bool write_out = true;
char *plot_dir = "Paper_plots/";



void draw_samesign_cmp(){
    setTDRStyle();
    init(year);

    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", 30, 150, 2000);

    TH1F *data_pt = new TH1F("data_pt", "MC signal", 40, 0, 1000);
    TH1F *diboson_pt = new TH1F("diboson_pt", "MC signal", 40, 0, 1000);
    TH1F *QCD_pt = new TH1F("QCD_pt", "MC signal", 40, 0, 1000);

    int xf_nbins = 16;
    TH1F *data_xf = new TH1F("data_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *diboson_xf = new TH1F("diboson_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *QCD_xf = new TH1F("QCD_xf", "MC signal", xf_nbins, 0, 0.8);



    TH1F *data_cost = new TH1F("data_cost", "Data", 10, -1.,1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", 10, -1.,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", 10, -1.,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", 10, -1.,1);




    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", 30, 150, 2000);

    TH1F *QCD_m = new TH1F("QCD_m", "QCD", 30, 150, 2000);

    TH1F *WJets_m = new TH1F("WJets_m", "WJets", 30, 150, 2000);

    TH1F *DY_m = new TH1F("DY_m", "DY (WW, WZ, ZZ)", 30, 150, 2000);
    TH1F *DY_cost = new TH1F("DY_cost", "DY (WW, WZ,ZZ)", 10, -1.,1);
    TH1F *DY_xf = new TH1F("DY_xf", "MC signal", xf_nbins, 0, 0.8);
    TH1F *DY_pt = new TH1F("DY_pt", "MC signal", 40, 0, 1000);

    TH1F *DY_phi = new TH1F("DY_phi", "", 20, -3.15, 3.15);
    TH1F *data_phi = new TH1F("data_phi", "", 20, -3.15, 3.15);
    TH1F *QCD_phi = new TH1F("QCD_phi", "", 20, -3.15, 3.15);
    TH1F *diboson_phi = new TH1F("diboson_phi", "", 20, -3.15, 3.15);

    TH1F * dummy = new TH1F("dummy", "", 100, 0, 100.);


    diboson_m->SetFillColor(kGreen+3);
    diboson_cost->SetFillColor(kGreen + 3);
    diboson_xf->SetFillColor(kGreen+3);
    diboson_pt->SetFillColor(kGreen+3);
    diboson_phi->SetFillColor(kGreen+3);

    QCD_xf->SetFillColor(kRed -7);
    QCD_m->SetFillColor(kRed -7);
    QCD_cost->SetFillColor(kRed -7);
    QCD_pt->SetFillColor(kRed -7);
    QCD_phi->SetFillColor(kRed -7);


    DY_xf->SetFillColor(kRed+1);
    DY_pt->SetFillColor(kRed+1);
    DY_m->SetFillColor(kRed+1);
    DY_cost->SetFillColor(kRed+1);
    DY_phi->SetFillColor(kRed+1);

    int m_low = 150.;
    int m_high = 10000.;
    bool ss = true;
    bool in_os_region = false;


    make_m_cost_pt_xf_hist(t_mumu_ss_data, data_m, data_cost, data_pt, data_xf, data_phi, dummy, true, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, dummy, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_ttbar, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, dummy, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_wt, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, dummy, false, type,   year, m_low, m_high, ss);
    make_m_cost_pt_xf_hist(t_mumu_ss_dy, DY_m, DY_cost, DY_pt, DY_xf, DY_phi, dummy, false, type,   year, m_low, m_high, ss);

    //Fakerate_est_mu(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, year, m_low, m_high, ss, in_os_region);
    make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, QCD_phi, dummy, type, year, m_low, m_high, ss, in_os_region);

    printf("Integrals of data, QCD, diboson, DY are %.2f %.2f %.2f %.2f \n", data_m->Integral(), QCD_m->Integral(), diboson_m->Integral(), DY_m->Integral());
    printf("Integrals of data, QCD, diboson, DY are %.2f %.2f %.2f %.2f \n", data_cost->Integral(), QCD_cost->Integral(), diboson_cost->Integral(), DY_cost->Integral());

    Double_t n_data = data_m->Integral();
    Double_t n_est = diboson_m->Integral() + QCD_m->Integral() + DY_m->Integral();
    bool normalize = false;
    
    if(normalize){
        Double_t n_data = data_m->Integral();
        Double_t n_mc = diboson_m->Integral() +  DY_m->Integral();
        Double_t n_QCD = QCD_m->Integral();
        Double_t qcd_ratio = (n_data - n_mc) / n_QCD;
        printf("Ratio of obs to expected QCD is %.2f \n", qcd_ratio);


        QCD_m->Scale(qcd_ratio);
        QCD_pt->Scale(qcd_ratio);
        QCD_cost->Scale(qcd_ratio);
        QCD_xf->Scale(qcd_ratio);
        QCD_phi->Scale(qcd_ratio);
    }



    bool scale_error=true;
    float qcd_err = 0.5;
    float diboson_err = 0.05;
    bool add_err = true;
    if(scale_error){
        setHistError(QCD_m, qcd_err, add_err);
        setHistError(QCD_cost, qcd_err, add_err);
        setHistError(QCD_xf, qcd_err, add_err);
        setHistError(QCD_phi, qcd_err, add_err);

        setHistError(diboson_m, diboson_err, add_err);
        setHistError(diboson_cost, diboson_err, add_err);
        setHistError(diboson_xf, diboson_err, add_err);
        setHistError(diboson_phi, diboson_err, add_err);

    }
    




    int nBins_x = QCD_m->GetXaxis()->GetNbins();
    int nBins_y = QCD_cost->GetYaxis()->GetNbins();
    //printf("Get size %i \n", nBins);





    

    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(DY_m);



    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; MuMu Cos(#theta)_{r}");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(DY_cost);

    THStack *pt_stack = new THStack("pt_stack", "Dimuon Pt Distribution: Data vs MC; Dimuon Pt (GeV)");
    pt_stack->Add(diboson_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(DY_pt);

    THStack *xf_stack = new THStack("xf_stack", "Dimuon x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(diboson_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(DY_xf);

    THStack *phi_stack = new THStack("phi_stack", "Dimuon phi Distribution: Data vs MC; #phi");
    phi_stack->Add(diboson_phi);
    phi_stack->Add(QCD_phi);
    phi_stack->Add(DY_phi);


    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.3, 0.3);
    leg1->AddEntry(data_pt, "data", "p");
    leg1->AddEntry(DY_pt, "DY (miss-sign)", "f");
    leg1->AddEntry(QCD_pt, "QCD + WJets", "f");
    leg1->AddEntry(diboson_m, "t#bar{t} + wt + WW + WZ + ZZ", "f");
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");
    TLegend *leg4 = (TLegend *) leg1->Clone("leg4");
    TLegend *leg5 = (TLegend *) leg1->Clone("leg5");



    TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    writeExtraText = false;
    char plt_file[100];


    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{#mu#mu} (GeV)", -1., true);
    CMS_lumi(p_m, year, 33 );
    sprintf(plt_file, "%sMuMu%i_ss_m_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "Cos(#theta_r)", 500., false);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sMuMu%i_ss_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dimuon pt (GeV)", -1., true);
    CMS_lumi(p_pt, year, 33);

    std::tie(c_xf, p_xf) = make_stack_ratio_plot(data_xf, xf_stack, leg4, "xf", "x_F (GeV)", -1., true);
    CMS_lumi(p_xf, year, 33);

    std::tie(c_phi, p_phi) = make_stack_ratio_plot(data_phi, phi_stack, leg5, "phi", "dimuon #phi", -1., true);
    CMS_lumi(p_phi, year, 33);



 
}

    
    
