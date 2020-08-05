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
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TLine.h"
#include "TTree.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "../../utils/root_files.h"
#include "../../utils/HistUtils.C"
#include "../../utils/PlotUtils.C"
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
using namespace std;

int make_amc_gen_cost(TTree *t_gen, TH2D *h_2d, TH2D *h_2d_up, TH2D *h_2d_down, TH1D * h_weights, TH1D *h_cost, float m_low, float m_high, int year = 2016){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm, *inc1(0), *inc2(0);
    float gen_weight, m, cost, cost_st;
    int inc_id1, inc_id2;
    Bool_t sig_event(1);
    t_gen->SetBranchAddress("gen_p", &gen_lep_p);
    t_gen->SetBranchAddress("gen_m", &gen_lep_m);
    t_gen->SetBranchAddress("inc1", &inc1);
    t_gen->SetBranchAddress("inc2", &inc2);
    t_gen->SetBranchAddress("inc_id1", &inc_id1);
    t_gen->SetBranchAddress("inc_id2", &inc_id2);
    //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
    //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
    t_gen->SetBranchAddress("m", &m);
    t_gen->SetBranchAddress("cost", &cost);
    t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);
    t_gen->SetBranchAddress("sig_event", &sig_event);


    string pdfname;
    if(year == 2016)
        pdfname = string("NNPDF30_nlo_nf_5_pdfas");
    else
        pdfname = string("NNPDF31_nnlo_hessian_pdfas");
    const PDF* my_pdf = mkPDF(pdfname.c_str(),0);



    int nEvents=0;
    double max_weight = 1.;
    double my_max = 0.;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        //bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 30.;
        if(m >= m_low && m <= m_high && sig_event){
            nEvents++;
            cm = *gen_lep_p + *gen_lep_m;
            m = cm.M();
            float pt = cm.Pt();
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */
            float my_cost = cost_st;
            double E_BEAM = 7500.;
            double x1 = std::abs(inc1->E() / E_BEAM);
            double x2 = std::abs(inc2->E() / E_BEAM);
            double Q2 = x1*x2*(2*E_BEAM)*(2*E_BEAM);

            float scale = m*m + pt*pt;
            double x1p = my_pdf->xfxQ2(inc_id1, x1, scale)/x1;
            double x2p = my_pdf->xfxQ2(inc_id2, x2, scale)/x2;
            //double x1p = my_pdf->xfxQ2(inc_id1, x1, scale);
            //double x2p = my_pdf->xfxQ2(inc_id2, x2, scale);

            double pdf_reweight = 1./(x1p * x2p);
            pdf_reweight = min(pdf_reweight, max_weight);
            //pdf_reweight = 1.;
            double fill_weight = 1000. * gen_weight * pdf_reweight;

            h_2d->Fill(m, my_cost, fill_weight);

            if(abs(inc_id1) == 1 && abs(inc_id2) == 1 && inc_id1 * inc_id2 < 0){ //u ubar
                my_max = max(pdf_reweight, my_max);
                h_cost->Fill(my_cost, fill_weight);
                h_2d_up->Fill(m, my_cost, fill_weight);
                h_weights->Fill(pdf_reweight);
                //printf("p1 %f p2 %f Q %.1f M %.1f Pz %.1f \n", x1 * E_BEAM, x2 * E_BEAM, sqrt(Q2), m, cm.Pz());
                //printf("x1 %f x1p %f x2 %f x2p %f \n", x1, x1p, x2, x2p);
            }
            if(abs(inc_id1) == 2 && abs(inc_id2) == 2 && inc_id1 * inc_id2 < 0){ //d dbar
                my_max = max(pdf_reweight, my_max);
                h_cost->Fill(my_cost, fill_weight);
                h_2d_down->Fill(m, my_cost, fill_weight);
                h_weights->Fill(pdf_reweight);
                //printf("p1 %f p2 %f Q %.1f M %.1f Pz %.1f \n", x1 * E_BEAM, x2 * E_BEAM, sqrt(Q2), m, cm.Pz());
                //printf("x1 %f x1p %f x2 %f x2p %f \n", x1, x1p, x2, x2p);
            }


        }
    }
    printf("selected %i events \n", nEvents);
    printf("Max weight %f \n", my_max);

    return nEvents;

}

void normalize(TH2D *h){
    for(int i=1; i<=h->GetNbinsX(); i++){
        for(int j=1; j<=h->GetNbinsY(); j++){
            float xw = h->GetXaxis()->GetBinWidth(i);
            float yw = h->GetYaxis()->GetBinWidth(j);
            float content = h->GetBinContent(i,j);
            float new_content = content/(xw*yw);
            
            float err = h->GetBinError(i,j);
            //printf("i,j xw, yw %i %i %.2f %.2f \n", i,j, xw,yw);

            h->SetBinContent(i,j,content/(xw*yw));
            h->SetBinError(i,j,err/(xw*yw));
            if(content < 1e-6 || new_content < 0. || content * 1.5 < err){
                printf("WARNING bin %i %i (m = %.0f) has content %.3e +/- %.3e \n", i,j, h->GetXaxis()->GetBinCenter(i), content, err);
            }

        }
    }
}




void LQ_rw_denom(int year, bool write_out = false){

    char out_file[100], in_file[100];
    sprintf(out_file, "../../analyze/SFs/%i/LQ_rw.root", year);
    sprintf(in_file, "../../analyze/output_files/DY%i_gen_level_aug4.root", year % 2000);

    TFile * f_out;
    if(write_out)
        f_out = TFile::Open(out_file, "RECREATE");
    

    TFile *f_gen = TFile::Open(in_file);
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/A0_fits/";



    int n_cost_bins = 20;
    //float cost_bins[] = {-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.0};
    float cost_bins[] = {-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};
    //float cost_bins[] = {0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};

    int n_LQ_pt_bins = 1;
    //float LQ_pt_bins[] = {0., 20., 60., 100., 10000};
    float LQ_pt_bins[] = {0., 10000};

    int n_LQ_m_bins = 19;
    float m_LQ_bins[] = {350., 400., 450., 500., 550., 600.,  650., 700.,  800.,  900.,  1000.,  1100.,  
        1200.,  1300., 1400.,  1600.,  1800.,  2000.,  2500.,   3000., };

    TH2D *h_mu = new TH2D("h_mu", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);
    TH2D *h_el = new TH2D("h_el", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);

    TH2D *h_mu_up = new TH2D("h_mu_up", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);
    TH2D *h_el_up = new TH2D("h_el_up", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);
    TH2D *h_mu_down = new TH2D("h_mu_down", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);
    TH2D *h_el_down = new TH2D("h_el_down", "", n_LQ_m_bins, m_LQ_bins,   n_cost_bins, cost_bins);

    TH1D *h_1d = new TH1D("h1", "", n_cost_bins, cost_bins);
    TH1D *h_cost = new TH1D("h1", "", n_cost_bins, cost_bins);
    TH1D *h_w = new TH1D("h_w", "", 1000., 0., 1.);


    int nEvents = 0;
    float m_low = m_LQ_bins[0];
    float m_high = 14000.;
    //float m_low = 700.;
    //float m_high = 1000.;

    make_amc_gen_cost(t_gen_mu,  h_mu, h_mu_up, h_mu_down, h_w, h_cost, m_low,m_high, year);
    make_amc_gen_cost(t_gen_el, h_el, h_el_up, h_el_down, h_w, h_cost, m_low, m_high, year);
    
    /*
    TH1D *h_el_cost_pre = h_el->ProjectionY("h_el_cost_pre", 10,10);
    TH1D *h_mu_cost_pre = h_mu->ProjectionY("h_mu_cost_pre", 10,10);
    make_ratio_plot("LQ_cos_theta_pre.png", h_el_cost_pre, "El",h_mu_cost_pre, "Mu", "El/Mu", "cos(#theta_{*})", false, true);
    */

    normalize(h_mu);
    normalize(h_el);

    normalize(h_mu_up);
    normalize(h_el_up);
    normalize(h_mu_down);
    normalize(h_el_down);
    
    //h_el->Print("range");


    TH1D *h_el_cost = h_el_up->ProjectionY("h_el_cost", 10,10);
    TH1D *h_mu_cost = h_mu_up->ProjectionY("h_mu_cost", 10,10);
    make_ratio_plot("LQ_cos_theta_up.png", h_el_cost, "El",h_mu_cost, "Mu", "El/Mu", "cos(#theta_{*})", false, true);

    TH1D *h_el_m = h_el_up->ProjectionX("h_el_m");
    TH1D *h_mu_m = h_mu_up->ProjectionX("h_mu_m");
    make_ratio_plot("LQ_m_up.png", h_el_m, "El",h_mu_m, "Mu", "El/Mu", "M (GeV)", false, true);

    h_el_cost = h_el_down->ProjectionY("h_el_cost", 10,10);
    h_mu_cost = h_mu_down->ProjectionY("h_mu_cost", 10,10);
    make_ratio_plot("LQ_cos_theta_down.png", h_el_cost, "El",h_mu_cost, "Mu", "El/Mu", "cos(#theta_{*})", false, true);

    h_el_m = h_el_down->ProjectionX("h_el_m");
    h_mu_m = h_mu_down->ProjectionX("h_mu_m");
    make_ratio_plot("LQ_down.png", h_el_m, "El",h_mu_m, "Mu", "El/Mu", "M (GeV)", false, true);

    h_el_cost->SetLineColor(kBlue);
    h_mu_cost->SetLineColor(kRed);
    //h_el_cost->Draw();
    //h_mu_cost->Draw("same");
    TCanvas * c_w = new TCanvas("c_w", "", 800, 800);
    h_cost->Draw();
    printf("h_cost integral %f \n", h_cost->Integral());
    c_w->Print("cost.png");
    h_w->Draw();
    h_w->GetXaxis()->SetTitle("1./pdf weight");
    c_w->Print("weights.png");

    if(write_out){
        f_out->cd();
        h_mu->Write();
        h_el->Write();
        h_mu_up->Write();
        h_el_up->Write();
        h_mu_down->Write();
        h_el_down->Write();
        f_out->Print();
        f_out->Close();
    }

}

int main(int argc, char * argv[]){
    if (argc != 2){
        printf("Takes 1 argument (year) \n");
        exit(1);
    }
    int year = atoi(argv[1]);
    bool write_out = false;
    LQ_rw_denom(year, write_out);
}

