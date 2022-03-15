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
#include "../../utils/root_files.h"
//#include "../../utils/HistUtils.C"
#include "../../utils/ScaleFactors.C"
#include "../../utils/PlotUtils.C"
//#include "../../utils/LQ_TemplateMaker_systematics.C"
#include "LQ_TemplateUtils.h"

#include "THStack.h"



void LQ_make_gen_templates(){

    gStyle->SetOptStat(0);
    gROOT->SetBatch(1);

        /*
     1 Madgraph5_aMC@NLO R2 LQ -> ee with RL couplings
  2
  3 LQ Mass  Y_ue  cross section(pb)  N_events  int Lum(fb-1)
  4 1 TeV    1.0   0.1866+-0.00072    41859       224.3
  5 2 TeV    1.0   0.08006+-0.00027   48939       611.3
  6 2 TeV    0.0   0.06132+-0.00020   50617       825.5 (SM events)
  */

//    for(int year=2016;year<=2018;year++){
        bool do_ptrw = false;
        //float m_LQ = 1000.;
        float m_LQ = 2000.;
        int year = 2017;
        char fout_name[200];
        sprintf(fout_name,"combine/templates/LQm%i_gen_templates%i_020222.root",int(m_LQ),year%2000);
        //sprintf(fout_name,"combine/templates/LQm%i_SM_gen_templates%i_020222.root",int(m_LQ),year%2000);
        string fout_n = string(fout_name, 200);

        char genfile_name[200];
        sprintf(genfile_name,"../analyze/output_files/DY%i_gen_level_aug4.root",year%2000);
        string genfile_n = string(genfile_name,200);
        TFile *f_gen = TFile::Open(genfile_n.c_str());
	char lhe_name[200];
	sprintf(lhe_name,"../analyze/root_files/LQ_m%i_test2_020122.root",int(m_LQ));
	string lhe_n = string(lhe_name,200);
        TFile *f_gen_data = TFile::Open(lhe_n.c_str());
        TFile *f_gen_data_SM = TFile::Open("../analyze/root_files/LQ_SM_test2_020122.root");
        //gROOT->SetBatch(1);

        //TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
        TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");
        TTree *t_gen_data = (TTree *) f_gen_data->Get("T_lhe");
        TTree *t_gen_data_SM = (TTree *) f_gen_data_SM->Get("T_lhe");
        //calculate the total gen_weights to scale the data temps 
        float gen_weight, sum_weights;
        //float xsec = 0.1866, xsec_SM = 0.06132;
        float xsec = 0.08, xsec_SM = 0.06132;
	int nevents = 48939, nevents_SM = 50617;
        TFile * fout = TFile::Open(fout_n.c_str(), "RECREATE");

        char dirname[40], title[300];

        string sys = "";
        printf("Starting year %i\n",year);
        
        sprintf(title, "ee%i_data_obs", year %2000);
        TH3F* h_data = new TH3F(title, "Data template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_data->SetDirectory(0);


        sprintf(title, "ee%i_data_obs", year %2000);
        TH3F* h_data_SM = new TH3F(title, "Data template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_data_SM->SetDirectory(0);

        sprintf(title, "ee%i_sym", year %2000);
        TH3F* h_sym = new TH3F(title, "sym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_sym->SetDirectory(0);

        sprintf(title, "ee%i_asym", year %2000);
        TH3F* h_asym = new TH3F(title, "asym template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_asym->SetDirectory(0);

        sprintf(title, "ee%i_alpha", year %2000);
        TH3F* h_alpha = new TH3F(title, "alpha template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_alpha->SetDirectory(0);
        sprintf(title, "ee%i_LQpure_u", year %2000);
        TH3F* h_LQpure_u = new TH3F(title, "LQpure_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQpure_u->SetDirectory(0);

        sprintf(title, "ee%i_LQpure_d", year %2000);
        TH3F* h_LQpure_d = new TH3F(title, "LQpure_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQpure_d->SetDirectory(0);

        sprintf(title, "ee%i_LQint_u", year %2000);
        TH3F* h_LQint_u = new TH3F(title, "LQint_u template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQint_u->SetDirectory(0);

        sprintf(title, "ee%i_LQint_d", year %2000);
        TH3F* h_LQint_d = new TH3F(title, "LQint_d template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_LQint_d->SetDirectory(0);

        sprintf(title, "ee%i_raw", year %2000);
        TH3F* h_raw = new TH3F(title, "raw template for gen level",
            n_lq_m_bins, lq_m_bins,n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_raw->SetDirectory(0);

        TH1F *h1_raw, *h1_data, *h1_data_SM, *h1_pl, *h1_mn, *h1_asym, *h1_sym, *h1_alpha, *h1_LQpure_u, *h1_LQint_u,*h1_LQpure_d, *h1_LQint_d, *h1_LQpure_test_u, *h1_LQint_test_u,*h1_LQpure_test_d, *h1_LQint_test_d;


        
           // printf("\n \n Start making gen level templates for LQ+DY\n");


        int nEvents = 0;
        int nEvents_data = 0;

    //nEvents += make_gen_temps(t_gen_mu, h_uncut, h_raw, h_sym, h_asym, h_alpha, m_low, m_high, do_ptrw, year, sys);
        sum_weights = make_gen_temps(t_gen_el, h_raw, h_sym, h_asym, h_alpha, h_LQpure_u, h_LQpure_d, h_LQint_u, h_LQint_d,  m_LQ, true, do_ptrw, year, sys);
    	//printf("Finished make_gen_temps, nEvents = %i\n",nEvents);
        nEvents_data += make_gen_data_temps(t_gen_data, h_data, xsec, nevents, year, sum_weights);
        nEvents_data+= make_gen_data_temps(t_gen_data_SM, h_data_SM, xsec_SM, nevents_SM, year, sum_weights);

        //ratio of NLO SM to LO SM
        // k factor extracted for all m, rap
        float k_factor = h_raw->Integral() / h_data_SM->Integral();

        printf("k-factor is %.3f \n", k_factor);
       // h_data->Scale(k_factor);
       // h_data_SM->Scale(k_factor);

        // k factor extracted for 3*3 mass*rap bins
       
        for(int i = 1; i <= h_raw->GetNbinsX(); i++){
            for( int j = 1; j <= h_raw->GetNbinsY(); j++){

                float sum_bins_LO = 0., sum_bins_NLO = 0.;
                for(int k = 1; k <= h_raw->GetNbinsZ(); k++){

                    sum_bins_LO+=h_data_SM->GetBinContent(i,j,k);
                    sum_bins_NLO+=h_raw->GetBinContent(i,j,k);
                }
                //extract k_factor for each mass and rap bin
                k_factor = sum_bins_NLO/sum_bins_LO;
                printf("k-factor %i = %f\n",i*n_lq_m_bins+j,k_factor);
                for(int k = 1; k <= h_data->GetNbinsZ(); k++){
                    float SM_content = h_data_SM->GetBinContent(i,j,k);
                    float SMLQ_content = h_data->GetBinContent(i,j,k);
                    h_data->SetBinContent(i,j,k, SMLQ_content*k_factor);
                    h_data_SM->SetBinContent(i,j,k,SM_content*k_factor);
                }
            }
        }

        printf("h_data integral %.2f, h_raw interal %.2f \n", h_data->Integral(), h_raw->Integral());
        printf("h_LQpure_u %.2f h_LQint_u %.2f \n", h_LQpure_u->Integral(), h_LQint_u->Integral());

        

        h_sym->Scale(0.5);
        h_asym->Scale(0.5);
        h_alpha->Scale(0.5);
        h_LQpure_u->Scale(0.5);
        h_LQpure_d->Scale(0.5);
        h_LQint_u->Scale(0.5);
        h_LQint_d->Scale(0.5);
        printf("Starting convert3d\n");

        h1_data = convert3d(h_data);
        h1_data_SM = convert3d(h_data_SM);
        h1_raw = convert3d(h_raw);
        h1_sym = convert3d(h_sym);
        h1_asym = convert3d(h_asym);
        h1_alpha = convert3d(h_alpha);
        h1_LQpure_u = convert3d(h_LQpure_u);
        h1_LQint_u = convert3d(h_LQint_u);
        h1_LQpure_d = convert3d(h_LQpure_d);
        h1_LQint_d = convert3d(h_LQint_d);
        //==============================================================================================================================

        TF1 *lq_check = new TF1("LQ shape",check_analytical(cost, m_LQ, 1.0, 500.*500. , Q_u, caq_u, cvq_u), -1.,1.);

        TCanvas *c_mumu6 = new TCanvas("c_mumu6", "Histograms", 200, 10, 900, 700);
        //h1_LQint_test_u->SetLineColor(kBlue);
        //h1_LQint_u->SetLineColor(kRed);
        //h1_LQint_test_u->SetLineWidth(2);
        //h1_LQint_u->SetLineWidth(2); 
        lq_check->SetTitle("LQpure_u+LQint_u, m_ll = 500, y_ue=1.0");
        lq_check->Draw();
        //h1_LQint_u->Draw("hist same ");
        //TLegend *leg6 = new TLegend(0.75, 0.75, 0.9, 0.9);
        //leg6->AddEntry(h1_LQint_u, "asym denom", "l");
        //leg6->AddEntry(h1_LQint_test_u, "sym denom", "l");
        //leg6->Draw();
        sprintf(title, "../generator_stuff/plots/LQ_shapecheck.png");
        c_mumu6->Print(title);
        delete c_mumu6;


        // testing the symmetric denominator
        /*
	bool only_sym = true;

        h_LQpure_u->Reset();
        h_LQpure_d->Reset();
        h_LQint_u->Reset();
        h_LQint_d->Reset();


        sum_weights = make_gen_temps(t_gen_el, h_raw, h_sym, h_asym, h_alpha, h_LQpure_u, h_LQpure_d, h_LQint_u, h_LQint_d,  m_LQ, only_sym, do_ptrw, year, sys);

        
        h_LQpure_u->Scale(0.5);
        h_LQpure_d->Scale(0.5);
        h_LQint_u->Scale(0.5);
        h_LQint_d->Scale(0.5);   
        

        h1_LQpure_test_u = convert3d(h_LQpure_u);
        h1_LQint_test_u = convert3d(h_LQint_u);
        h1_LQpure_test_d = convert3d(h_LQpure_d);
        h1_LQint_test_d = convert3d(h_LQint_d);

        delete h_sym,h_asym,h_alpha,h_LQpure_u,h_LQpure_d,h_LQint_u,h_LQint_d,h_data,h_data_SM;

        TCanvas *c_mumu5 = new TCanvas("c_mumu5", "Histograms", 200, 10, 900, 700);
        h1_LQpure_test_u->SetLineColor(kBlue);
        h1_LQpure_u->SetLineColor(kRed);
        h1_LQpure_test_u->SetLineWidth(2);
        h1_LQpure_u->SetLineWidth(2); 
        h1_LQpure_test_u->SetTitle("LQpure_u");
        h1_LQpure_test_u->Draw("hist");
        h1_LQpure_u->Draw("hist same ");
        TLegend *leg5 = new TLegend(0.75, 0.75, 0.9, 0.9);
        leg5->AddEntry(h1_LQpure_u, "asym denom", "l");
        leg5->AddEntry(h1_LQpure_test_u, "sym denom", "l");
        leg5->Draw();
        sprintf(title, "../generator_stuff/plots/LQpure_u_symd_vs_asymd.png");
        c_mumu5->Print(title);
        delete c_mumu5;

        


*/        
        //===================================================================================================================================
        printf("Starting make_pl_mn\n");
            //h1_sym->Print("range");
    	//h1_asym->Print("range");

        // n_y_bins -= 1;
        int n_1d_bins = get_n_1d_bins(n_y_bins, n_cost_bins);
        sprintf(title, "ee%i_fpl", year %2000);
        h1_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_pl->SetDirectory(0);
        sprintf(title, "ee%i_fmn", year %2000);
        h1_mn = new TH1F(title, "Minus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mn->SetDirectory(0);


        make_pl_mn_templates(h1_sym, h1_asym, h1_pl, h1_mn); 
        printf("Finished make_pl_mn\n");
      
        
        

        fout->cd();
        snprintf(dirname, 10, "LQ");
        gDirectory->mkdir(dirname);
        gDirectory->cd(dirname);

            //h_uncut->Write();
       // h1_data->Scale(.5);
        h1_data->Write();
        h1_pl->Write();
        h1_mn->Write();
        h1_alpha->Write();
        h1_sym->Write();
        h1_asym->Write();
        h1_LQpure_u->Write();
        h1_LQpure_d->Write();
        h1_LQint_u->Write();
        h1_LQint_d->Write();

        fout->Close();

        float a0 = 0.05, afb = 0.62;
        float alph = 2.*a0/(2.- a0);//alph=2/39
        float norm = 3./(4.*(2.+alph));//norm=0.365625
        TH1F *h1_total = (TH1F *) h1_alpha->Clone("h1_total");
        h1_total->SetDirectory(0);
        TH1F *h1_total_SM_NLO = (TH1F *) h1_alpha->Clone("h1_total_SM_NLO");
        h1_total_SM_NLO->SetDirectory(0);
        h1_total->Scale(alph*norm);//0.01875
        h1_total_SM_NLO->Scale(alph*norm);//0.01875
        h1_pl->Scale(norm+afb);//0.965625
        h1_mn->Scale(norm-afb);//-.234375


        h1_total->Add(h1_pl);
        h1_total->Add(h1_mn);
        h1_total->Add(h1_LQpure_u);
        h1_total->Add(h1_LQint_u);
        h1_total_SM_NLO->Add(h1_pl);
        h1_total_SM_NLO->Add(h1_mn);
       // h1_total->Write();
	
        TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
        h1_data->SetLineColor(kBlue);
        h1_total->SetLineColor(kRed);
        h1_data->SetLineWidth(2);
        h1_total->SetLineWidth(2); 
        h1_total->SetTitle("Fake data (SM+LQ) vs DY+LQ templates (AFB=0.6,A0=0.05,y_ue=1.0,mLQ=1000)");
        h1_data->Draw("hist");
        h1_total->Draw("hist same ");
        TLegend *leg1 = new TLegend(0.75, 0.75, 0.9, 0.9);
        leg1->AddEntry(h1_total, "DY+LQ_ue templates", "l");
        leg1->AddEntry(h1_data, "Fake data (y_ue=1.0)", "l");
        leg1->Draw();
        sprintf(title, "../generator_stuff/plots/data_vs_total_%i.png", year %2000);
        c_mumu1->Print(title);
        delete c_mumu1;

        TCanvas *c_mumu3 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
        h1_data->SetLineColor(kBlue);
        h1_data_SM->SetLineColor(kRed);
        h1_data->SetLineWidth(2);
        h1_data_SM->SetLineWidth(2); 
        h1_data_SM->SetTitle("Fake data (SM+LQ) vs Fake data (SM) (AFB=0.6,A0=0.05,y_ue=1.0,mLQ=1000)");
        h1_data_SM->Draw("hist");
        h1_data->Draw("hist same ");
        TLegend *leg3 = new TLegend(0.75, 0.75, 0.9, 0.9);
        leg3->AddEntry(h1_data_SM, "Fake data only SM", "l");
        leg3->AddEntry(h1_data, "Fake data (y_ue=1.0)", "l");
        leg3->Draw();
        sprintf(title, "../generator_stuff/plots/data_vs_noLQ_%i.png", year %2000);
        c_mumu3->Print(title);
        delete c_mumu3;

        //checking NLO vs LO


        TCanvas *c_mumu4 = new TCanvas("c_mumu4", "Histograms", 200, 10, 900, 700);
        h1_data_SM->SetLineColor(kBlue);
        h1_total_SM_NLO->SetLineColor(kRed);
        h1_data_SM->SetLineWidth(2);
        h1_total_SM_NLO->SetLineWidth(2); 
        h1_data_SM->SetTitle("Fake SM@NLO vs Fake SM (AFB=0.6,A0=0.05)");
        //h1_total_SM_NLO->Draw("hist");
	h1_data_SM->Draw("hist ");
        h1_total_SM_NLO->Draw("hist same ");
        TLegend *leg4 = new TLegend(0.75, 0.75, 0.9, 0.9);
        //leg2->AddEntry(h1_LQpure_u, "LQ_ue templates", "l");
        leg4->AddEntry(h1_data_SM, "Fake SM ", "l");
        leg4->AddEntry(h1_total_SM_NLO, "Fake SM@NLO ", "l");
        leg4->Draw();
        sprintf(title, "../generator_stuff/plots/SM_vs_SM@NLO_%i.png", year %2000);
        c_mumu4->Print(title);
        delete c_mumu4;

        //checking LQ temps vs fake data minus SM

        h1_data_SM->Scale(-1.);
        h1_data->Add(h1_data_SM);
        //h1_LQint_u->Scale(-1.);
        h1_LQpure_u->Add(h1_LQint_u);

        TCanvas *c_mumu2 = new TCanvas("c_mumu2", "Histograms", 200, 10, 900, 700);
        h1_data->SetLineColor(kBlue);
        h1_LQpure_u->SetLineColor(kRed);
        h1_data->SetLineWidth(2);
        h1_LQpure_u->SetLineWidth(2); 
        h1_data->SetTitle("Fake data minus Fake SM (AFB=0.6,A0=0.05,y_ue=1.0,mLQ=1000)");
        h1_data->Draw("hist");
        //h1_LQpure_u->Draw("hist same ");
        TLegend *leg2 = new TLegend(0.75, 0.75, 0.9, 0.9);
        //leg2->AddEntry(h1_LQpure_u, "LQ_ue templates", "l");
        leg2->AddEntry(h1_data, "Fake data minus SM (y_ue=1.0)", "l");
        leg2->Draw();
        sprintf(title, "../generator_stuff/plots/data_minus_SM_%i.png", year %2000);
        c_mumu2->Print(title);
        delete c_mumu2;
/*
        THStack *hs = new THStack("hs","");
        THStack->SetDirectory(0);
        h1_pl->SetLineColor(kBlue);
        h1_pl->SetFillColor(kBlue);
        *hs->Add(h1_pl);
        h1_mn->SetLineColor(kRed);
        h1_mn->SetFillColor(kRed);
        *hs->Add(h1_mn);
        h1_alpha->SetLineColor(kGreen);
        h1_alpha->SetFillColor(kGreen);
        *hs->Add(h1_alpha);

         
*/
        

            //h_uncut->Reset();
            //h_raw->Reset();
        h1_data->Reset();
        h1_data_SM->Reset();
        h1_pl->Reset();
        h1_mn->Reset();
        h1_alpha->Reset();
        h1_asym->Reset();
        h1_sym->Reset();
        h1_LQpure_u->Reset();
        h1_LQpure_d->Reset();
        h1_LQint_u->Reset();
        h1_LQint_d->Reset();
        h1_total->Reset();
        



        printf("Templates written to %s \n", fout_n.c_str());

    }
//}
