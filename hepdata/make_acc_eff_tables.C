#define STAND_ALONE
#include "../utils/root_files.h"
#include "../analyze/combine/LQ_TemplateUtils.h"
#include <iostream>

int fill_acc_eff(TTree *t_gen, TH3F *h_dy, TH3F *h_dy_m, TH3F *h_dy_pt, TH3F *h_dy_eta, TH3F *h_dy_all, bool flag_muon){
    
        TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
        float gen_weight, m, cost, cost_st;
        float evt_weight;
        Bool_t sig_event(1);
        t_gen->SetBranchAddress("gen_p", &gen_lep_p);
        t_gen->SetBranchAddress("gen_m", &gen_lep_m);
        //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
        //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
        t_gen->SetBranchAddress("m", &m);
        t_gen->SetBranchAddress("cost", &cost);
        t_gen->SetBranchAddress("cost_st", &cost_st);
        t_gen->SetBranchAddress("gen_weight", &gen_weight);
        t_gen->SetBranchAddress("sig_event", &sig_event);

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);

            cm = *gen_lep_p + *gen_lep_m;

            float gen_cost = cost_st;
            float pt = cm.Pt();
            float rap = abs(cm.Rapidity());
            float evt_weight = gen_weight * 100000.;
            float eta_p = abs(gen_lep_p->Eta());
            float eta_m = abs(gen_lep_m->Eta());
            float ang = gen_lep_p->Angle(gen_lep_m->Vect());

            h_dy->Fill(m, rap, gen_cost, evt_weight); 
	    
            if(m >= 500.) h_dy_m->Fill(m, rap, gen_cost, evt_weight);

            float leading = std::fmax(gen_lep_m->Pt(), gen_lep_p->Pt());
            float subleading = std::fmin(gen_lep_m->Pt(), gen_lep_p->Pt());

            if (leading > 40 and subleading > 15){
                h_dy_pt->Fill(m, rap, gen_cost, evt_weight); 
            }
            if (flag_muon){
                
                if( eta_p < 2.4 and eta_m < 2.4 and ang < (TMath::Pi() - 0.005))
                    h_dy_eta->Fill(m, rap, gen_cost, evt_weight);
            }
            else{
                if ((eta_p < 1.442 or (eta_p > 1.566 and eta_p < 2.5)) and (eta_m < 1.442 or (eta_m > 1.566 and eta_m < 2.5)))
                    h_dy_eta->Fill(m, rap, gen_cost, evt_weight);
            }
            bool pass = m >= 500 and leading > 40 and subleading > 15;
            if(flag_muon) pass = pass and eta_p < 2.4 and eta_m < 2.4 and ang < (TMath::Pi() - 0.005);
            else pass = pass and (eta_p < 1.442 or (eta_p > 1.566 and eta_p < 2.5)) and (eta_m < 1.442 or (eta_m > 1.566 and eta_m < 2.5));

            if(pass){
                h_dy_all->Fill(m, rap, gen_cost, evt_weight);
            }
            
        }
        return 1;

}


void make_acc_eff_tables(){

    gStyle->SetOptStat(0);
    gROOT->SetBatch(1);


    for(int year = 16; year <= 18; year++){
        char title[200];
        sprintf(title, "mumu%i_DY", year);
        auto h_mumu_dy = new TH3F(title, "mumu_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_mumu_dy->SetDirectory(0);
        sprintf(title, "mumu%i_DY_m", year);
        auto h_mumu_dy_m = new TH3F(title, "mumu_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_mumu_dy_m->SetDirectory(0);
        sprintf(title, "mumu%i_DY_pt", year);
        auto h_mumu_dy_pt = new TH3F(title, "mumu_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_mumu_dy_pt->SetDirectory(0);
        sprintf(title, "mumu%i_DY_eta", year);
        auto h_mumu_dy_eta = new TH3F(title, "mumu_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_mumu_dy_eta->SetDirectory(0);
        sprintf(title, "mumu%i_DY_all", year);
        auto h_mumu_dy_all = new TH3F(title, "mumu_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_mumu_dy_all->SetDirectory(0);

        sprintf(title, "elel%i_DY", year);
        auto h_elel_dy = new TH3F(title, "elel_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_elel_dy->SetDirectory(0);
	sprintf(title, "elel%i_DY_m", year);
        auto h_elel_dy_m = new TH3F(title, "elel_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_elel_dy_m->SetDirectory(0);
        sprintf(title, "elel%i_DY_pt", year);
        auto h_elel_dy_pt = new TH3F(title, "elel_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_elel_dy_pt->SetDirectory(0);
        sprintf(title, "elel%i_DY_eta", year);
        auto h_elel_dy_eta = new TH3F(title, "elel_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_elel_dy_eta->SetDirectory(0);
        sprintf(title, "elel%i_DY_all", year);
        auto h_elel_dy_all = new TH3F(title, "elel_dy",
             n_lq_m_bins, lq_m_bins, n_y_bins, y_bins, n_cost_bins, cost_bins);
        h_elel_dy_all->SetDirectory(0);
        
        char fin[180];
        sprintf(fin,"../analyze/output_files/DY%i_gen_level_april11.root", year);
        TFile *f1 = TFile::Open(fin, "READ");
        TTree *t_gen_mu = (TTree *) f1->Get("T_gen_mu");
        TTree *t_gen_el = (TTree *) f1->Get("T_gen_el");
        
        
        fill_acc_eff(t_gen_mu, h_mumu_dy, h_mumu_dy_m, h_mumu_dy_pt, h_mumu_dy_eta, h_mumu_dy_all, true);
        fill_acc_eff(t_gen_el, h_elel_dy, h_elel_dy_m, h_elel_dy_pt, h_elel_dy_eta, h_elel_dy_all, false);
        std::cout << "Filled 3D hists for elel and mumu for year " << year + 2000 << std::endl;

        auto h1_mumu_dy = convert3d(h_mumu_dy);
        auto h1_mumu_dy_m = convert3d(h_mumu_dy_m);
        auto h1_mumu_dy_pt = convert3d(h_mumu_dy_pt);
        auto h1_mumu_dy_eta = convert3d(h_mumu_dy_eta);
        auto h1_mumu_dy_all = convert3d(h_mumu_dy_all);

        auto h1_elel_dy = convert3d(h_elel_dy);
        auto h1_elel_dy_m = convert3d(h_elel_dy_m);
        auto h1_elel_dy_pt = convert3d(h_elel_dy_pt);
        auto h1_elel_dy_eta = convert3d(h_elel_dy_eta);
        auto h1_elel_dy_all = convert3d(h_elel_dy_all);
	std::cout << "converted 3D hists to 1D" << std::endl;
        
        char outfile1[200],outfile2[200] ;
        FILE *fout_mumu;
	    sprintf(outfile1,"mumu_acc_eff_y%i.txt",year);
	    fout_mumu = fopen(outfile1, "w");
        FILE *fout_elel;
	    sprintf(outfile2,"elel_acc_eff_y%i.txt",year);
	    fout_elel = fopen(outfile2, "w");
        
        float pass_m, pass_pt, pass_eta, pass_all;
        for(int i = 1; i <= h1_mumu_dy->GetNbinsX(); i++){
		
	    pass_m = h1_mumu_dy_m->GetBinContent(i)/h1_mumu_dy->GetBinContent(i);
	    pass_pt = h1_mumu_dy_pt->GetBinContent(i)/h1_mumu_dy->GetBinContent(i);
            pass_eta = h1_mumu_dy_eta->GetBinContent(i)/h1_mumu_dy->GetBinContent(i);
	    pass_all = h1_mumu_dy_all->GetBinContent(i)/h1_mumu_dy->GetBinContent(i);
            fprintf(fout_mumu,"%i %.4f %.4f %.4f %.4f\n", i, pass_m, pass_pt, pass_eta, pass_all);

	    pass_m = h1_elel_dy_m->GetBinContent(i)/h1_elel_dy->GetBinContent(i);
            pass_pt = h1_elel_dy_pt->GetBinContent(i)/h1_elel_dy->GetBinContent(i);
            pass_eta = h1_elel_dy_eta->GetBinContent(i)/h1_elel_dy->GetBinContent(i);
            pass_all = h1_elel_dy_all->GetBinContent(i)/h1_elel_dy->GetBinContent(i);    
            fprintf(fout_elel,"%i %.4f %.4f %.4f %.4f\n", i, pass_m, pass_pt, pass_eta, pass_all);
        }
	std::cout << "Saved tables to files" << std::endl;

    }
    }
    



    
