#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/LQ_TemplateUtils.h"
#include <iostream>





void LQ_draw_bkg_templates(){


	bool scramble_data =false ;
    bool fake_data =true; //use mc instead of data
    bool use_xF = false;
    const string sys_label = "";
    char *plot_dir = "AN_plots/Bkg_templates/";
    int i_start = 1; 
    int i_end = 4;
    float x_start = 0.7;
    float x_end = 0.9;
    float y_start = 0.75;
    float y_end = 0.9;
   
   // std::cout << "enter m_LQ:";        
   // std::cin >> m_LQ; 
    gStyle->SetOptStat(0);
    gROOT->SetBatch(1);


    for(int year = 2016; year<=2018; year++){

    printf("Initializing files \n");
    init(year);
    //init_emu(year);
    printf("Setting up SFs... ");
    setup_all_SFs(year);

    TH1F *h1_elel_asym, *h1_elel_sym; 
    TH1F *h1_mumu_asym, *h1_mumu_sym; 
    TH1F *h1_elel_pl, *h1_elel_mn, *h1_elel_alpha;
    TH1F *h1_mumu_pl, *h1_mumu_mn, *h1_mumu_alpha;
    TH1F *h1_elel_db, *h1_elel_top,  *h1_elel_tautau, *h1_elel_data, *h1_elel_mc, *h1_elel_qcd, *h1_elel_gam;
	TH1F *h1_mumu_db, *h1_mumu_top, *h1_mumu_tautau, *h1_mumu_data, *h1_mumu_mc, *h1_mumu_qcd, *h1_mumu_gam;

	int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }	
	char title[100];
    //titles taken care of in conversion

    sprintf(title, "ee%i_qcd%s", year %2000, sys_label.c_str());
    TH3F* h_elel_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_elel_qcd->SetDirectory(0);
    sprintf(title, "mumu%i_qcd%s", year %2000, sys_label.c_str());
    TH3F* h_mumu_qcd = new TH3F(title, "Fakes template",
             n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_mumu_qcd->SetDirectory(0);
    bool incl_ss = true;
    bool ss_binning = false;
    float elel_sign_scaling, elel_err, mumu_sign_scaling, mumu_err;
    printf("making ElEl fakes template \n");
    gen_fakes_template(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, h_elel_qcd, year, 
                FLAG_ELECTRONS, incl_ss, ss_binning, use_xF, sys_label);
    printf("making MuMu fakes template \n");
    gen_fakes_template(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, h_mumu_qcd, year, FLAG_MUONS, 
                incl_ss, ss_binning, use_xF, sys_label);
   


    printf("Integral of QCD templates are %.2f %.2f \n", h_elel_qcd->Integral(), h_mumu_qcd->Integral());

    symmetrize3d(h_mumu_qcd);
    symmetrize3d(h_elel_qcd);

    h1_mumu_qcd = convert3d(h_mumu_qcd);
    h1_elel_qcd = convert3d(h_elel_qcd);

    printf("Made qcd templates \n");
    delete h_elel_qcd, h_mumu_qcd;

    bool ss= false;
    
   sprintf(title, "mumu%i_top%s", year %2000, sys_label.c_str());
        auto h_mumu_top = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins,n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_top->SetDirectory(0);
        sprintf(title, "mumu%i_db%s", year %2000, sys_label.c_str());
        auto h_mumu_db = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_db->SetDirectory(0);
        sprintf(title, "mumu%i_tautau%s", year %2000, sys_label.c_str());
        auto h_mumu_tautau = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_tautau->SetDirectory(0);
        sprintf(title, "mumu%i_gam%s", year %2000, sys_label.c_str());
        auto h_mumu_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_gam->SetDirectory(0);

        TTree *mumu_ts[2] = {t_mumu_ttbar, t_mumu_wt};
        printf("Making mumu back \n");
        bool emu_costrw = true;
        gen_combined_background_template(2, mumu_ts, h_mumu_top, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);


        mumu_ts[0] = t_mumu_diboson;
        gen_combined_background_template(1, mumu_ts, h_mumu_db, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);
        emu_costrw = false;

        mumu_ts[0] = t_mumu_tautau;
        printf("Making mumu tautau \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_tautau, year, FLAG_MUONS,  ss, use_xF,  emu_costrw, sys_label);

        mumu_ts[0] = t_mumu_gamgam;
        printf("Making mumu gamgam \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_gam, year, FLAG_MUONS,  ss, use_xF, emu_costrw, sys_label);

         symmetrize3d(h_mumu_top);
        //symmetrize3d(h_mumu_db);

        
        h1_mumu_top = convert3d(h_mumu_top);
        h1_mumu_db = convert3d(h_mumu_db);
        h1_mumu_tautau = convert3d(h_mumu_tautau);
        h1_mumu_gam = convert3d(h_mumu_gam);
        delete h_mumu_top, h_mumu_db, h_mumu_tautau, h_mumu_gam;


        sprintf(title, "ee%i_top%s", year %2000, sys_label.c_str());
        auto h_elel_top = new TH3F(title, "Combined background template",
               n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_top->SetDirectory(0);
        sprintf(title, "ee%i_db%s", year %2000, sys_label.c_str());
        auto h_elel_db = new TH3F(title, "Combined background template",
                n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_db->SetDirectory(0);
        sprintf(title, "ee%i_tautau%s", year %2000, sys_label.c_str());
        auto h_elel_tautau = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_tautau->SetDirectory(0);
        sprintf(title, "ee%i_gam%s", year %2000, sys_label.c_str());
        auto h_elel_gam = new TH3F(title, "Combined background template",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_gam->SetDirectory(0);


        TTree *elel_ts[2] = {t_elel_ttbar, t_elel_wt};
         emu_costrw = true;
        gen_combined_background_template(2, elel_ts, h_elel_top, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);

        elel_ts[0] = t_elel_diboson;
        gen_combined_background_template(1, elel_ts, h_elel_db, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);
        emu_costrw = false;

        elel_ts[0] = t_elel_tautau;
        gen_combined_background_template(1, elel_ts, h_elel_tautau, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);

        elel_ts[0] = t_elel_gamgam;
        printf("Making ee gamgam \n");
        gen_combined_background_template(1, elel_ts, h_elel_gam, year, FLAG_ELECTRONS, ss, use_xF, emu_costrw, sys_label);

        symmetrize3d(h_elel_top);
        //symmetrize3d(h_elel_db);
        
       
        h1_elel_top = convert3d(h_elel_top);
        h1_elel_db = convert3d(h_elel_db);
        h1_elel_tautau = convert3d(h_elel_tautau);
        h1_elel_gam = convert3d(h_elel_gam);

        delete h_elel_top, h_elel_db, h_elel_tautau, h_elel_gam;

        h1_mumu_qcd->SetLineWidth(2);
        h1_mumu_top->SetLineWidth(2);
        h1_mumu_db->SetLineWidth(2);
        h1_mumu_tautau->SetLineWidth(2);
        h1_mumu_gam->SetLineWidth(2);

        h1_mumu_qcd->SetLineColor(kRed);
        h1_mumu_top->SetLineColor(kBlue);
        h1_mumu_db->SetLineColor(kGreen);
        h1_mumu_tautau->SetLineColor(kOrange);
        h1_mumu_gam->SetLineColor(kYellow);

        h1_elel_qcd->SetLineWidth(2);
        h1_elel_top->SetLineWidth(2);
        h1_elel_db->SetLineWidth(2);
        h1_elel_tautau->SetLineWidth(2);
        h1_elel_gam->SetLineWidth(2);

        h1_elel_qcd->SetLineColor(kRed);
        h1_elel_top->SetLineColor(kBlue);
        h1_elel_db->SetLineColor(kGreen);
        h1_elel_tautau->SetLineColor(kOrange);
        h1_elel_gam->SetLineColor(kYellow);

        char mu_title[100], el_title[100];
        char mu_fname1[100], el_fname1[100];
        sprintf(mu_fname1, "%s/Mu%i_bkgs.png", plot_dir, year%2000);
        sprintf(el_fname1, "%s/El%i_bkgs.png", plot_dir, year%2000);
            

        sprintf(mu_title, "Channel : Muons %i, All backgrounds", year);
        TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
         h1_mumu_top->SetTitle(mu_title);
         h1_mumu_top->Draw("hist");
         h1_mumu_qcd->Draw("hist same");
         h1_mumu_db->Draw("hist same");
         h1_mumu_tautau->Draw("hist same");
         h1_mumu_gam->Draw("hist same");

          TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_mumu_top,"ttbar+tW","l");
            leg1->AddEntry(h1_mumu_qcd,"W+jets+QCD","l");
            leg1->AddEntry(h1_mumu_db,"WW+ZZ+WZ","l");
            leg1->AddEntry(h1_mumu_tautau,"DY #tau#tau","l");
            leg1->AddEntry(h1_mumu_gam,"#gamma#gamma -> ll","l");

            leg1->Draw();
            
            c_mumu1->Print(mu_fname1);
            delete c_mumu1;

         sprintf(el_title, "Channel : Electrons %i, All backgrounds", year);
        TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
         h1_elel_top->SetTitle(el_title);
         h1_elel_top->Draw("hist");
         h1_elel_qcd->Draw("hist same");
         h1_elel_db->Draw("hist same");
         h1_elel_tautau->Draw("hist same");
         h1_elel_gam->Draw("hist same");

         TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_elel_top,"ttbar+tW","l");
            leg2->AddEntry(h1_elel_qcd,"W+jets+QCD","l");
            leg2->AddEntry(h1_elel_db,"WW+ZZ+WZ","l");
            leg2->AddEntry(h1_elel_tautau,"DY #tau#tau","l");
            leg2->AddEntry(h1_elel_gam,"#gamma#gamma -> ll","l");

            leg2->Draw();
            
            c_elel1->Print(el_fname1);
            delete c_elel1;



        printf("making muon mc templates \n");
        sprintf(title, "mumu%i_sym%s", year %2000, sys_label.c_str());
        auto h_mumu_sym = new TH3F(title, "Symmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_sym->SetDirectory(0);
        sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
        auto h_mumu_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_alpha->SetDirectory(0);
        sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
        auto h_mumu_asym = new TH3F(title, "Asymmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_mumu_asym->SetDirectory(0);
		

        printf("Making mumu mc \n");
        //gen_mc_template includes m_LQ
        gen_mc_SM_template(t_mumu_mc,  h_mumu_sym, h_mumu_asym, h_mumu_alpha, year, FLAG_MUONS, use_xF, sys_label );
        
        h1_mumu_sym = convert3d(h_mumu_sym);
        h1_mumu_asym = convert3d(h_mumu_asym);
        h1_mumu_alpha = convert3d(h_mumu_alpha);

        printf("making electron mc templates \n");
        sprintf(title, "ee%i_sym%s", year %2000, sys_label.c_str());
        auto h_elel_sym = new TH3F(title, "Symmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_sym->SetDirectory(0);
        sprintf(title, "ee%i_alpha%s", year %2000, sys_label.c_str());
        auto h_elel_alpha = new TH3F(title, "Gauge boson polarization template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_alpha->SetDirectory(0);
        sprintf(title, "ee%i_asym%s", year %2000, sys_label.c_str());
        auto h_elel_asym = new TH3F(title, "Asymmetric template of mc",
                 n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        h_elel_asym->SetDirectory(0);

        printf("starting elel dy \n");
        //gen_mc_template includes m_LQ
        gen_mc_SM_template(t_elel_mc,  h_elel_sym, h_elel_asym, h_elel_alpha, year, FLAG_ELECTRONS, use_xF, sys_label );
        
        h1_elel_sym = convert3d(h_elel_sym);
        h1_elel_asym = convert3d(h_elel_asym);
        h1_elel_alpha = convert3d(h_elel_alpha);

        int n_1d_bins = get_n_1d_bins(n_var1_bins, n_cost_bins);
  

        //char title[100];
        sprintf(title, "mumu%i_fpl%s", year%2000, sys_label.c_str());
        h1_mumu_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mumu_pl->SetDirectory(0);
        sprintf(title, "mumu%i_fmn%s", year%2000, sys_label.c_str());
        h1_mumu_mn = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_mumu_mn->SetDirectory(0);
        make_pl_mn_templates(h1_mumu_sym, h1_mumu_asym, h1_mumu_pl, h1_mumu_mn);


        h1_mumu_sym->Reset();
        h1_mumu_asym->Reset();



        //char title[100];
        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_elel_pl->SetDirectory(0);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn = new TH1F(title, "Plus template of DY", n_1d_bins, 0, n_1d_bins);
        h1_elel_mn->SetDirectory(0);
        make_pl_mn_templates(h1_elel_sym, h1_elel_asym, h1_elel_pl, h1_elel_mn);

        sprintf(title, "ee%i_fpl%s", year%2000, sys_label.c_str());
        h1_elel_pl->SetName(title);
        sprintf(title, "ee%i_fmn%s", year%2000, sys_label.c_str());
        h1_elel_mn->SetName(title);

        h1_elel_sym->Reset();
        h1_elel_asym->Reset();

        h1_mumu_pl->SetLineWidth(2);
        h1_mumu_mn->SetLineWidth(2);
        h1_mumu_alpha->SetLineWidth(2);
        //h1_mumu_tautau->SetLineWidth(2);
        //h1_mumu_gam->SetLineWidth(2);

        h1_mumu_pl->SetLineColor(kRed);
        h1_mumu_mn->SetLineColor(kBlue);
        h1_mumu_alpha->SetLineColor(kGreen);
        //h1_mumu_tautau->SetLineColor(kOrange);
        //h1_mumu_gam->SetLineColor(kYellow);

        h1_elel_pl->SetLineWidth(2);
        h1_elel_mn->SetLineWidth(2);
        h1_elel_alpha->SetLineWidth(2);
        //h1_elel_tautau->SetLineWidth(2);
        //h1_elel_gam->SetLineWidth(2);

        h1_elel_pl->SetLineColor(kRed);
        h1_elel_mn->SetLineColor(kBlue);
        h1_elel_alpha->SetLineColor(kGreen);

        char mu_title2[100], el_title2[100];
        char mu_fname2[100], el_fname2[100];
        sprintf(mu_fname2, "%s/Mu%i_DYbkgs.png", plot_dir, year%2000);
        sprintf(el_fname2, "%s/El%i_DYbkgs.png", plot_dir, year%2000);

        sprintf(el_title2, "Channel : Electrons %i; DY background",year);
        TCanvas *c_elel2 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
         h1_elel_alpha->SetTitle(el_title2);
         h1_elel_alpha->Draw("hist");
         h1_elel_pl->Draw("hist same");
         h1_elel_mn->Draw("hist same");
         

         TLegend *leg3 = new TLegend(x_start, y_start, x_end, y_end);
            leg3->AddEntry(h1_elel_pl,"f_{+} template","l");
            leg3->AddEntry(h1_elel_mn,"f_{-} template","l");
            leg3->AddEntry(h1_elel_alpha,"f_{#alpha} template","l");

            leg3->Draw();
            
            c_elel2->Print(el_fname2);
            delete c_elel2;

            sprintf(mu_title2, "Channel : Muons %i; DY background",year);
        TCanvas *c_mumu2 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
         h1_mumu_alpha->SetTitle(mu_title2);
         h1_mumu_alpha->Draw("hist");
         h1_mumu_pl->Draw("hist same");
         h1_mumu_mn->Draw("hist same");
         

         TLegend *leg4 = new TLegend(x_start, y_start, x_end, y_end);
            leg4->AddEntry(h1_mumu_pl,"f_{+} template","l");
            leg4->AddEntry(h1_mumu_mn,"f_{-} template","l");
            leg4->AddEntry(h1_mumu_alpha,"f_{#alpha} template","l");

            leg4->Draw();
            
            c_mumu2->Print(mu_fname2);
            delete c_mumu2;
}
}
