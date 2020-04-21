#include "../../analyze/combine/make_templates.C"


void check_template(){
    
        int year = 2017;
        init(year);
        //setup_all_SFs(year);
        string sys_label = string("");

        char title[100];

        TH2F * h_mumu_plain = new TH2F("mumu_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_mumu_mc = new TH2F("mumu_dy", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);

        TH2F * h_elel_plain = new TH2F("elel_plain", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);
        TH2F * h_elel_mc = new TH2F("elel_dy", "", n_xf_bins, xf_bins, n_cost_bins, cost_bins);

        double m_low = 150.;
       
        double m_high = 171.;
        Double_t alpha_denom = amc_alpha[4];

        bool ss = false;



        TTree *mumu_ts[1] = {t_mumu_mc};
        printf("Making mumu \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_plain, year, m_low, m_high, FLAG_MUONS,   ss,  use_xf, "");
        one_mc_template(t_mumu_mc, alpha_denom, afb, h_mumu_mc, year, m_low, m_high, FLAG_MUONS,  use_xf, ""); 
        auto h1_mumu_back = convert2d(h_mumu_plain);
        auto h1_mumu_templ = convert2d(h_mumu_dy);

        TTree *elel_ts[1] = {t_elel_mc};
        printf("Making elel back \n");
        gen_combined_background_template(1, elel_ts, h_elel_plain, year, m_low, m_high, FLAG_ELECTRONS,  ss, use_xf, "");
        one_mc_template(t_elel_mc, alpha_denom, afb, h_elel_mc, year, m_low, m_high, FLAG_ELECTRONS, use_xf, ""); 
        auto h1_elel_back = convert2d(h_elel_plain);
        auto h1_elel_templ = convert2d(h_elel_dy);




        printf("MuMu Back is %.2f.Comb Temp is %.2f  \n", 
                h1_mumu_back->Integral(),  h1_mumu_templ->Integral());


        h1_mumu_back->SetLineColor(kGreen +3);
        h1_mumu_templ->SetLineColor(kBlack);
        h1_mumu_back->Print("all");
        h1_mumu_templ->Print("all");


        TCanvas *c_mumu = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
        h1_mumu_templ->Draw("hist");
        //h1_mumu_lo->Draw("hist same");
        h1_mumu_back->Draw("hist same");

        TLegend *leg = new TLegend(0.15, 0.3, 0.3, 0.15);
        leg->AddEntry(h1_mumu_templ, "Template", "f");
        //leg->AddEntry(h1_mumu_lo, "LO Template", "f");
        leg->AddEntry(h1_mumu_back, "Straight Hist", "f");
        leg->Draw();

        c_mumu->Print("mumu17_template_checka.png");


        
        printf("elel Back is %.2f.Comb Temp is %.2f   \n", 
                h1_elel_back->Integral(),  h1_elel_templ->Integral());

        h1_elel_back->SetLineColor(kGreen +3);
        h1_elel_templ->SetLineColor(kBlack);



        TCanvas *c_elel = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
        h1_elel_templ->Draw("hist");
        //h1_elel_lo->Draw("hist same");
        h1_elel_back->Draw("hist same");

        leg->Draw();

        c_elel->Print("elel17_template_checka.png");


        /*
        */

}




