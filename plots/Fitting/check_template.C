#define STAND_ALONE
#include "../../analyze/combine/make_templates.C"
#include "../../utils/PlotUtils.C"


void check_template(){
    
        int year = 2017;
        double m_low = 150.;
        double m_high = 171.;
        bool use_xf = false;

        init(year);
        setup_all_SFs(year);
        string sys_label = string("");

        char title[100];

        int n_var1_bins = n_y_bins;
        float *var1_bins = y_bins;
        if(use_xF){
            n_var1_bins = n_xf_bins;
            var1_bins = xf_bins;
        }

        TH2F * h_mumu_plain = new TH2F("mumu_plain", "", n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH2F * h_mumu_mc = new TH2F("mumu_dy", "", n_var1_bins, var1_bins, n_cost_bins, cost_bins);

        TH2F * h_elel_plain = new TH2F("elel_plain", "", n_var1_bins, var1_bins, n_cost_bins, cost_bins);
        TH2F * h_elel_mc = new TH2F("elel_dy", "", n_var1_bins, var1_bins, n_cost_bins, cost_bins);


        bool ss = false;

        float afb = 0.616;
        float A0 = 0.062;



        TTree *mumu_ts[1] = {t_mumu_mc};
        printf("Making mumu \n");
        gen_combined_background_template(1, mumu_ts, h_mumu_plain, year, m_low, m_high, FLAG_MUONS,   ss,  use_xf, "");
        one_mc_template(t_mumu_mc, A0, afb, h_mumu_mc, year, m_low, m_high, FLAG_MUONS,  use_xf, ""); 
        //one_mc_template(t_mumu_mc, A0, afb, h_mumu_mc, year, m_low, m_high, FLAG_MUONS,  use_xf, ""); 
        auto h1_mumu_back = convert2d(h_mumu_plain);
        auto h1_mumu_templ = convert2d(h_mumu_mc);

        TTree *elel_ts[1] = {t_elel_mc};
        printf("Making elel back \n");
        gen_combined_background_template(1, elel_ts, h_elel_plain, year, m_low, m_high, FLAG_ELECTRONS,  ss, use_xf, "");
        one_mc_template(t_elel_mc, A0, afb, h_elel_mc, year, m_low, m_high, FLAG_ELECTRONS, use_xf, ""); 
        auto h1_elel_back = convert2d(h_elel_plain);
        auto h1_elel_templ = convert2d(h_elel_mc);



        printf("MuMu Back is %.2f.Comb Temp is %.2f  \n", 
                h1_mumu_back->Integral(),  h1_mumu_templ->Integral());
        
        TCanvas *c_mumu = make_ratio_plot("Template Test", h1_mumu_back, "Raw", h1_mumu_templ, "Sum of Templates", "ratio", "cos(#theta)", false, false, 0.9, 1.1);


        //c_mumu->Print("mumu17_template_checka.png");


        
        printf("elel Back is %.2f.Comb Temp is %.2f   \n", 
                h1_elel_back->Integral(),  h1_elel_templ->Integral());

        TCanvas *c_elel = make_ratio_plot("Template Test", h1_elel_back, "Raw", h1_elel_templ, "Sum of Templates", "ratio", "cos(#theta)", false, false, 0.9, 1.1);

        //c_elel->Print("elel17_template_checka.png");


        /*
        */

}




