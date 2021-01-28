#include "find_kl_limit.C"



TCanvas* test_Zp_can(FILE *f1, int M_Zp, Double_t cpl, Double_t *AFB_test){
    //hard code measured AFB's + uncertainties and SM predictions for AFB
    //return pvalue for a given Zp mass


    Double_t AFB_Zp[n_bins];

    read_AFBs(f1, AFB_Zp, M_Zp, cpl);
    /*
    printf("AFBs(%i, %.2f) are: ", M_Zp, cpl);
    for(int i = 0; i<n_bins; i++){
        printf("%.2f ", AFB_Zp[i]);
    }
    printf("\n");
    */

    //float scale = max(abs(test_stat(AFB_test, AFB_Zp)) + 10., abs(test_stat(AFB_Zp, AFB_Zp)) + 10);
    //printf("%.1f \n", scale);
    float scale = 40;

    int n_trials = 50000;

    int n_bins = 1000;
    float h_min = -scale;
    float h_max = scale;

    char h_title [100];

    sprintf(h_title, "Mass Z' = %i #Kappa_{L} = %.2f; test statistic ", M_Zp, cpl);

    TH1D * h_dist = new TH1D("h_dist1", h_title, n_bins,h_min,h_max);
    TRandom3 *r3 = new TRandom3();
    Double_t x[n_bins]; //randomly generated x vals
    //Find distribution of test statistic under Null hypothesis (ZPrime is there)
    for(int i=0; i<n_trials; i++){
        for(int j=0; j<n_bins; j++){
            x[j] = r3->Gaus(AFB_Zp[j], AFB_unc[j]);
        }
        Double_t test_val = test_stat(x, AFB_Zp);
        h_dist->Fill(test_val);
    }

    Double_t t_obs = test_stat(AFB_test, AFB_Zp);
    //prin2.2"T_obs is %.2f \n", t_obs);

    Double_t pval = get_pval(h_dist, t_obs);
    //if pval < 0.05 we will reject Null hypothesis (of Zprime existiing)
    
    //find test stat that is at p = 0.05
    float p_05 = 1.;
    float t_05 = h_max;
    for(t_05 = h_max - 1e-4; t_05 > h_min; t_05 -= (h_max - h_min)/n_bins){
        p_05 = get_pval(h_dist, t_05);
        if(p_05 > 0.05) break;
    }

    printf("t_obs  %.2f p-val %.3f \n", t_obs, pval);
    printf("t_05  %.2f p-05 %.3f \n", t_05, p_05);
   
    TCanvas *c1 = new TCanvas("c1", "", 800, 800);
    h_dist->Draw("hist");
    float draw_max = h_dist->GetMaximum();

    TLine *t_line = new TLine(t_obs, 0., t_obs, draw_max);
    t_line->SetLineColor(kRed);
    t_line->SetLineWidth(2);
    t_line->Draw("same");

    TLine *t05_line = new TLine(t_05, 0., t_05, draw_max);
    t05_line->SetLineColor(kBlue);
    t05_line->SetLineWidth(2);
    t05_line->SetLineStyle(kDashed);
    t05_line->Draw("same");

    auto latex = TLatex();
    latex.SetTextSize(0.025);
    int align;
    if(pval < 0.001) align = 33;
    else align = 13;
    latex.SetTextAlign(13);
    latex.SetNDC(true);
    char tt[100];
    sprintf(tt, "t_{exp} is %.1f, p-value %.2f", t_obs, pval);
    scale = (t_obs - h_min)/(h_max - h_min);
    latex.DrawLatex(scale, 0.75, tt);
    

    return c1;
}

void make_test_stat_plot(){

    FILE *f1 = fopen("AFBs.txt", "r");
    //test_Zp_can(f1, 2500., 1.0, AFB_SM);
    test_Zp_can(f1, 3000., 1.0, AFB_SM);
    //test_Zp_can(f1, 3500., 1.0, AFB_SM);
    //test_Zp_can(f1, 4000., 1.0, AFB_SM);
}

