#include "../plots/tdrstyle.C"
#include "../plots/CMS_lumi.C"
#include "find_kl_limit.C"

Double_t get_kl_limit(FILE *f1, int M_Zp, Double_t kl_start, Double_t *AFB_test, bool print = false){
    //printf("kl_start %.3f \n", kl_start);
    Double_t kl_min = 0.05;
    Double_t kl_max = 2.05;
    Double_t kl_step = 0.05;
    Double_t alpha = 0.05;
    Double_t pval = 0.;
    Double_t kl;
    Double_t eps = 0.0001;


    for(kl = kl_start; kl >= kl_min-eps && kl <= kl_max + eps; kl-=kl_step){
        //printf("%.2f kl \n", kl);
        pval = test_Zp(f1, M_Zp, kl, AFB_test);
        if(print) printf("kl, kl_start, pval = %.2f %.2f %.3f \n", kl, kl_start, pval);
        if (pval > alpha && kl < kl_start){
            if(print) printf("return \n");
            return kl + kl_step;
        }
        else if(pval  > alpha){
            if(print) printf("Started too low, changing kl start from %.2f to %.2f\n", kl_start, kl_start +2*kl_step);
            kl_start = kl_start + 2*kl_step;
            kl = kl_start + kl_step;
            //printf("%i %.2f \n", m, kl);
        }
    }

    //cant find limit, return minimum
    Double_t ret_kl;
    if( kl >= kl_max) ret_kl = kl_max;
    if( kl <= kl_min) ret_kl = kl_min;
    printf("Couldn't find limit for M=%i, last tried %.2f. Returning %.2f  \n", M_Zp, kl, ret_kl);
    return ret_kl;
}
Double_t get_stddev(int n_vals, Double_t *vals){
    Double_t mean, var;
    int n_entries = n_vals;
    mean = 0.;
    var = 0.;
    for(int i =0; i< n_vals; i++){
        //printf("%.2f \n", vals[i]);
        if(vals[i] < 0. || std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    //printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vals; i++){
        if(vals[i] < 0. || std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    //printf("std %.3f \n\n\n", sqrt(var));
    return sqrt(var);
}
TGraph *makeAFillGraph(int len, Double_t *x, Double_t *y1, Double_t *y2, int linecolor=1, int fillcolor=0, int fillstyle=0){
    vector<Double_t> g_x;
    vector<Double_t> g_y;
    for(int i =0; i<len; i++){
        printf("filling in graph for m=%.2f \n", x[i]);
        g_x.push_back(x[i]);
        g_y.push_back(y1[i]);
    }
    for(int i =len -1; i>=0; i--){
        g_x.push_back(x[i]);
        g_y.push_back(y2[i]);
    }
    TGraph * gr = new TGraph(2*len, g_x.data(), g_y.data());
    gr->SetLineColor(linecolor);
    gr->SetFillColor(fillcolor);
    gr->SetFillStyle(fillstyle);
    return gr;
}

void make_limit_plot(){
    int n_m_bins = 100;
    int n_trials = 400;
    int m_start = 1500;
    int m_max = 4500;
    //int m_max = 2000;
    int m_step = 100;
    int m;
    int i=0;
    bool inc_meas_limit = true;
    Double_t kl_start = 1.;
    Double_t kl_min = 0.05;
    Double_t kl_step = 0.05;
    FILE *f1 = fopen("AFBs.txt", "r");
    TRandom3 *r3 = new TRandom3();

    Double_t masses[n_m_bins], kl_limit[n_m_bins], kl_expected_lim[n_m_bins], kl_expected_lim_up[n_m_bins], 
             kl_expected_lim_upup[n_m_bins], kl_expected_lim_down[n_m_bins],    
             kl_expected_lim_downdown[n_m_bins], kl_lim_trials[n_m_bins][n_trials],
             kl_limit_stds[n_m_bins], kl_limit_2stds[n_m_bins];
    Double_t meas_lim =1.0;

    for(m = m_start; m <=m_max; m+=m_step){
        printf("trying m=%i \n", m);
        if(inc_meas_limit){
            meas_lim = get_kl_limit(f1, m, kl_start, AFB_measured, true);
            if(meas_lim <= 0.04 || meas_lim >=2.){
                printf("bad lim for m=%i meas = %.2f \n", m, meas_lim); 
                continue;
            }
            if(i>0) meas_lim = max(meas_lim, kl_limit[i-1]);
            kl_limit[i] = meas_lim;
        }
        //Double_t exp_lim  = get_kl_limit(f1, m, kl_start, AFB_SM);
        Double_t exp_lim  = get_kl_limit(f1, m, kl_start, AFB_exp);
        if(i>0) exp_lim = max(exp_lim, kl_expected_lim[i-1]);
        printf("EXP LIM IS %.2f \n", exp_lim);
        masses[i] = m;
        kl_expected_lim[i] = exp_lim;
        //generate random expected limits
        Double_t AFB_rand[n_bins];
        for(int j=0; j<n_trials; j++){
            for(int k=0; k<n_bins; k++){
                AFB_rand[k] = r3->Gaus(AFB_SM[k], AFB_unc[k]);
            }
            kl_lim_trials[i][j] = get_kl_limit(f1, m, kl_expected_lim[i]+2*kl_step, AFB_rand);
        }
        Double_t std_dev = get_stddev(n_trials, kl_lim_trials[i]);
        kl_limit_stds[i] = std::min(exp_lim, std_dev);
        kl_limit_2stds[i] = std::min(exp_lim, 2. * std_dev);
        

        //printf("std dev, std dev %.3f %.3f \n", get_stddev(n_trials, kl_lim_trials[i]), get_stddev(n_trials, kl_lim_trials[i]));
        //printf("std dev %.3f, 1sig %.3f  2sig %.3f  \n", std_dev, kl_limit_stds[i], kl_limit_2stds[i]);
        i++;
    }
    setTDRStyle();
    TGraph *exp = new TGraph(i, masses, kl_expected_lim);    
    TGraphErrors *one_sig = new TGraphErrors(i, masses, kl_expected_lim, nullptr, kl_limit_stds);
    TGraphErrors *two_sig = new TGraphErrors(i, masses, kl_expected_lim, nullptr, kl_limit_2stds);
    //TGraph *one_sig = makeAFillGraph(i, masses, kl_expected_lim_down, kl_expected_lim_up, 0, kGreen+1, 1001);
    //TGraph *two_sig = makeAFillGraph(i, masses, kl_expected_lim_downdown, kl_expected_lim_upup, 0, kOrange, 1001);
    //
    one_sig->SetFillColor(kGreen+1);
    two_sig->SetFillColor(kOrange);

    //exp->Print();
    //one_sig->Print();
    //two_sig->Print();


    exp->SetLineColor(1);
    exp->SetLineStyle(2);
    exp->SetLineWidth(4);
    auto legbb = new TLegend(0.45,0.55,0.90,0.85);
	legbb->SetFillStyle(1001);
	legbb->SetFillColor(0);    
	legbb->SetBorderSize(0);  
	legbb->SetTextSize(0.040);
	legbb->SetTextFont(42);
	legbb->AddEntry(exp,"Expected","l");
	legbb->AddEntry(one_sig,"#pm 1 std. deviation","f");
	legbb->AddEntry(two_sig,"#pm 2 std. deviation","f");



    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    two_sig->SetTitle("Z' Limits");
    two_sig->GetYaxis()->SetTitle("#kappa_{L}(g_{Z'}/g_{Z})");
    two_sig->GetYaxis()->CenterTitle();
    two_sig->GetYaxis()->SetTitleOffset(1.0);
    two_sig->GetXaxis()->SetTitle("Z' mass (GeV)");

    two_sig->Draw("AL3SAME");
    one_sig->Draw("3SAME");
    exp->Draw("LSAME");

    TGraph *meas = new TGraph(i, masses, kl_limit);    
    if(inc_meas_limit){
        meas->SetLineColor(1);
        meas->SetLineStyle(1);
        meas->SetLineWidth(4);
        legbb->AddEntry(meas,"Observed","l");
        meas->Draw("LSAME");
    }

    printf("idx, mass, expected limit, expected std dev obs limit \n");
    for(int idx=0; idx <i; idx++){
        printf("%i %.2f %.2f %.2f %.2f \n", idx, masses[idx], kl_expected_lim[idx], kl_limit_stds[idx], kl_limit[idx]);
    }


    legbb->Draw();
    
    writeExtraText = false;
    extraText = "Preliminary";
    int iPeriod = -1; 
    CMS_lumi( c1, iPeriod, 0 );


    //c1->Print("limit_obs_test.png");
    //c1->Print("limit_obs_test.pdf");

    TFile *f_out = TFile::Open("limit_obs.root", "RECREATE");
    f_out->cd();
    exp->SetName("exp_limit");
    two_sig->SetName("two_sigma_band");
    one_sig->SetName("one_sigma_band");
    meas->SetName("obs_limit");
    two_sig->Write();
    one_sig->Write();
    exp->Write();
    meas->Write();
    f_out->Write();
    f_out->Close();

}
