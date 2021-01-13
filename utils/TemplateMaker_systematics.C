//construct templates from reconstructed events
//
//
//
#include "TempMaker.C"


using namespace std;


void symmetrize2d(TH2F *h_2d){
    int n_xf_bins = h_2d->GetNbinsX();
    int n_cost_bins = h_2d->GetNbinsY();

    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins/2; j++){
            float content = h_2d->GetBinContent(i,j);
            float error = h_2d->GetBinError(i,j);

            int opp_j = (n_cost_bins + 1) -j;
            float content2 = h_2d->GetBinContent(i,opp_j);
            float error2 = h_2d->GetBinError(i,opp_j);

            float new_content = (content + content2)/2.0;
            float new_error = pow((error*error + error2*error2), 0.5)/2.0;
            h_2d->SetBinContent(i,j, new_content);
            h_2d->SetBinContent(i,opp_j, new_content);

            h_2d->SetBinError(i,j, new_error);
            h_2d->SetBinError(i,opp_j, new_error);


        }
    }
}


void fixup_template_sum(TH2F *h_sym, TH2F *h_asym){
    //avoid negative pdfs by fixing up template edges
    int n_xf_bins = h_sym->GetNbinsX();
    int n_cost_bins = h_sym->GetNbinsY();

    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins/2; j++){
            //printf("%i %i \n", i,j);
            float val_sym = h_sym->GetBinContent(i,j);
            float val_asym = h_asym->GetBinContent(i,j);
            int opp_j = (n_cost_bins + 1) -j;
            if(val_sym - 2*abs(val_asym) <= 0.){
                //at LO needs to sym/2, give some cushion b/c of alpha term
                float val_asym_new = -val_sym/2.8;
                h_asym->SetBinContent(i, j, val_asym_new);
                h_asym->SetBinContent(i, opp_j, -val_asym_new);
                printf("Fixed up bin %i %i. Old asym val was %.2f, new is %.2f \n\n", i, j, val_asym, val_asym_new);
            }
            if(val_sym < 1e-5){
                val_sym = 1e-5;
                h_sym->SetBinContent(i, j, val_sym);
                h_sym->SetBinContent(i, opp_j, val_sym);
            }
        }

    }
}



void cleanup_template(TH2F *h){
    //printf("%i %i \n", h->GetNbinsX(), h->GetNbinsY());
    for(int i=0; i<= h->GetNbinsX()+1; i++){
        for(int j=0; j<= h->GetNbinsY()+1; j++){
            //printf("%i %i \n", i,j);
            float min_val = 1e-6;
            float val = h->GetBinContent(i,j);
            float err = h->GetBinError(i,j);
            float max_err = 0.7; //percent
            if(val< min_val){
                val = min_val;
                h->SetBinContent(i,j,min_val);
                h->SetBinError(i,j,err);
            }
            if(err > max_err * val){
                //prevent val froming being close to fit boundary at 0
                //By setting max stat error
                err = val * max_err;
                h->SetBinError(i,j, err);
            }



        }
    }
}


void cleanup_fakes_template(TH2F *h){
    //printf("%i %i \n", h->GetNbinsX(), h->GetNbinsY());
    TH2F *h_copy = (TH2F *) h->Clone("clone");
    for(int i=1; i<= h->GetNbinsX(); i++){
        for(int j=1; j<= h->GetNbinsY(); j++){
            float small_val = 0.5;
            float min_val = 1e-6;
            float val = h->GetBinContent(i,j);
            float err = h->GetBinError(i,j);
            float max_err = 0.7; //percent
            if(val< min_val){
                if(err > 0.5 && abs(val)/err < 1.){
                    //bin content is negative/small but uncertainty overlaps with 0
                    // change bin content to be positive so that this
                    // uncertainty can be in the template (can't put an unc. on zero)
                    val = err / 2.;
                    h->SetBinContent(i,j, val);
                }

                else if(val < min_val){
                    //zero this bin
                    val = min_val;
                    h->SetBinContent(i,j,min_val);
                }
            }

            if(err > max_err * val){
                //prevent val froming being close to fit boundary at 0
                //By setting max error as a fraction of the bin content
                err = val * max_err;
                h->SetBinError(i,j, err);
            }

        }
    }
    h_copy->Delete();
}


//set error on dest histogram to match fractional error on source histogram
void set_frac_error(TH1 *h_source, TH1 *h_dest){
    for(int i=0; i<= h_source->GetNbinsX() * h_source->GetNbinsY(); i++){ //works for 1d and 2d histograms
        double err_rate = h_source->GetBinError(i) / h_source->GetBinContent(i);
        h_dest->SetBinError(i,err_rate * h_dest->GetBinContent(i));
    }
}


void make_pl_mn_templates(TH1* h_sym, TH1* h_asym, TH1* h_pl, TH1 *h_mn){

        h_pl->Add(h_sym, h_asym, 1., 1.);
        h_mn->Add(h_sym, h_asym, 1., -1.);
        h_pl->Scale(0.5);
        h_mn->Scale(0.5);
        set_frac_error(h_sym, h_pl);
        set_frac_error(h_sym, h_mn);
}


void print_hist(TH2 *h){
    printf("\n");
    for(int i=1; i<= h->GetNbinsX(); i++){
        for(int j=1; j<= h->GetNbinsY(); j++){
            printf("%.2e ",   (float) h->GetBinContent(i,j));
        }
        printf("\n");
    }

}
//
//static type means functions scope is only this file, to avoid conflicts

int gen_data_template(TTree *t1, TH2F* h,  
        int year, Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool scramble_data = true, bool ss = false, bool use_xF = false){
    h->Sumw2();
    int nEvents = 0;
    int n=0;
    TempMaker tm(t1, true, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.setup();

    TRandom *rand;
    float sign = 1.;
    if(scramble_data) rand = new TRandom3();
    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);

        if(tm.m >= m_low && tm.m <= m_high && tm.met_pt < met_cut && tm.has_no_bjets && tm.not_cosmic){
            tm.doCorrections();
            float var1 = abs(tm.cm.Rapidity());
            if(use_xF)  var1 = tm.xF;

            n++;
            if(!ss){
                if(scramble_data){
                    //switch + and - back and forth
                    tm.cost = sign * std::fabs(tm.cost); 
                    sign *= -1.;
                    //randomly flip data events
                    //if(rand->Uniform(1.) > 0.5) tm.cost = std::fabs(tm.cost);
                    //else tm.cost = -std::fabs(tm.cost);
                }
                h->Fill(var1, tm.cost, 1); 
            }
            else{
                h->Fill(var1, -abs(tm.cost), 1);
            }
            nEvents++;
        }


    }

    tm.finish();
    return nEvents;
}



int gen_mc_template(TTree *t1, TH2F* h_sym, TH2F *h_asym, TH2F *h_alpha,
        int year, Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool use_xF = false,
        const string &sys_label = "" ){
    printf("Making mc template for sys %s \n", sys_label.c_str());

    h_sym->Sumw2();
    h_asym->Sumw2();
    h_alpha->Sumw2();

    int n = 0;

    TempMaker tm(t1, false, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.is_gen_level = true;
    tm.do_ptrw = true;

    tm.setup_systematic(sys_label);
    tm.setup();

    float max_obs = 0.;

    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;
        if(pass){

            tm.doCorrections();
            tm.getEvtWeight();
            n++;
            float gen_cost = tm.cost_st;
            //float denom = 3./8.*(1.+gen_cost*gen_cost + 0.5 * alpha_denom * (1. - 3. *gen_cost*gen_cost));
            float denom = tm.getReweightingDenom();
            float reweight_a = gen_cost/ denom;
            float reweight_s = (1 + gen_cost*gen_cost)/denom;
            float reweight_alpha = (1 - gen_cost*gen_cost)/denom;


            float var1 = abs(tm.cm.Rapidity());
            if(use_xF)  var1 = tm.xF;

            

            h_sym->Fill(var1, tm.cost, reweight_s * tm.evt_weight); 
            h_sym->Fill(var1, -tm.cost, reweight_s * tm.evt_weight); 

            h_asym->Fill(var1, tm.cost, reweight_a * tm.evt_weight);
            h_asym->Fill(var1, -tm.cost, -reweight_a * tm.evt_weight);

            h_alpha->Fill(var1, tm.cost, reweight_alpha * tm.evt_weight); 
            h_alpha->Fill(var1, -tm.cost, reweight_alpha * tm.evt_weight); 

        }
    }
    //printf("Max obs is %.3f \n", max_obs);

    tm.finish();
    int mbin = find_bin(m_bins, m_low + 0.1);
    tm.fixRFNorm(h_sym, mbin, year);
    tm.fixRFNorm(h_asym, mbin, year);
    tm.fixRFNorm(h_alpha, mbin, year);

    h_sym->Scale(0.5);
    h_asym->Scale(0.5);
    h_alpha->Scale(0.5);
    //cleanup_template(h_sym);
    fixup_template_sum(h_sym, h_asym);
    t1->ResetBranchAddresses();
    printf("MC templates generated from %i events. Sym integral is %.1f \n \n", n, h_sym->Integral());
    return 0;
}

int gen_combined_background_template(int nTrees, TTree **ts, TH2F* h,  
        int year, Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, 
        bool ss =false, bool use_xF = false, bool emu_reweight = false,  const string &sys_label = ""){
    h->Sumw2();

    int nEvents = 0;

    for(int i=0; i<nTrees; i++){
        printf("Tree %i \n", i);
        TTree *t1 = ts[i];

        TempMaker tm(t1, false, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.is_gen_level = false;
        tm.do_emu_costrw = emu_reweight;

        tm.setup_systematic(sys_label);
        tm.setup();

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

            if(pass){
                tm.doCorrections();
                tm.getEvtWeight();

                float var1 = abs(tm.cm.Rapidity());
                if(use_xF)  var1 = tm.xF;
                nEvents++;
                if(!ss) h->Fill(var1, tm.cost, tm.evt_weight);
                else{
                    h->Fill(var1, -abs(tm.cost), tm.evt_weight);
                }
            }
        }
        tm.finish();
    }


    printf("Performing templ. cleanup (removing neg. bins) \n");
    cleanup_template(h);
    printf("Tot Weight is %.2f from %i events \n", h->Integral(), nEvents);
    return 0;
}

int one_mc_template(TTree *t1, Double_t a0, Double_t afb, TH2F* h_dy, 
        int year, Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool use_xF = false,
        const string &sys_label = "" ){

    int n_var1_bins = n_y_bins;
    float *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }

    TH2F h_sym = TH2F("h_sym", "Symmetric template of mc",
            n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_sym.SetDirectory(0);
    TH2F h_alpha = TH2F("h_alpha", "Gauge boson polarization template of mc",
            n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_alpha.SetDirectory(0);
    TH2F h_asym = TH2F("h_asym", "Asymmetric template of mc",
            n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_asym.SetDirectory(0);

    gen_mc_template(t1, &h_sym, &h_asym, &h_alpha, year, m_low, m_high, flag1,  use_xF, sys_label);


    float alpha = 2.* a0/ (2. - a0);
    double norm = 3./4./(2.+alpha);

    TH2F *h_pl = (TH2F *) h_sym.Clone("h_pl");
    TH2F *h_mn = (TH2F *) h_sym.Clone("h_pl");
    h_pl->Reset();
    h_mn->Reset();
    make_pl_mn_templates(&h_sym, &h_asym, h_pl, h_mn);


    

    h_dy->Add(h_pl, h_mn, (norm + afb), (norm - afb));
    h_alpha.Scale(norm * alpha);
    h_dy->Add(&h_alpha);

    return 1;
}






void gen_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree* t_QCD_contam, TH2F *h, 
        int year,  Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool incl_ss = true, bool ss_binning = false, bool use_xF = false, string sys_label = ""){
    h->Sumw2();

    int shape_sys = 0;

    if(!sys_label.empty()){
        printf("Making fakes template for sys %s \n", sys_label.c_str());
        int sys_shift = 0;

        if(sys_label.find("Up") != string::npos){
            sys_shift = 1;
        }
        else if(sys_label.find("Down") != string::npos){
            sys_shift = -1;
        }
        if(sys_label.find("fakesrw") != string::npos){
            int foo;
            
            if(sys_shift > 0){
                if(flag1 == FLAG_MUONS) sscanf(sys_label.c_str(), "_mufakesrw%ib%iUp", &shape_sys, &foo);
                else sscanf(sys_label.c_str(), "_elfakesrw%ib%iUp", &shape_sys, &foo);

            }
            else{
                if(flag1 == FLAG_MUONS) sscanf(sys_label.c_str(), "_mufakesrw%ib%iDown", &shape_sys, &foo);
                else sscanf(sys_label.c_str(), "_elfakesrw%ib%iDown", &shape_sys, &foo);
                shape_sys *= -1;
            }
            printf("Doing fakes costrw sys %i \n", shape_sys);
        }
    }
    


    TH2D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    double tot_weight_os = 0.;
    double tot_weight_ss = 0.;
    if(flag1 == FLAG_MUONS) setup_new_mu_fakerate(&FR, year);
    else setup_new_el_fakerate(&FR, year);
    //TH2D *FR;
    h_err = (TH2D *) FR.h->Clone("h_err");
    h_err->Reset();
    for (int l=0; l<=3; l++){
        TTree *t;
        bool is_data, is_one_iso;
        if (l==0){
            t = t_WJets;
            is_data =true;
            is_one_iso = true;
        }
        if (l==1){
            t = t_QCD;
            is_data =true;
            is_one_iso = false;
        }
        if (l==2){
            t = t_WJets_contam;
            is_data =false;
            is_one_iso = true;
        }
        if (l==3){
            t = t_QCD_contam;
            is_data =false;
            is_one_iso = false;
        }
        TempMaker tm(t, is_data, year);
        tm.is_one_iso = is_one_iso;
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.setup();

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);

            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

            bool opp_sign = false;
            if(flag1 == FLAG_MUONS) opp_sign = ((abs(tm.mu1_charge - tm.mu2_charge)) > 0.01);
            else opp_sign = ((abs(tm.el1_charge - tm.el2_charge)) > 0.01);

            if(!incl_ss) pass = pass && opp_sign; 
            if(pass){
                double evt_reweight = 0.;

                if(flag1 == FLAG_MUONS){
                    Double_t mu1_fakerate, mu2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = mu1_fakerate/(1-mu1_fakerate);
                    }

                    if(l==1){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                    }
                    if(l==2){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate )/(1-mu1_fakerate);
                    }
                    if(l==3){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                    }
                    if((l==0) && tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), tm.getEvtWeight());
                    if((l==0) && tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), tm.getEvtWeight());
                }
                else{
                    Double_t el1_fakerate, el2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = el1_fakerate/(1-el1_fakerate);
                    }
                    if(l==1){
                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                    }
                    if(l==2){

                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate )/(1-el1_fakerate);
                    }
                    if(l==3){

                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));
                    }
                    if((l==0)  && tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), tm.getEvtWeight());
                    if((l==0)  && tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), tm.getEvtWeight());
                }
                double tot_evt_weight = 0.;
                if(is_data) tot_evt_weight = evt_reweight; 
                else{
                    tot_evt_weight = evt_reweight * tm.getEvtWeight();
                }

                float var1 = abs(tm.cm.Rapidity());
                if(use_xF)  var1 = tm.xF;

                if(!ss_binning) h->Fill(var1, tm.cost, tot_evt_weight);
                else{
                    h->Fill(var1, -abs(tm.cost), tot_evt_weight);
                }
                if(opp_sign) tot_weight_os += tot_evt_weight;
                else tot_weight_ss += tot_evt_weight;


            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());

    }

    //scale to just OS normalization
    tot_weight_os = std::max(1e-4, tot_weight_os);
    tot_weight_ss = std::max(1e-4, tot_weight_ss);
    float scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
    if(!incl_ss) scaling = 1.;

    h->Scale(scaling);
    //symmetrize to avg statistical fluctuations before removing negative bins
    symmetrize2d(h);
    //remove negative bins
    printf("Before:\n");
    h->Print("range");
    cleanup_fakes_template(h);
    printf("After: \n");
    h->Print("range");


    fakes_costrw_helper h_rw;
    setup_fakes_costrw_helper(&h_rw, year);
    if(flag1 == FLAG_MUONS) fakes_cost_reweight(h, h_rw.mu_rw, shape_sys);
    else fakes_cost_reweight(h, h_rw.el_rw, shape_sys);
    
    set_fakerate_errors(h_err, FR.h, h);
    Double_t err;
    Double_t integ = h->IntegralAndError(1, h->GetNbinsX(), 1, h->GetNbinsY(), err);
    printf("Total fakerate est is %.0f +/- %.0f \n", integ, err);
    return;
}

