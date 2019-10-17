#ifndef TemplateUtils
#define TemplateUtils
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "../../utils/TemplateMaker_systematics.C"
#include "../../utils/root_files.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"

Double_t m_low;
Double_t m_high;

RooWorkspace *w;
RooRealVar *var = new RooRealVar("var", "var", 0.,n_cost_bins * n_xf_bins);
RooRealVar *var_ss = new RooRealVar("var_ss", "var_ss", 0.,(n_cost_bins/2) * n_xf_bins);
RooRealVar *Rqcd_ee_ss = new RooRealVar("Rqcd_ee_ss", "ee QCD normalization", 1, 0., 10.);
RooRealVar *Rqcd_mumu_ss = new RooRealVar("Rqcd_mumu_ss", "mumu QCD normalization", 1, 0., 10.);

TH1F *h_dummy = new TH1F("h_dummy", "", n_cost_bins*n_xf_bins, 0, n_cost_bins*n_xf_bins);



bool do_emu_scale = false;
bool do_RC = true;

void write_roo_hist(TH1F *h, RooRealVar *my_var){
    my_var->Print();
    h->Print("all");
    RooDataHist r(h->GetName(), h->GetName(), *my_var, h);
    r.Print();
    w->import(r);
    delete h;
}
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

TH1F* convert2d(TH2F *h_2d){
    int n_xf_bins = h_2d->GetNbinsX();
    int n_cost_bins = h_2d->GetNbinsY();

    TH1F *h_1d = new TH1F(h_2d->GetName(), "",  n_xf_bins * n_cost_bins, 0, n_xf_bins*n_cost_bins);
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins; j++){
            float content = h_2d->GetBinContent(i,j);
            float error = h_2d->GetBinError(i,j);
            int gbin = (i-1)*n_cost_bins + j;
            //printf("gbin %i: i j %i %i \n", gbin, i, j);
            h_1d->SetBinContent(gbin, content);
            h_1d->SetBinError(gbin, error);
        }
    }
    return h_1d;
}

int TwoDToOneDIdx(int nYbins, int i, int j){
    nYbins = nYbins/2;
    j = ((j-1) % nYbins) +1;
    int g_idx = (i-1)*nYbins + j;
    return g_idx;
}
int mirrorIdx(int nYbins, int j){
    j = ((j-1) % nYbins);
    int m_idx = nYbins - j;
    return m_idx;
}

void TwoDToSymIdxs(int nYbins, int i, int j, int &sym1_idx, int &sym2_idx){
    sym1_idx = (i-1)*nYbins + j;
    sym2_idx = (i-1)*nYbins + (nYbins - (j -1));
}


bool replace(std::string& str, const std::string& from, const std::string& to) {
//from https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void convert_qcd_to_param_hist(TH2F *h, FILE *f_log, float sign_scaling, int flag){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();
    RooArgList *bin_list_os = new RooArgList();
    RooArgList *bin_list_ss = new RooArgList();
    //h->Print("all");

    TH1F *h1 = convert2d(h);

    char h_name[40];
    char h_ss_name[40];
    char R_sign_param[40];
    sprintf(h_name, "%s_qcd_param", h->GetName());
    RooRealVar *R_qcd_sign_fraction;
    fprintf(f_log, "\n");
    sprintf(h_ss_name, "%s_ss_qcd_param", h->GetName());
    sprintf(R_sign_param, "R_%s_os_fakes", h->GetName());
    if(flag == FLAG_ELECTRONS){
        R_qcd_sign_fraction = new RooRealVar(R_sign_param, "Fraction of os fakes events", sign_scaling , 0., 1.);
        fprintf(f_log, "%s param %.4f 0.05 \n", R_sign_param, sign_scaling);
    }
    else{
        R_qcd_sign_fraction = new RooRealVar(R_sign_param, "Fraction of os fakes events", 1.0, 0.99, 1.01);
        fprintf(f_log, "%s param %.4f 0.0001 \n", R_sign_param, 1.0);
    }
    for(int i=1; i <= n_xf_bins; i++){
        for(int j=1; j <= n_cost_bins; j++){



            int g_idx = TwoDToOneDIdx(n_cost_bins, i, j);
            int sym1_idx, sym2_idx;
            TwoDToSymIdxs(n_cost_bins, i,j, sym1_idx, sym2_idx);
            //printf("i,j: %i %i ", i,j);
            //printf("g_idx, sym1, sym2: %i %i %i  \n", g_idx, sym1_idx, sym2_idx);

            double content = h1->GetBinContent(g_idx);
            double error = h1->GetBinError(g_idx);
            if(content<0) printf("Bin %i Content is %.0f \n", j, content);

            //printf("Bin %.1f error %.1f \n", content,error);
            char bin_name[40];
            char form_name_ss[40], form_name1_os[40], form_name2_os[40];
            sprintf(bin_name, "%s_bin%i",h_name, g_idx); 
            sprintf(form_name_ss, "%s_form_%i",h_ss_name, g_idx); 
            sprintf(form_name1_os, "%s_form_%i",h_name, sym1_idx); 
            sprintf(form_name2_os, "%s_form_%i",h_name, sym2_idx); 
            //prevent underflowing by fixing super small bins
            content = max(content, 0.001);
            if(flag == FLAG_MUONS){
                //muons only use opposite sign fakes
                content *= 2.;
                error *= 2.;
            }

            if (content < error){
                content = error/2.;
                error = 0.1*content;
            }
            else if(content < 2.5 * error){
                error = 0.3*content;
            }
            error = max(error, 0.0001);
            if(j<=(n_cost_bins/2)){
                //printf("first fill \n");
                RooRealVar *bin = new RooRealVar(bin_name, bin_name, content, 0., 10000.);
                fprintf(f_log, "%s param %.4f %.4f \n", bin_name, content, error);
                RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                RooFormulaVar *form_ss = new RooFormulaVar(form_name_ss, form_name_ss, "@0*(1.0 - @1)", RooArgList(*bin, *R_qcd_sign_fraction));
                //form1->Print();
                //form_ss->Print();
                bin_list->add(*bin);
                bin_list_ss->add(*form_ss);
                bin_list_os->add(*form1);
            }

            else{
                //printf("2nd fill \n");
                int old_j = sym2_idx % n_cost_bins;
                int old_g_idx = TwoDToOneDIdx(n_cost_bins, i, old_j);
                sprintf(bin_name, "%s_bin%i",h_name, old_g_idx); 
                //printf("Looking for bin %s \n", bin_name);
                RooRealVar *bin = (RooRealVar *) bin_list->find(bin_name);
                if(bin==nullptr) printf("NULL lookup of %s from bin list \n", bin_name);
                RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                //form1->Print();
                bin_list_os->add(*form1);
            }
        }
    
    }
    bin_list_ss->Print();
    bin_list_os->Print();
    char norm_ss_name[40], norm_name[40];
    sprintf(norm_ss_name, "%s_norm", h_ss_name);
    sprintf(norm_name, "%s_norm", h_name);
    RooAddition *norm_ss = new RooAddition(norm_ss_name, norm_ss_name, *bin_list_ss);

    RooAddition *norm = new RooAddition(norm_name, norm_name, *bin_list_os);

    RooParametricHist *p= new RooParametricHist (h_name, h_name, *var, *bin_list_os, *h_dummy);
    RooParametricHist *p_ss= new RooParametricHist (h_ss_name, h_ss_name, *var_ss, *bin_list_ss, *h1);

    
    w->import(*p_ss);
    w->import(*p, RooFit::RecycleConflictNodes());
    w->import(*norm,RooFit::RecycleConflictNodes());
    w->import(*norm_ss,RooFit::RecycleConflictNodes());
}



#endif
