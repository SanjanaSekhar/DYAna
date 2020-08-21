#ifndef LQ_TemplateUtils
#define LQ_TemplateUtils
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "../../utils/LQ_TemplateMaker_systematics.C"
#include "../../utils/root_files.h"

#ifndef STAND_ALONE
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#endif
//root files are accessed in root_files.h
//this file uses no of xf and cost bins so they are accessed from elsewhere
Double_t m_low;
Double_t m_high;
bool use_xF = false;

RooWorkspace *w;
RooRealVar *var = new RooRealVar("var", "var", 0.,n_cost_bins * n_xf_bins * n_m_bins);
RooRealVar *var_ss = new RooRealVar("var_ss", "var_ss", 0.,(n_cost_bins/2) * n_xf_bins *n_m_bins);
RooRealVar *Rqcd_ee_ss = new RooRealVar("Rqcd_ee_ss", "ee QCD normalization", 1, 0., 10.);
RooRealVar *Rqcd_mumu_ss = new RooRealVar("Rqcd_mumu_ss", "mumu QCD normalization", 1, 0., 10.);

TH1F *h_dummy = new TH1F("h_dummy", "", n_cost_bins*n_xf_bins*n_m_bins, 0, n_cost_bins*n_xf_bins*n_m_bins);



bool do_emu_scale = false;
bool do_RC = true;



void write_roo_hist(TH1F *h, RooRealVar *my_var){
    my_var->Print();
    //h->Print("all");
    RooDataHist r(h->GetName(), h->GetName(), *my_var, h);
    //r.Print();
    w->import(r);
    delete h;
} 
//changed
void symmetrize3d(TH3F *h_3d){ //this function is called in make_templates.C on 8 2D hists that we create
    int n_m_bins = h_3d->GetNbinsX();
    int n_xf_bins = h_3d->GetNbinsY();
    int n_cost_bins = h_3d->GetNbinsZ();

    for(int k=1; k<=n_m_bins; k++){	
	    for(int i=1; i<=n_xf_bins; i++){
	        for(int j=1; j<= n_cost_bins/2; j++){
	            float content = h_3d->GetBinContent(k,i,j);
	            float error = h_3d->GetBinError(k,i,j);

	            int opp_j = (n_cost_bins + 1) -j;
	            float content2 = h_3d->GetBinContent(k,i,opp_j);
	            float error2 = h_3d->GetBinError(k,i,opp_j);

	            float new_content = (content + content2)/2.0;
	            float new_error = pow((error*error + error2*error2), 0.5)/2.0;
	            h_3d->SetBinContent(k,i,j, new_content);
	            h_3d->SetBinContent(k,i,opp_j, new_content);

	            h_3d->SetBinError(k,i,j, new_error);
	            h_3d->SetBinError(k,i,opp_j, new_error);
	        }
	    }
	}    
}

//changed but doubt
TH1F* convert3d(TH3F *h_3d){
    int n_m_bins = h_3d->GetNbinsX();
    int n_xf_bins = h_3d->GetNbinsY();
    int n_cost_bins = h_3d->GetNbinsZ();

    TH1F *h_1d = new TH1F(h_3d->GetName(), "",  n_xf_bins * n_cost_bins * n_m_bins, 0, n_xf_bins*n_cost_bins*n_m_bins);// 0 is the 1st numbering of the bin
	for(int k=1; k<=n_m_bins; k++){    
	    for(int i=1; i<=n_xf_bins; i++){
	        for(int j=1; j<= n_cost_bins; j++){
	            float content = h_3d->GetBinContent(k,i,j);
	            float error = h_3d->GetBinError(k,i,j);
	            int gbin = (k-1) * n_xf_bins*n_cost_bins + (i-1) * n_cost_bins + j; //trying to convert 3 indices to 1
	            //printf("gbin %i: i j %i %i \n", gbin, i, j);
	            h_1d->SetBinContent(gbin, content);
	            h_1d->SetBinError(gbin, error);
	        }
	    }
	}
    return h_1d;
}
/*
int one_idx(int i, int j, int n_binsx, int n_binsy){
   //lose 2 bins for each row above mid-row
   
   if(i <= n_binsx/2) return (k-1) * n_xf_bins*n_cost_bins + (i-1) * n_binsy + j;
   if(j == n_binsy) j-=1;
   if(j>1) j-=1;

   int base = (n_binsx/2) * n_binsy;
   return base + std::max((i - n_binsx/2 -1), 0)* (n_binsy-2) + j;

}

TH1F* convert3d(TH3F *h_3d){
    int n_m_bins = h_3d->GetNbinsX();
    float n_binsx = h_3d->GetNbinsY();
    float n_binsy = h_3d->GetNbinsZ();
  //  int n_1d_bins = get_n_1d_bins(n_binsx, n_binsy);

    int n_1d_bins = std::round(std::ceil(n_binsx/2.) * n_binsy + std::floor(n_binsx/2.) * (n_binsy-2));

    TH1F *h_1d = new TH1F(h_3d->GetName(), "",  n_1d_bins*n_m_bins, 0, n_1d_bins*n_m_bins);// 0 is the 1st numbering of the bin
    for(int k=1; k<=n_m_bins; k++){    
        for(int i=1; i<=n_binsx; i++){
            for(int j=1; j<= n_binsy; j++){
            float content = h_3d->GetBinContent(k,i,j);
            float error = h_3d->GetBinError(k,i,j);
            int gbin = one_idx(i,j, n_binsx, n_binsy);
            
            //add in case previous bin filled
            float content_1d = h_1d->GetBinContent(gbin); 
            float error_1d = h_1d->GetBinError(gbin); 
            h_1d->SetBinContent(gbin, content_1d + content);
            h_1d->SetBinError(gbin, std::pow(error_1d*error_1d + error*error, 0.5));
            }
        }
    }
    return h_1d;
}
*/
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
/*
#ifndef STAND_ALONE
void convert_qcd_to_param_hist(TH3F *h, FILE *f_log, float sign_scaling, float sign_scale_err, int flag){
    //convert a hist to a parametric hist 
    RooArgList *bin_list = new RooArgList();
    RooArgList *bin_list_os = new RooArgList();
    RooArgList *bin_list_ss = new RooArgList();
    //h->Print("all");

    TH1F *h1 = convert3d(h);

    char h_name[40];
    char h_ss_name[40];
    char R_sign_param[40];
    sprintf(h_name, "%s_qcd_param", h->GetName());
    RooRealVar *R_qcd_sign_fraction;
    fprintf(f_log, "\n");
    sprintf(h_ss_name, "%s_ss_qcd_param", h->GetName());
    sprintf(R_sign_param, "R_%s_os_fakes", h->GetName());
    if(sign_scaling > 0.99){
        sign_scaling = 0.995;
        sign_scale_err = 0.0001;
    }
    if(flag == FLAG_ELECTRONS){
        R_qcd_sign_fraction = new RooRealVar(R_sign_param, "Fraction of os fakes events", sign_scaling , 0., 0.9995);
        fprintf(f_log, "%s param %.4f %.5f \n", R_sign_param, sign_scaling, sign_scale_err);
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
                bin_list->add(*bin);
                fprintf(f_log, "%s param %.4f %.4f \n", bin_name, content, error);
                //for os region multiply by 0.5 because there are two bins (both signs)
                if(flag == FLAG_ELECTRONS){
                    RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                    RooFormulaVar *form_ss = new RooFormulaVar(form_name_ss, form_name_ss, "@0*(1.0 - @1)", RooArgList(*bin, *R_qcd_sign_fraction));
                    bin_list_os->add(*form1);
                    bin_list_ss->add(*form_ss);
                }
                else{
                    RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0", RooArgList(*bin));
                    bin_list_os->add(*form1);
                }

                //form1->Print();
                //form_ss->Print();
            }

            else{
                //printf("2nd fill \n");
                int old_j = sym2_idx % n_cost_bins;
                int old_g_idx = TwoDToOneDIdx(n_cost_bins, i, old_j);
                sprintf(bin_name, "%s_bin%i",h_name, old_g_idx); 
                //printf("Looking for bin %s \n", bin_name);
                RooRealVar *bin = (RooRealVar *) bin_list->find(bin_name);
                if(bin==nullptr) printf("NULL lookup of %s from bin list \n", bin_name);

                if(flag == FLAG_ELECTRONS){
                    RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0*@1", RooArgList(*bin, *R_qcd_sign_fraction));
                    bin_list_os->add(*form1);
                }
                else{
                    RooFormulaVar *form1 = new RooFormulaVar(form_name1_os, form_name1_os, "0.5*@0", RooArgList(*bin));
                    bin_list_os->add(*form1);
                }

                //form1->Print();
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
#endif*/



#endif