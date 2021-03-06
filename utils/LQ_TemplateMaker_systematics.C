//construct templates from reconstructed events
//
//
//
#define _USE_MATH_DEFINES
#include <cmath>
#include "TempMaker.C"
//#include "ScaleFactors.C"

using namespace std;

//define constants

Double_t alpha = 1/127.9;
Double_t m_Z0 = 91.1875;
Double_t sin2_thetaw = 0.231; //sin^2(theta_W) (weak mixing angle)
Double_t G_F = 1.166e-5;
Double_t g_z = 2.4952; //width of Z0

//use coupling definitions from Quigg edition 1

//up quark
Double_t Q_q = 2./3. ;
Double_t I_3 = 1./2. ;

//down quark
//Q_q = -1./3.
//I_3 = -1./2.

Double_t crq = -2 *Q_q * sin2_thetaw;
Double_t clq = 2*I_3- 2. *Q_q * sin2_thetaw;

Double_t crl = 2 * sin2_thetaw;
Double_t cll = 2 * sin2_thetaw - 1;

Double_t cvq = crq + clq;
Double_t caq = -(clq - crq);

Double_t cvl = crl + cll;
Double_t cal = cll - crl;




//get this checked

void fixup_template_sum(TH3F *h_sym, TH3F *h_asym){
    //avoid negative pdfs by fixing up template edges
    int n_m_bins = h_sym->GetNbinsX();
    int n_xf_bins = h_sym->GetNbinsY();
    int n_cost_bins = h_sym->GetNbinsZ();
    for(int k=1; k<=n_m_bins; k++){
        for(int i=1; i<=n_xf_bins; i++){
            for(int j=1; j<= n_cost_bins/2; j++){
                //printf("%i %i \n", i,j);
                Double_t val_sym = h_sym->GetBinContent(k,i,j);
                Double_t val_asym = h_asym->GetBinContent(k,i,j);
                int opp_j = (n_cost_bins + 1) -j; //get this checked
                if(val_sym - 2*abs(val_asym) <= 0.){
                    //at LO needs to sym/2, give some cushion b/c of alpha term
                    Double_t val_asym_new = -val_sym/2.8;
                    h_asym->SetBinContent(k, i, j, val_asym_new);
                    h_asym->SetBinContent(k, i, opp_j, -val_asym_new);
                    printf("Fixed up bin %i %i %i. Old asym val was %.2f, new is %.2f \n\n", k, i, j, val_asym, val_asym_new);
                }
                if(val_sym < 1e-5){
                    val_sym = 1e-5;
                    h_sym->SetBinContent(k, i, j, val_sym);
                    h_sym->SetBinContent(k, i, opp_j, val_sym);
                }
            }

        }
    }
}

//get this checked
void cleanup_template(TH3F *h){
    //printf("%i %i \n", h->GetNbinsX(), h->GetNbinsY());
    for(int k=0; k<= h->GetNbinsX()+1; k++)
        for(int i=0; i<= h->GetNbinsY()+1; i++){
            for(int j=0; j<= h->GetNbinsZ()+1; j++){
                //printf("%i %i \n", i,j);
                Double_t min_val = 1e-6;
                Double_t val = h->GetBinContent(k,i,j);
                Double_t err = h->GetBinError(k,i,j);
                Double_t max_err = 0.7; //percent
                if(val< min_val){
                    h->SetBinContent(k,i,j,min_val);
                    h->SetBinError(k,i,j,err);
                }
                if(err > max_err * val){
                    //prevent val froming being close to fit boundary at 0
                    //By setting max stat error
                    err = val * max_err;
                    h->SetBinError(k,i,j, err);
                }



            }
        }
}
//make the f+ and f- templates
void make_pl_mn_templates(TH1* h_sym, TH1* h_asym, TH1* h_pl, TH1 *h_mn){

        h_pl->Add(h_sym, h_asym, 1., 1.);
        h_mn->Add(h_sym, h_asym, 1., -1.);
        h_pl->Scale(0.5);
        h_mn->Scale(0.5);
        for(int i=0; i<= h_sym->GetNbinsX() * h_sym->GetNbinsY(); i++){ //works for 1d an 2d histograms
            double err_rate = h_sym->GetBinError(i) / h_sym->GetBinContent(i);
            h_pl->SetBinError(i,err_rate * h_pl->GetBinContent(i));
            h_mn->SetBinError(i,err_rate * h_mn->GetBinContent(i));
        }
}
//changed
void print_hist(TH3 *h){
    printf("\n");

    for(int k=1; k<= h->GetNbinsX();k++){
        for(int i=1; i<= h->GetNbinsY(); i++){
            for(int j=1; j<= h->GetNbinsZ(); j++){
                printf("%.2e ",   (Double_t) h->GetBinContent(k,i,j));
            }
            printf("\n");
        }
    }
}
//
//static type means functions scope is only this file, to avoid conflicts
//changed
int gen_data_template(TTree *t1, TH3F* h,  
        int year, int flag1 = FLAG_MUONS, bool scramble_data = true, bool ss = false, bool use_xF = false){
    h->Sumw2();
    int nEvents = 0;
    int n=0;
    TempMaker tm(t1, true, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.setup();

    TRandom *rand;
    Double_t sign = 1.;
    //if(scramble_data) rand = new TRandom3();
    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);

        if( tm.met_pt < met_cut && tm.has_no_bjets && tm.not_cosmic){
            tm.doCorrections();
            Double_t var1 = abs(tm.cm.Rapidity());
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
                h->Fill(tm.m, var1, tm.cost, 1); //is tm.m ok here
            }
            else{
                h->Fill(tm.m, var1, -abs(tm.cost), 1);
            }
            nEvents++;
        }


    }

    tm.finish();
    return nEvents;
}

//input m_LQ in make_templates.C
int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F *h_asym, TH3F *h_alpha, TH3F *h_LQpure, TH3F *h_LQint,
        int year, Double_t m_LQ, int flag1 = FLAG_MUONS, bool use_xF = false, const string &sys_label = "" ){

    //printf("Setting up LQ rw helper... ");
    //LQ_rw_helper h_LQ;
    //setup_LQ_rw_helper(&h_LQ, year);
    
   // TH1F h_LQpure_weight = TH1F("h_LQpure_wt", "Weights distribution of LQpure", 200, 0., 450.)
    //TH1F h_LQpure_weight = TH1F("h_LQpure_wt", "Weights distribution of LQpure", 200, 0., 450.)
   //         n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    

    printf("Making mc template for sys %s \n", sys_label.c_str());

    h_sym->Sumw2(); //what is sumw2 -> create structure to store the sum of the squares of weights
    h_asym->Sumw2();
    h_alpha->Sumw2();
    h_LQpure->Sumw2();
    h_LQint->Sumw2();

    int n = 0;

    TempMaker tm(t1, false, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.is_gen_level = true;
   tm.do_ptrw = true;

    tm.setup_systematic(sys_label);
    tm.setup();

    Double_t max_obs = 0.;

   // Double_t reweight_LQpure_pos_min=10000.;
   // Double_t reweight_LQpure_pos_max=0.;
   // Double_t reweight_LQint_pos_min=10000.;
   // Double_t reweight_LQint_pos_max=0.;
   // int n_rogue=0;
   // int n_rogue_pos=0;
   // int n_rogue_neg=0;
    //printf("\n================denom/LQ_denom>=1000. for following entries================\n");
    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        //printf("%f", tm.getEvent(i));
        bool pass = (tm.m>=lq_m_bins[0]) && (tm.met_pt < met_cut)  && tm.has_no_bjets && tm.not_cosmic;
        if(pass){

            tm.doCorrections();
            tm.getEvtWeight();
            n++;
            Double_t gen_cost = tm.cost_st;
            Double_t var1 = abs(tm.cm.Rapidity());
            if(use_xF)  var1 = tm.xF;
            //Double_t denom = 3./8.*(1.+gen_cost*gen_cost + 0.5 * alpha_denom * (1. - 3. *gen_cost*gen_cost));
            Double_t denom = tm.getReweightingDenom();
            Double_t LQ_denom = tm.getLQReweightingDenom();
            if(LQ_denom ==0.) continue;
            /*
            if(denom/LQ_denom>=1000.)
            {
                printf("tm.m = %f, var1 = %f, tm.cost = %f \n",tm.m, var1, tm.cost);
            }
            */
            //if(LQ_denom==0.){printf("\n LQ denom is zero for m = %f",tm.m); LQ_denom = 1e-8;}
            //if(LQ_denom<0) printf("\n LQ_denom is negative : %f", LQ_denom);
            Double_t reweight_a = gen_cost/ denom;
            Double_t reweight_s = (1 + gen_cost*gen_cost)/denom;
            Double_t reweight_alpha = (1 - gen_cost*gen_cost)/denom;
            //for LQ, 2 terms-> pure and interference
            Double_t s = tm.m*tm.m;
            //(1./2.56819)*1e9 -> conversion of GeV^-2 to pb
            Double_t reweight_LQpure_norm = ((1./2.56819)*1e9/(128*M_PI*s));
            Double_t reweight_LQpure_num1 = ((1 - gen_cost)*(1 - gen_cost));
            Double_t reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1-gen_cost)* ((2*m_LQ*m_LQ/s)+1-gen_cost));
            Double_t reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            Double_t reweight_LQpure_pos;
            reweight_LQpure_pos = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);
          //  else reweight_LQpure_pos = (reweight_LQpure_norm*reweight_LQpure_num/denom);
          //  if(abs(reweight_LQpure_pos*tm.evt_weight)>30.)
            //if(LQ_denom==0.000001)
           //{
           // printf("for tm.m = %f, var1 = %f, tm.cost = %f \n",tm.m, var1, tm.cost);
           // printf("LQpure_wt = %f for tm.m = %f, var1 = %f, tm.cost = %f \n",(reweight_LQpure_pos*tm.evt_weight),tm.m, var1, tm.cost);
           // printf("reweight_LQpure_norm = %0.12f, reweight_LQpure_num = %0.12f, LQ_denom = %0.12f, reweight_LQpure_pos = %0.12f\n",reweight_LQpure_norm, reweight_LQpure_num, LQ_denom, reweight_LQpure_pos);
            //printf("==============\n");
            //n_rogue_pos++;
          //  reweight_LQpure_pos=0.;
           //}
           // if(reweight_LQpure_pos>50.)reweight_LQpure_pos = 0.;
            reweight_LQpure_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1+gen_cost)* ((2*m_LQ*m_LQ/s)+1+gen_cost));
            reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            Double_t reweight_LQpure_neg;
            reweight_LQpure_neg = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);
          //  else reweight_LQpure_neg = (reweight_LQpure_norm*reweight_LQpure_num/denom);
          // if(abs(reweight_LQpure_neg*tm.evt_weight)>30.)
            //if(LQ_denom==1e-6)
           //{
           // printf("LQpure_wt (-cost)= %f for tm.m = %f, var1 = %f, tm.cost = %f \n",(reweight_LQpure_neg*tm.evt_weight),tm.m, var1, tm.cost);
            //printf("reweight_LQpure_neg = %0.12f, tm.evt_weight = %.12f, LQpure_wt (-cost)= %.12f \n",reweight_LQpure_neg, tm.evt_weight, (reweight_LQpure_neg*tm.evt_weight));
           // printf("reweight_LQpure_norm = %0.12f, reweight_LQpure_num = %0.12f, LQ_denom = %0.12f, reweight_LQpure_neg = %0.12f\n",reweight_LQpure_norm, reweight_LQpure_num, LQ_denom, reweight_LQpure_neg);
           // printf("==============\n");
           // n_rogue_neg++;
          //  reweight_LQpure_neg=0.;
           //}
            //Double_t reweight_LQpure = (reweight_LQpure_num/LQ_denom);
            // Double_t reweight_LQpure = (reweight_LQpure_num/denom);
            // 
           // if(reweight_LQpure == 0.)printf("\n for m = %f reweight_LQpure = %f",tm.m, reweight_LQpure );
           // printf("\n LQ_denom = %0.12f",LQ_denom);
            Double_t reweight_LQint_norm1 = ((alpha*Q_q)/(16*s));
             //printf("\n LQint N1 = %0.12f",reweight_LQint_norm1);
            Double_t reweight_LQint_norm2_num = ((m_Z0*m_Z0-s)*(cal+cvl)*(caq-cvq)*G_F*m_Z0*m_Z0);
            //printf("\n LQint N2 num = %0.12f",reweight_LQint_norm2_num);
            Double_t reweight_LQint_norm2_denom = (128*1.4142*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            //printf("\n LQint N2 denom = %0.12f",reweight_LQint_norm2_denom);
            Double_t reweight_LQint_norm2 = (reweight_LQint_norm2_num/reweight_LQint_norm2_denom);
           // printf("\n LQint N2 = %0.12f",reweight_LQint_norm2);
            Double_t reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2)*(1./2.56819)*1e9;
             //printf("\n LQint N1+N2 = %0.12f",reweight_LQint_norm);
            Double_t reweight_LQint_num1 = ((1 - gen_cost)*(1 - gen_cost));
            Double_t reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1-gen_cost);
            Double_t reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            Double_t reweight_LQint_pos;
            reweight_LQint_pos = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);
            //else reweight_LQint_pos = (reweight_LQint_norm*reweight_LQint_num/denom);
          //  if(abs(reweight_LQint_pos*tm.evt_weight)>30.)
           // if(LQ_denom==1e-6)
          // {
            //printf("LQint_wt = %f for tm.m = %f, var1 = %f, tm.cost = %f \n",(reweight_LQint_pos*tm.evt_weight),tm.m, var1, tm.cost);
            //printf("reweight_LQint_pos = %0.12f, tm.evt_weight = %.12f, LQint_wt = %.12f \n",reweight_LQint_pos, tm.evt_weight, (reweight_LQint_pos*tm.evt_weight));
           // printf("reweight_LQint_norm = %0.12f, reweight_LQint_num = %0.12f, LQ_denom = %0.12f, reweight_LQint_pos = %0.12f\n",reweight_LQint_norm, reweight_LQint_num, LQ_denom, reweight_LQint_pos);
           // printf("==============\n");
           // n_rogue_pos++;
           // reweight_LQint_pos=0.;
          // }
           
          //  if(reweight_LQint_pos>50.)reweight_LQint_pos = 0.;
            reweight_LQint_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1+gen_cost);
            reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            Double_t reweight_LQint_neg;
            reweight_LQint_neg = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);
           // else reweight_LQint_neg = (reweight_LQint_norm*reweight_LQint_num/denom);
          // if(abs(reweight_LQint_neg*tm.evt_weight)>30.)
           // if(LQ_denom==1e-6)
           //{
           // printf("LQint_wt (-cost)= %f for tm.m = %f, var1 = %f, tm.cost = %f \n",(reweight_LQint_neg*tm.evt_weight),tm.m, var1, tm.cost);
           // printf("reweight_LQint_neg = %0.12f, tm.evt_weight = %.12f, LQint_wt (-cost)= %.12f \n",reweight_LQint_neg, tm.evt_weight, (reweight_LQint_neg*tm.evt_weight));
           // printf("reweight_LQint_norm = %0.12f, reweight_LQint_num = %0.12f, LQ_denom = %0.12f, reweight_LQint_neg = %0.12f\n",reweight_LQint_norm, reweight_LQint_num, LQ_denom, reweight_LQint_neg);
           // printf("==============\n");
           // n_rogue_neg++;
          //  reweight_LQint_neg=0.;
          // }
               
           // Double_t reweight_LQint = (reweight_LQint_num/LQ_denom);
            // Double_t reweight_LQint = (reweight_LQint_num/denom);
           //if(reweight_LQint == 0.) printf("\n for m = %f reweight_LQint = %f",tm.m, reweight_LQint);
           

       //   if((reweight_LQpure_pos*tm.evt_weight)<reweight_LQpure_pos_min)reweight_LQpure_pos_min = reweight_LQpure_pos*tm.evt_weight;
       //   if((reweight_LQpure_pos*tm.evt_weight)>reweight_LQpure_pos_max)reweight_LQpure_pos_max = reweight_LQpure_pos*tm.evt_weight;
       //   if((reweight_LQint_pos*tm.evt_weight)<reweight_LQint_pos_min)reweight_LQint_pos_min = reweight_LQint_pos*tm.evt_weight;
       //   if((reweight_LQint_pos*tm.evt_weight)>reweight_LQint_pos_max)reweight_LQint_pos_max = reweight_LQint_pos*tm.evt_weight;
       
           
           
            
          
            
            //modified these
            h_sym->Fill(tm.m, var1, tm.cost, reweight_s * tm.evt_weight); 
            h_sym->Fill(tm.m, var1, -tm.cost, reweight_s * tm.evt_weight); 

            h_asym->Fill(tm.m, var1, tm.cost, reweight_a * tm.evt_weight);
            h_asym->Fill(tm.m, var1, -tm.cost, -reweight_a * tm.evt_weight);

            h_alpha->Fill(tm.m, var1, tm.cost, reweight_alpha * tm.evt_weight); 
            h_alpha->Fill(tm.m, var1, -tm.cost, reweight_alpha * tm.evt_weight); 
            //LQ terms
            h_LQpure->Fill(tm.m, var1, tm.cost, reweight_LQpure_pos * tm.evt_weight); 
            h_LQpure->Fill(tm.m, var1, -tm.cost, reweight_LQpure_neg * tm.evt_weight);
           // h_LQpure_wt->Fill(reweight_LQpure_pos*tm.evt_weight);
           // h_LQpure_wt->Fill(reweight_LQpure_neg*tm.evt_weight);
         /*
            Int_t binx = h_LQpure->GetXaxis()->FindBin(tm.m);
            Int_t biny = h_LQpure->GetYaxis()->FindBin(var1);
            Int_t binz = h_LQpure->GetZaxis()->FindBin(tm.cost);
            if((h_LQpure->GetBinContent(binx,biny,binz))<0)
          {
            printf("LQpure_content negative for tm.m = %f, var1 = %f, tm.cost = %f \n",tm.m, var1, tm.cost);
            printf("reweight_LQpure = %.12f\n", reweight_LQpure);
            printf("tm.evt_weight = %.12f\n", tm.evt_weight);
            printf("LQ_denom = %.12f\n", LQ_denom);

          }
        */
            h_LQint->Fill(tm.m, var1, tm.cost, reweight_LQint_pos * tm.evt_weight); 
            h_LQint->Fill(tm.m, var1, -tm.cost, reweight_LQint_neg * tm.evt_weight);
          //  h_LQint_wt->Fill(reweight_LQint_pos*tm.evt_weight);
          //  h_LQint_wt->Fill(reweight_LQint_neg*tm.evt_weight);
          
        }
    }
/*
    printf("No. of +cost rogue events in %i = %i \n", year,n_rogue_pos);
    printf("No. of -cost rogue events in %i = %i \n", year,n_rogue_neg);
    printf("No. of +- cost rogue events = %i, n_rogue = %i\n",(n_rogue_pos+n_rogue_neg),n_rogue );
*/
    tm.finish();
    //int mbin = find_bin(m_bins, m_low + 0.1);
   // tm.fixRFNorm(h_sym, mbin); //not done for LQ
    //tm.fixRFNorm(h_asym, mbin);
    //tm.fixRFNorm(h_alpha, mbin);


    h_sym->Scale(0.5);
    h_asym->Scale(0.5);
    h_alpha->Scale(0.5);
    h_LQpure->Scale(0.5);
    h_LQint->Scale(0.5);

    //cleanup_template(h_sym);
    fixup_template_sum(h_sym, h_asym);
    t1->ResetBranchAddresses();
    printf("MC templates generated from %i events. Sym integral is %.1f \n \n", n, h_sym->Integral()); // what is this

   // printf("================LQpure template negative for following entries================\n");
    for(int k=1; k<=n_lq_m_bins; k++){    
        for(int i=1; i<=n_y_bins; i++){
            for(int j=1; j<= n_cost_bins; j++){
                Double_t pure_content = h_LQpure->GetBinContent(k,i,j);
               // Double_t int_content = h_LQint->GetBinContent(k,i,j);
               // int gbin = (k-1) * n_xf_bins*n_cost_bins + (i-1) * n_cost_bins + j;
               if(pure_content<0.){ 
               // printf("i_lqm_bin = %i, i_rap_bin = %i, i_cost_bin = %i, converted_index = %i, pure_content= %0.12f\n",k,i,j,gbin,pure_content );
                h_LQpure->SetBinContent(k,i,j,abs(pure_content));
                }
                Double_t int_content = h_LQint->GetBinContent(k,i,j);
               // int gbin = (k-1) * n_xf_bins*n_cost_bins + (i-1) * n_cost_bins + j;
               if(int_content<=0.){
              //  printf("i_lqm_bin = %i, i_rap_bin = %i, i_cost_bin = %i, converted_index = %i, int_content= %0.12f\n",k,i,j,gbin,int_content );
               h_LQint->SetBinContent(k,i,j,abs(int_content));
            }
        }
      }
    }

//  printf("year = %i, reweight_LQpure_pos_min = % 0.12f\n", year, reweight_LQpure_pos_min );
   // printf("year = %i, reweight_LQpure_pos_max = % 0.12f\n", year, reweight_LQpure_pos_max );
    //printf("year = %i, reweight_LQint_pos_min = % 0.12f\n", year, reweight_LQint_pos_min );
    //printf("year = %i, reweight_LQint_pos_max = % 0.12f\n", year, reweight_LQint_pos_max );
  
    return 0;
}
//changed, get checked?
int gen_combined_background_template(int nTrees, TTree **ts, TH3F* h,  
        int year,  int flag1 = FLAG_MUONS, 
        bool ss =false, bool use_xF = false,  const string &sys_label = ""){
    h->Sumw2();

    int nEvents = 0;

    for(int i=0; i<nTrees; i++){
        printf("Tree %i \n", i);
        TTree *t1 = ts[i];

        TempMaker tm(t1, false, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.is_gen_level = false;

        tm.setup_systematic(sys_label);
        tm.setup();


        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass =  tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

            if(pass){
                tm.doCorrections();
                tm.getEvtWeight();

                Double_t var1 = abs(tm.cm.Rapidity());
                if(use_xF)  var1 = tm.xF;
                nEvents++;
                if(!ss) h->Fill(tm.m, var1, tm.cost, tm.evt_weight);
                else{
                    h->Fill(tm.m, var1, -abs(tm.cost), tm.evt_weight);
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
//get checked
int one_mc_template(TTree *t1, Double_t afb, TH3F* h_dy, 
        int year, Double_t m_LQ, int flag1 = FLAG_MUONS, bool use_xF = false, bool use_LQ_denom=true,
        const string &sys_label = "" ){

    int n_var1_bins = n_y_bins;
    Float_t *var1_bins = y_bins;
    if(use_xF){
        n_var1_bins = n_xf_bins;
        var1_bins = xf_bins;
    }

    TH3F h_sym = TH3F("h_sym", "Symmetric template of mc",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_sym.SetDirectory(0);
    TH3F h_alpha = TH3F("h_alpha", "Gauge boson polarization template of mc",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_alpha.SetDirectory(0);
    TH3F h_asym = TH3F("h_asym", "Asymmetric template of mc",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_asym.SetDirectory(0);
    TH3F h_LQpure = TH3F("h_LQpure", "LQpure template of mc",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_LQpure.SetDirectory(0);
    TH3F h_LQint = TH3F("h_LQint", "LQint template of mc",
            n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
    h_LQint.SetDirectory(0);
    //includes m_LQ
  //  gen_mc_template(t1, &h_sym, &h_asym, &h_alpha, &h_LQpure, &h_LQint, year, m_LQ, flag1,  use_xF, use_LQ_denom,sys_label);


   // double norm = 3./4./(2.+alpha);

    auto h_pl = h_sym + h_asym;
    auto h_mn = h_sym - h_asym;
    h_pl.Scale(0.5);
    h_mn.Scale(0.5);

    //h_alpha.Scale(norm * alpha);
    //h_LQpure.Scale(norm * alpha);
    //h_LQint.Scale(norm * alpha);
    

    h_dy->Add(&h_pl, &h_mn);
    h_dy->Add(&h_alpha);
    h_dy->Add(&h_LQpure);
    h_dy->Add(&h_LQint);

    return 1;
}




/*
//what to change in this
std::pair<Double_t, Double_t> gen_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree* t_QCD_contam, TH3F *h, 
        int year,  Double_t m_low, Double_t m_high, int flag1 = FLAG_MUONS, bool incl_ss = true, bool ss_binning = false, bool use_xF = false){
    h->Sumw2();
    TH3D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    double tot_weight_os = 0.;
    double tot_weight_ss = 0.;
    if(flag1 == FLAG_MUONS) setup_new_mu_fakerate(&FR, year);
    else setup_new_el_fakerate(&FR, year);
    //TH3D *FR;
    h_err = (TH3D *) FR.h->Clone("h_err");
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

                Double_t var1 = abs(tm.cm.Rapidity());
                if(use_xF)  var1 = tm.xF;

                if(!ss_binning) h->Fill(tm.m, var1, tm.cost, tot_evt_weight);
                else{
                    h->Fill(tm.m, var1, -abs(tm.cost), tot_evt_weight);
                }
                if(opp_sign) tot_weight_os += tot_evt_weight;
                else tot_weight_ss += tot_evt_weight;


            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());

    }
    //remove negative bins
    cleanup_template(h);

    set_fakerate_errors(h_err, FR.h, h);
    printf("Total Fakerate weight Weight is %.2f \n", h->Integral());
    // remove outliers
    tot_weight_os = std::max(1e-4, tot_weight_os);
    tot_weight_ss = std::max(1e-4, tot_weight_ss);
    Double_t scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
    if(!incl_ss) scaling = 1.;
    
    Double_t err;
    Double_t integ = h->IntegralAndError(1, h->GetNbinsX(), 1, h->GetNbinsY(), err);
    printf("Total fakerate est is %.0f +/- %.0f \n", integ, err);
    return std::make_pair(scaling, err/integ);
}

void gen_emu_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TH3F *h, 
        int year, Double_t m_low = 150., Double_t m_high = 10000., bool ss = false){
    //need to change this to TempMaker class at some point...
    FakeRate el_FR, mu_FR;
    //TH3D *FR;
    setup_new_el_fakerate(&el_FR, year);
    setup_new_mu_fakerate(&mu_FR, year);
    //FR.h->Print();
    for (int l=0; l<=2; l++){
        TTree *t;
        if (l==0) t = t_WJets;
        if (l==1) t = t_QCD;
        if (l==2) t = t_WJets_MC;
        Double_t m, xF, cost, jet1_btag, jet2_btag, gen_weight;
        Double_t jet1_pt, jet2_pt, pu_SF;

        Double_t el_id_SF, el_reco_SF;
        Double_t era1_iso_SF, era1_id_SF;
        Double_t era2_iso_SF, era2_id_SF;
        Double_t evt_fakerate, lep1_fakerate, lep2_fakerate, el1_eta, el1_pt, mu1_eta, mu1_pt;
        TLorentzVector *el = 0;
        TLorentzVector *mu = 0;
        Int_t iso_lep, no_bjets;
        Double_t_t met_pt, el_charge, mu_charge;
        Int_t nJets;
        nJets = 2;
        pu_SF=1;
        t->SetBranchAddress("el1_charge", &el_charge);
        t->SetBranchAddress("mu1_charge", &mu_charge);
        t->SetBranchAddress("met_pt", &met_pt);
        t->SetBranchAddress("jet2_btag", &jet2_btag);
        t->SetBranchAddress("jet1_btag", &jet1_btag);
        t->SetBranchAddress("jet1_pt", &jet1_pt);
        t->SetBranchAddress("jet2_pt", &jet2_pt);
        //t1->SetBranchAddress("evt_fakerate", &evt_fakerate);
        //t1->SetBranchAddress("el_fakerate", &el1_fakerate);
        t->SetBranchAddress("has_nobjets", &no_bjets);
        t->SetBranchAddress("el1_pt", &el1_pt);
        t->SetBranchAddress("mu1_pt", &mu1_pt);
        t->SetBranchAddress("el1_eta", &el1_eta);
        t->SetBranchAddress("mu1_eta", &mu1_eta);
        t->SetBranchAddress("nJets", &nJets);
        t->SetBranchAddress("el", &el);
        t->SetBranchAddress("mu", &mu);

        if(l==0 || l==2 ){
            t->SetBranchAddress("iso_lep", &iso_lep);
        }
        if(l==2){
            t->SetBranchAddress("el_id_SF", &el_id_SF);
            t->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t->SetBranchAddress("gen_weight", &gen_weight);
            t->SetBranchAddress("pu_SF", &pu_SF);
            t->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t->SetBranchAddress("era2_id_SF", &era2_id_SF);
        }

        Long64_t size  =  t->GetEntries();
        for (int i=0; i<size; i++) {
            t->GetEntry(i);
            if(l==0){
                //iso_lep: 0 = muons, 1 electrons
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                if(iso_lep ==0) lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = lep1_fakerate/(1-lep1_fakerate);
            }
            if(l==1){
                lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                lep2_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                evt_fakerate = -(lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));
            }
            if(l==2){
                Double_t era1_SF =  era1_id_SF;
                Double_t era2_SF =  era2_id_SF;
                Double_t mu_SF, mu_lumi;
                if(year ==2016){
                    mu_SF = (era1_SF *bcdef_lumi16 + era2_SF * gh_lumi16)/(bcdef_lumi16 + gh_lumi16);
                    mu_lumi=mu_lumi16;
                }
                if(year==2017){
                    mu_SF = era1_SF;
                    mu_lumi=mu_lumi17;
                }
                if(year==2018){
                    mu_SF = era1_SF;
                    mu_lumi=mu_lumi18;
                }

                Double_t mc_weight = gen_weight * el_id_SF * el_reco_SF *pu_SF * mu_SF  * 1000. * mu_lumi;
                if(iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(mu1_pt, mu1_eta, mu_FR.h);
                if(iso_lep ==0) lep1_fakerate = get_new_fakerate_prob(el1_pt, el1_eta, el_FR.h);
                evt_fakerate = -(lep1_fakerate * mc_weight)/(1-lep1_fakerate);
            }

            TLorentzVector cm = *el + *mu;
            bool opp_sign = ((abs(el_charge - mu_charge)) > 0.01);
            Double_t m = cm.M();
            bool pass = m>= m_low && m <= m_high && met_pt < met_cut  && no_bjets && opp_sign;
            if(pass){
                Double_t xF = compute_xF(cm); 
                cost = get_cost(*el, *mu);
                if(ss) h->Fill(xF, -abs(cost), evt_fakerate);
                else{
                    h->Fill(xF, cost, evt_fakerate);
                }
            }
        }

        printf("After iter %i current fakerate est is %.0f \n", l, h->Integral());
    }
    cleanup_template(h);
    printf("Total fakerate est is %.0f \n", h->Integral());
}


void gen_emu_template(TTree *t1, TH3F *h, 
        bool is_data = false, int year=2016, Double_t m_low = 150., Double_t m_high = 999999., bool ss= true){
    Long64_t size  =  t1->GetEntries();
    Double_t m, xF, cost, mu1_pt, mu2_pt, jet1_btag, jet2_btag, gen_weight;
    Double_t era1_HLT_SF, era1_iso_SF, era1_id_SF;
    Double_t era2_HLT_SF, era2_iso_SF, era2_id_SF, el_id_SF, el_reco_SF, el_HLT_SF;
    Double_t jet1_pt, jet2_pt, jet1_b_weight, jet2_b_weight, pu_SF;
    jet1_b_weight = jet2_b_weight =1.;
    Double_t_t met_pt, el_charge, mu_charge;
    Int_t nJets, no_bjets;
    nJets = 2;
    pu_SF=1;
    TLorentzVector *el=0;
    TLorentzVector *mu=0;
    int nEvents =0;

    t1->SetBranchAddress("el", &el);
    t1->SetBranchAddress("mu", &mu);
    t1->SetBranchAddress("el1_charge", &el_charge);
    t1->SetBranchAddress("mu1_charge", &mu_charge);
    t1->SetBranchAddress("met_pt", &met_pt);
    t1->SetBranchAddress("jet2_btag", &jet2_btag);
    t1->SetBranchAddress("jet1_btag", &jet1_btag);
    t1->SetBranchAddress("jet1_pt", &jet1_pt);
    t1->SetBranchAddress("jet2_pt", &jet2_pt);
    t1->SetBranchAddress("nJets", &nJets);
    t1->SetBranchAddress("has_nobjets", &no_bjets);
    if(!is_data){
        t1->SetBranchAddress("gen_weight", &gen_weight);
        t1->SetBranchAddress("pu_SF", &pu_SF);
        t1->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
        t1->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
        t1->SetBranchAddress("era1_id_SF", &era1_id_SF);
        t1->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
        t1->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
        t1->SetBranchAddress("era2_id_SF", &era2_id_SF);
        t1->SetBranchAddress("el_id_SF", &el_id_SF);
        t1->SetBranchAddress("el_reco_SF", &el_reco_SF);     
    }

    for (int i=0; i<size; i++) {
        t1->GetEntry(i);
        bool opp_sign = ((abs(el_charge - mu_charge)) > 0.01);
        TLorentzVector cm;
        cm = *el + *mu;
        m = cm.M();
        double cost = get_cost(*el, *mu);
        xF  = compute_xF(cm); 

        bool pass = m>= m_low && m <= m_high && met_pt < met_cut  && no_bjets && opp_sign;

        if(pass){
            if(is_data){
                if(ss) h->Fill(xF, -abs(cost), 1.);
                else h->Fill(xF, cost, 1.);
            }
            else{
                nEvents++;
                Double_t evt_weight = gen_weight * pu_SF * el_id_SF * el_reco_SF;
                Double_t era1_weight = era1_iso_SF * era1_id_SF ;
                Double_t era2_weight = era2_iso_SF * era2_id_SF ;
                era1_weight *= era1_HLT_SF;
                era2_weight *= era2_HLT_SF;

                if (nJets >= 1){
                    evt_weight *= jet1_b_weight;
                }
                if (nJets >= 2){
                    evt_weight *= jet2_b_weight;
                }
                //printf(" %.2e %.2e %.2e \n", evt_weight, evt_weight *era1_weight, evt_weight *era2_weight);
                Double_t tot_weight;
                if(year ==2016) tot_weight = 1000 * (evt_weight * era2_weight * gh_lumi16 + evt_weight * era1_weight * bcdef_lumi16);
                if(year ==2017) tot_weight = 1000 * evt_weight * era1_weight * mu_lumi17;
                if(year ==2018) tot_weight = 1000 * evt_weight * era1_weight * mu_lumi18;
                if(ss) h->Fill(xF, -abs(cost), tot_weight);
                else h->Fill(xF, cost, tot_weight);

            }



        }
    }
    cleanup_template(h);
} */
