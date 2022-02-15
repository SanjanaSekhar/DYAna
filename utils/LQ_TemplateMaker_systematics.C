//construct templates from reconstructed events
//
//
//
#define _USE_MATH_DEFINES
#include <cmath>
#include "TempMaker.C"
//#include "ScaleFactors.C"
//#include "HistUtils.C"
using namespace std;

//define constants

float alpha = 1/127.9;
float m_Z0 = 91.1875;
float sin2_thetaw = 0.231; //sin^2(theta_W) (weak mixing angle)
float G_F = 1.166e-5;
float g_z = 2.4952; //width of Z0

//use coupling definitions from Quigg edition 1
float crl = 2 * sin2_thetaw;
float cll = 2 * sin2_thetaw - 1;

float cvl = crl + cll;
float cal = cll - crl;

//up quark
float Q_u = 2./3. ;
float I3_u = 1./2. ;

float crq_u = -2 *Q_u * sin2_thetaw;
float clq_u = 2*I3_u- 2. *Q_u * sin2_thetaw;

float cvq_u = crq_u + clq_u;
float caq_u = -(clq_u - crq_u);

//down quark
float Q_d = -1./3. ;
float I3_d = -1./2. ;

float crq_d = -2 *Q_d * sin2_thetaw;
float clq_d = 2*I3_d- 2. *Q_d * sin2_thetaw;

float cvq_d = crq_d + clq_d;
float caq_d = -(clq_d - crq_d);

float Q_q, caq, cvq;

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

void fixup_template_sum(TH3F *h_sym, TH3F *h_asym){
    //avoid negative pdfs by fixing up template edges
  int n_m_bins = h_sym->GetNbinsX();
  int n_xf_bins = h_sym->GetNbinsY();
  int n_cost_bins = h_sym->GetNbinsZ();
  for(int k=1; k<=n_m_bins; k++){
    for(int i=1; i<=n_xf_bins; i++){
      for(int j=1; j<= n_cost_bins/2; j++){
                //printf("%i %i \n", i,j);
        float val_sym = h_sym->GetBinContent(k,i,j);
        float val_asym = h_asym->GetBinContent(k,i,j);
                int opp_j = (n_cost_bins + 1) -j; //get this checked
                if(val_sym - 2*abs(val_asym) <= 0.){
                    //at LO needs to sym/2, give some cushion b/c of alpha term
                  float val_asym_new = -val_sym/2.8;
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
                float min_val = 1e-6;
                float val = h->GetBinContent(k,i,j);
                float err = h->GetBinError(k,i,j);
                float max_err = 0.7; //percent
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
          void cleanup_fakes_template(TH3F *h){
    //printf("%i %i \n", h->GetNbinsX(), h->GetNbinsY());
            TH3F *h_copy = (TH3F *) h->Clone("clone");
            for(int k=0; k<= h->GetNbinsX(); k++){
              for(int i=1; i<= h->GetNbinsY(); i++){
                for(int j=1; j<= h->GetNbinsZ(); j++){
                  float small_val = 0.5;
                  float min_val = 1e-6;
                  float val = h->GetBinContent(k,i,j);
                  float err = h->GetBinError(k,i,j);
            float max_err = 0.7; //percent
            if(val< min_val){
              if(err > 0.5 && abs(val)/err < 1.){
                    //bin content is negative/small but uncertainty overlaps with 0
                    // change bin content to be positive so that this
                    // uncertainty can be in the template (can't put an unc. on zero)
                val = err / 2.;
                h->SetBinContent(k,i,j, val);
              }

              else if(val < min_val){
                    //zero this bin
                val = min_val;
                h->SetBinContent(k,i,j,min_val);
              }
            }

            if(err > max_err * val){
                //prevent val froming being close to fit boundary at 0
                //By setting max error as a fraction of the bin content
              err = val * max_err;
              h->SetBinError(k,i,j, err);
            }
          }
          }
        }
        h_copy->Delete();
      }
//set error on dest histogram to match fractional error on source histogram
      void set_frac_error(TH1 *h_source, TH1 *h_dest){
    for(int i=0; i<= h_source->GetNbinsX() * h_source->GetNbinsY() * h_source->GetNbinsZ(); i++){ //works for 1d and 2d histograms
      double err_rate = h_source->GetBinError(i) / h_source->GetBinContent(i);
      h_dest->SetBinError(i,err_rate * h_dest->GetBinContent(i));
    }
  }
  void make_pl_mn_templates(TH1* h_sym, TH1* h_asym, TH1* h_pl, TH1 *h_mn){

    h_pl->Add(h_sym, h_asym, 1., 1.);
    h_mn->Add(h_sym, h_asym, 1., -1.);
    h_pl->Scale(0.5);
    h_mn->Scale(0.5);
    //set_frac_error(h_sym, h_pl);
    //set_frac_error(h_sym, h_mn);
  }
//changed
  void print_hist(TH3 *h){
    printf("\n");

    for(int k=1; k<= h->GetNbinsX();k++){
      for(int i=1; i<= h->GetNbinsY(); i++){
        for(int j=1; j<= h->GetNbinsZ(); j++){
          printf("%.2e ",   (float) h->GetBinContent(k,i,j));
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
    float sign = 1.;
    if(scramble_data) rand = new TRandom3();
    for (int i=0; i<tm.nEntries; i++) {
      tm.getEvent(i);

      if( tm.met_pt < met_cut && tm.has_no_bjets && tm.not_cosmic){
        tm.doCorrections();
        float var1 = abs(tm.cm.Rapidity());
        if(use_xF)  var1 = tm.xF;

        n++;
        if(!ss){
          if(scramble_data){
                    //switch + and - back and forth
                    //tm.cost = sign * std::fabs(tm.cost); 
                    //sign *= -1.;
                    //randomly flip data events
            if(rand->Uniform(1.) > 0.5) tm.cost = std::fabs(tm.cost);
            else tm.cost = -std::fabs(tm.cost);
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

float get_LQ_denom(float gen_cost,float s,float Q_q, float caq, float cvq){
   float XS1 = (M_PI*pow(alpha,2)*pow(Q_q,2)*(pow(gen_cost,2)+1))/(2*s);
            //pure Z0 term
          float XS2_num = ((((cal*caq*pow(gen_cost,2)+ cal*caq+ 8*gen_cost*cvl*cvq)*caq +(pow(gen_cost,2)+1)*cal*pow(cvq,2))*cal+(pow(caq,2)+pow(cvq,2))*(pow(gen_cost,2)+1)*pow(cvl,2))*pow(G_F,2)*pow(m_Z0,4)*s);
          float XS2_denom = (256*M_PI*(pow((m_Z0*m_Z0-s),2) + pow(g_z*m_Z0,2)));
          float XS2 = XS2_num/ XS2_denom;
            //Z0 gamma interference
          float XS45_num =  - ((gen_cost*gen_cost+1)*cvl*cvq + 2*cal*caq*gen_cost) * (m_Z0*m_Z0-s) * alpha*G_F*m_Z0*m_Z0*Q_q;
          float XS45_denom = (8*sqrt(2)*(pow((m_Z0*m_Z0-s),2)+pow((g_z*m_Z0),2)));
          float XS45 = XS45_num/XS45_denom;
            float LQ_denom = (XS1 + XS2 + XS45); //new LQdenom is basically just LO SM
            return LQ_denom;
}

//input m_LQ in make_templates.C
        int gen_mc_template(TTree *t1, TH3F* h_sym, TH3F *h_asym, TH3F *h_alpha, TH3F *h_LQpure_u, TH3F *h_LQint_u,TH3F *h_LQpure_d, TH3F *h_LQint_d,
          int year, float m_LQ, int flag1 = FLAG_MUONS, bool use_xF = false, const string &sys_label = "" ){

          printf("Making mc template for sys %s \n", sys_label.c_str());

    h_sym->Sumw2(); //what is sumw2 -> create structure to store the sum of the squares of weights
    h_asym->Sumw2();
    h_alpha->Sumw2();
    h_LQpure_u->Sumw2();
    h_LQint_u->Sumw2();
    h_LQpure_d->Sumw2();
    h_LQint_d->Sumw2();

    int n = 0;

    TempMaker tm(t1, false, year);
    if(flag1 == FLAG_MUONS) tm.do_muons = true;
    else tm.do_electrons = true;
    tm.is_gen_level = true;
    tm.do_ptrw = true;

    tm.btag_mc_eff_idx = 1; //idx for DY MC btag effs

    tm.setup_systematic(sys_label);
    tm.setup();

    float max_obs = 0.;

    for (int i=0; i<tm.nEntries; i++) {
      tm.getEvent(i);
      bool pass = (tm.m >= lq_m_bins[0]) && (tm.met_pt < met_cut)  && tm.has_no_bjets && tm.not_cosmic;
      if(pass){

        tm.doCorrections();
        tm.getEvtWeight();
        n++;
        float gen_cost = tm.cost_st;
        float var1 = abs(tm.cm.Rapidity());
        if(use_xF)  var1 = tm.xF;

        

            // SM terms
          
            //(1./2.56819)*1e9 -> conversion of GeV^-2 to pb
            // jacobian for conversion = 2sqrt(x1x2/s); x1x2s=gen_m^2; sqrt(s) = CM energy -> 13 TeV
          float s = tm.gen_m*tm.gen_m;
          float n_conv = (1./2.56819)*1e9;
          float LQ_jacobian = (2*tm.gen_m*1e-6)/(13.*13.);

           // float LQ_denom = tm.getLQReweightingDenom(flag_q);
            //if(LQ_denom==0.) {
              //printf("\nhello flag_q = %i, tm.m = %f, rap = %f, cost = %f\n",flag_q,tm.m,var1,tm.cost); 
            //continue;}


          float denom = tm.getReweightingDenom();
          float reweight_a = gen_cost/denom;
          float reweight_s = (1 + gen_cost*gen_cost)/denom;
          float reweight_alpha = (1 - gen_cost*gen_cost)/denom;

            //fill SM temps

          h_sym->Fill(tm.m, var1, tm.cost, reweight_s* tm.evt_weight ); 
          h_sym->Fill(tm.m, var1, -tm.cost, reweight_s * tm.evt_weight ); 

          h_asym->Fill(tm.m, var1, tm.cost, reweight_a * tm.evt_weight );
          h_asym->Fill(tm.m, var1, -tm.cost, -reweight_a * tm.evt_weight );

          h_alpha->Fill(tm.m, var1, tm.cost, reweight_alpha * tm.evt_weight ); 
          h_alpha->Fill(tm.m, var1, -tm.cost, reweight_alpha * tm.evt_weight ); 

            /* 
            
            float reweight_s_norm1 = (M_PI*alpha*alpha*Q_q*Q_q)/(2*s);
            float reweight_s_norm2 = ((cal*cal + cvl*cvl)*(caq*caq + cvq*cvq)*G_F*G_F*pow(m_Z0,4)*s)/(256*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_s_norm3 = (cvl*cvq*(m_Z0*m_Z0-s)*alpha*G_F*m_Z0*m_Z0*Q_q)/(8*sqrt(2)*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_s_norm = (reweight_s_norm1 + reweight_s_norm2 - reweight_s_norm3)*n_conv*LQ_jacobian;
            reweight_s = reweight_s_norm*(1 + gen_cost*gen_cost)/LQ_denom;

            
            float reweight_a_norm1 = (8*cvl*cvq*cal*caq*G_F*G_F*pow(m_Z0,4)*s)/(256*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_a_norm2 = (2*cal*caq*(m_Z0*m_Z0-s)*alpha*G_F*m_Z0*m_Z0*Q_q)/(8*sqrt(2)*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_a_norm = (reweight_a_norm1 - reweight_a_norm2)*n_conv*LQ_jacobian;
            reweight_a = reweight_a_norm*gen_cost/LQ_denom;

            float reweight_dy = reweight_s + reweight_a;
            if(!old){
                //using h_sym for full dy temp to check new method
            h_sym->Fill(tm.m, var1, tm.cost, reweight_dy * tm.evt_weight *tm.evt_pdfweight); 
           

            h_asym->Fill(tm.m, var1, tm.cost, reweight_a * tm.evt_weight*tm.evt_pdfweight );
            h_asym->Fill(tm.m, var1, -tm.cost, -reweight_a * tm.evt_weight*tm.evt_pdfweight );
            }
            

            */






              //for LQ, 2 terms-> pure and interference

              // LQ terms: LO LQ/LO SM
            // Need to modify LQ_denom
              //flag_q=1 for d-dbar, 2 for u-ubar, 3 for s-sbar, 4 for c-cbar, 0 for everything
        int flag_q=0;
       // if((tm.inc_id1 == 1 && tm.inc_id2 == -1)||(tm.inc_id1 == -1 && tm.inc_id2 == 1)) flag_q=1;
       // if((tm.inc_id1 == 2 && tm.inc_id2 == -2)||(tm.inc_id1 == -2 && tm.inc_id2 == 2)) flag_q=2;
        if((tm.inc_id1 == 3 && tm.inc_id2 == -3)||(tm.inc_id1 == -3 && tm.inc_id2 == 3)) flag_q=1;
        if((tm.inc_id1 == 4 && tm.inc_id2 == -4)||(tm.inc_id1 == -4 && tm.inc_id2 == 4)) flag_q=2;
        if(flag_q!=0){ 

          if(flag_q==1){
            Q_q=Q_d;
            caq=caq_d;
            cvq=cvq_d;
          }

          if(flag_q==2){
            Q_q=Q_u;
            caq=caq_u;
            cvq=cvq_u;
          }

         
            float LQ_denom = get_LQ_denom(gen_cost,s,Q_q,caq,cvq);

              //float reweight_LQpure_norm = (n_conv*LQ_jacobian/(128*M_PI*s));
            float reweight_LQpure_norm = (1/(128*M_PI*s));
              //weight(cost)
            float reweight_LQpure_num1 = ((1 - gen_cost)*(1 - gen_cost));
            float reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1-gen_cost)* ((2*m_LQ*m_LQ/s)+1-gen_cost));
            float reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            float reweight_LQpure_pos;
            reweight_LQpure_pos = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);
              //weight(-cost)
            reweight_LQpure_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1+gen_cost)* ((2*m_LQ*m_LQ/s)+1+gen_cost));
            reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            float reweight_LQpure_neg;
            reweight_LQpure_neg = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);



            float reweight_LQint_norm1 = ((alpha*Q_q)/(16*s));
            float reweight_LQint_norm2_num = ((m_Z0*m_Z0-s)*(cal+cvl)*(caq-cvq)*G_F*m_Z0*m_Z0);
            float reweight_LQint_norm2_denom = (128*1.4142*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_LQint_norm2 = (reweight_LQint_norm2_num/reweight_LQint_norm2_denom);
             // float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2)*n_conv*LQ_jacobian;
            float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2);
               //weight(cost)
            float reweight_LQint_num1 = ((1 - gen_cost)*(1 - gen_cost));
            float reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1-gen_cost);
            float reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            float reweight_LQint_pos;
            reweight_LQint_pos = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);
              //weight(-cost)
            reweight_LQint_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1+gen_cost);
            reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            float reweight_LQint_neg;
            reweight_LQint_neg = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);

              //printf("%f\n",reweight_LQpure_pos);

              //dLQ temps
            if(flag_q==1){
              h_LQpure_d->Fill(tm.m, var1, tm.cost, reweight_LQpure_pos * tm.evt_weight ); 
              //h_LQpure_d->Fill(tm.m, var1, -tm.cost, reweight_LQpure_neg * tm.evt_weight );
              h_LQint_d->Fill(tm.m, var1, tm.cost, reweight_LQint_pos * tm.evt_weight ); 
              //h_LQint_d->Fill(tm.m, var1, -tm.cost, reweight_LQint_neg * tm.evt_weight);
            }
              //uLQ temps
            if(flag_q==2){
              h_LQpure_u->Fill(tm.m, var1, tm.cost, reweight_LQpure_pos * tm.evt_weight ); 
              //h_LQpure_u->Fill(tm.m, var1, -tm.cost, reweight_LQpure_neg * tm.evt_weight);
              h_LQint_u->Fill(tm.m, var1, tm.cost, reweight_LQint_pos * tm.evt_weight); 
              //h_LQint_u->Fill(tm.m, var1, -tm.cost, reweight_LQint_neg * tm.evt_weight);
            }

          }
        }
      }
      tm.finish();
    //int mbin = find_bin(m_bins, m_low + 0.1);
   // tm.fixRFNorm(h_sym, mbin); //not done for LQ
    //tm.fixRFNorm(h_asym, mbin);
   // tm.fixRFNorm(h_alpha, mbin,year);


      h_sym->Scale(0.5);
      h_asym->Scale(0.5);
      h_alpha->Scale(0.5);
      //h_LQpure_u->Scale(0.5);
      //h_LQint_u->Scale(0.5);
      //h_LQpure_d->Scale(0.5);
      //h_LQint_d->Scale(0.5);

    //cleanup_template(h_sym);
      fixup_template_sum(h_sym, h_asym);
      t1->ResetBranchAddresses();
      printf("MC templates generated from %i events. Sym integral is %.1f \n \n", n, h_sym->Integral()); 

      return 0;
    }

    int gen_combined_background_template(int nTrees, TTree **ts, TH3F* h,  
      int year,  int flag1 = FLAG_MUONS, 
      bool ss =false, bool use_xF = false,  bool emu_reweight = false, const string &sys_label = ""){
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

         auto hname = h->GetName();
        if(string(hname).find("top") != string::npos) tm.btag_mc_eff_idx = 0; //idx for ttbar MC btag effs
        else if(string(hname).find("gam") != string::npos) tm.btag_mc_eff_idx = 1; //idx for DY MC btag effs
        else tm.btag_mc_eff_idx = 2; //idx for diboson MC btag effs

        printf("h %s : btag mc idx %i \n", h->GetName(), tm.btag_mc_eff_idx);

        tm.setup_systematic(sys_label);
        tm.setup();


        for (int i=0; i<tm.nEntries; i++) {
          tm.getEvent(i);
          bool pass =  tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

          if(pass){
            tm.doCorrections();
            tm.getEvtWeight();

            float var1 = abs(tm.cm.Rapidity());
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

    int one_mc_template(TTree *t1, float a0, float afb, TH3F* h_dy, 
      int year, float m_LQ, int flag1 = FLAG_MUONS, bool use_xF = false, 
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
    // uLQ
      TH3F h_LQpure_u = TH3F("h_LQpure_u", "LQpure template of mc",
        n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
      h_LQpure_u.SetDirectory(0);
      TH3F h_LQint_u = TH3F("h_LQint_u", "LQint template of mc",
        n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
      h_LQint_u.SetDirectory(0);
    //dLQ
      TH3F h_LQpure_d = TH3F("h_LQpure_d", "LQpure template of mc",
        n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
      h_LQpure_d.SetDirectory(0);
      TH3F h_LQint_d = TH3F("h_LQint_d", "LQint template of mc",
        n_lq_m_bins, lq_m_bins, n_var1_bins, var1_bins, n_cost_bins, cost_bins);
      h_LQint_d.SetDirectory(0);
    //includes m_LQ
    gen_mc_template(t1, &h_sym, &h_asym, &h_alpha, &h_LQpure_u, &h_LQint_u, &h_LQpure_d, &h_LQint_d, year, m_LQ, flag1,  use_xF,sys_label);


        float alpha = 2.* a0/ (2. - a0);
    double norm = 3./4./(2.+alpha);

    TH3F *h_pl = (TH3F *) h_sym.Clone("h_pl");
    TH3F *h_mn = (TH3F *) h_sym.Clone("h_pl");
    h_pl->Reset();
    h_mn->Reset();
    make_pl_mn_templates(&h_sym, &h_asym, h_pl, h_mn);


    

    h_dy->Add(h_pl, h_mn, (norm + afb), (norm - afb));
    h_alpha.Scale(norm * alpha);
    h_dy->Add(&h_alpha);
    //fix this if running syscheck
      h_dy->Add(&h_LQpure_u);
      h_dy->Add(&h_LQint_u);

      return 1;
    }




void gen_fakes_template(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree* t_QCD_contam, TH3F *h, 
        int year, int flag1 = FLAG_MUONS, bool incl_ss = true, bool ss_binning = false, bool use_xF = false, string sys_label = ""){
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
   
    float eta1,eta2,pt1,pt2;


    TH2D *h_err;
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    float pt;
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

            bool pass = (tm.m >= lq_m_bins[0]) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

            bool opp_sign = false;
            if(flag1 == FLAG_MUONS) opp_sign = ((abs(tm.mu1_charge - tm.mu2_charge)) > 0.01);
            else opp_sign = ((abs(tm.el1_charge - tm.el2_charge)) > 0.01);

            if(!incl_ss) pass = pass && opp_sign; 
            if(pass){
                double evt_reweight = 0.;

                if(flag1 == FLAG_MUONS){
                    float mu1_fakerate, mu2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);

                        evt_reweight = mu1_fakerate/(1-mu1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 1);
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), 1);
                    }

                    if(l==1){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));

                        h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f),  -0.5*(mu2_fakerate/(1-mu2_fakerate)));
                        h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f),  -0.5*(mu1_fakerate/(1-mu1_fakerate)));
                    }
                    if(l==2){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate )/(1-mu1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), -tm.getEvtWeight());
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), -tm.getEvtWeight());
                    }
                    if(l==3){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                        h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 0.5*tm.getEvtWeight() * (mu2_fakerate/(1-mu2_fakerate)));
                        h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), 0.5*tm.getEvtWeight() * (mu1_fakerate/(1-mu1_fakerate)));
                    }
                }
                else{
                    float el1_fakerate, el2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = el1_fakerate/(1-el1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), 1);
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), 1);
                    }
                    if(l==1){
                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));

                        h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), -0.5*(el2_fakerate/(1-el2_fakerate)));
                        h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), -0.5*(el1_fakerate/(1-el1_fakerate)));

                    }
                    if(l==2){

                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate )/(1-el1_fakerate);


                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), -tm.getEvtWeight());
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), -tm.getEvtWeight());
                    }
                    if(l==3){

                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));

                        h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f),  0.5*tm.getEvtWeight() * (el2_fakerate/(1-el2_fakerate)));
                        h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f),  0.5*tm.getEvtWeight() * (el1_fakerate/(1-el1_fakerate)));
                    }
                }
                double tot_evt_weight = 0.;
                if(is_data) tot_evt_weight = evt_reweight; 
                else{
                    tot_evt_weight = evt_reweight * tm.getEvtWeight();
                }

                float var1 = abs(tm.cm.Rapidity());
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


    //symmetrize to avg statistical fluctuations before removing negative bins
    symmetrize3d(h);
    //remove negative bins
    cleanup_fakes_template(h);


    //scale to just OS normalization
    //
    printf("Tot weight OS is %.1f, ss %.1f, hist integral %.1f \n", tot_weight_os, tot_weight_ss, h->Integral());
    if(tot_weight_os > 10.){
        tot_weight_ss = std::max(1e-4, tot_weight_ss);
        tot_weight_os = std::max(1e-4, tot_weight_os);
    }
    else{ //low stats, use weight after negative bin cleanups
        tot_weight_ss = std::max(1e-4, tot_weight_ss);
        tot_weight_os = h->Integral() - tot_weight_ss;
    }
    float scaling = tot_weight_os / (tot_weight_ss + tot_weight_os);
    if(!incl_ss) scaling = 1.;
    h->Scale(scaling);
    h_err->Scale(scaling);

    set_fakerate_errors(h_err, FR.h, h);


    fakes_costrw_helper h_rw;
    setup_fakes_costrw_helper(&h_rw, year);
    if(flag1 == FLAG_MUONS) fakes_cost_reweight(h, h_rw.mu_rw, shape_sys);
    else fakes_cost_reweight(h, h_rw.el_rw, shape_sys);
    
    Double_t err;
    Double_t integ = h->IntegralAndError(1, h->GetNbinsX(), 1, h->GetNbinsY(), 1, h->GetNbinsZ(), err);


    //printf("After Scale: \n");
    //h->Print("range");
    printf("Total fakerate est is %.0f +/- %.0f \n", integ, err);
    return;
}

//make templates based on generator level samples (used for assessing impact of
//fiducial cuts on AFB
float make_gen_temps(TTree *t_gen, TH3F *h_raw, TH3F *h_sym, TH3F *h_asym, TH3F *h_alpha,  TH3F *h_LQpure_u, TH3F *h_LQpure_d,  TH3F *h_LQint_u, TH3F *h_LQint_d,
        float m_LQ, bool do_ptrw = false, int year = 2016, string sys_label = ""){

    printf("Making LQ+DY generator level templates\n");
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    int inc_id1, inc_id2;
    float mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_up, mu_RF_down;
    float evt_weight;
    float pdf_weights[60];
    int  do_ptrw_sys = 0;
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
    t_gen->SetBranchAddress("mu_R_up", &mu_R_up);
    t_gen->SetBranchAddress("mu_R_down", &mu_R_down);
    t_gen->SetBranchAddress("mu_F_up", &mu_F_up);
    t_gen->SetBranchAddress("mu_F_down", &mu_F_down);
    t_gen->SetBranchAddress("mu_RF_up", &mu_RF_up);
    t_gen->SetBranchAddress("mu_RF_down", &mu_RF_down);
    //t_gen->SetBranchAddress("pdf_weights", &pdf_weights);
    t_gen->SetBranchAddress("inc_id1", &inc_id1);
    t_gen->SetBranchAddress("inc_id2", &inc_id2);

    A0_helpers A0_helper; 
    setup_A0_helper(&A0_helper, year);
    /*
    ptrw_helper ptrw_SFs; 
    setup_ptrw_helper(&ptrw_SFs, year);


    if(sys_label.find("ptrw") != string::npos){
        int foo;
        if(  sys_label.find("Up") != string::npos){
            sscanf(sys_label.c_str(), "ptrw%ibUp", &do_ptrw_sys);
        }
        else{
            sscanf(sys_label.c_str(), "ptrw%ibDown", &do_ptrw_sys);
            do_ptrw_sys *= -1;
        }
    printf("ptrw sys %i \n", do_ptrw_sys);
    }

  */
    //float pt_cut = 26.;
    float pt_cut = 30.;
    float sum_weights = 0.;

    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
      //  if(sys_label.find("ptcutUp") != string::npos) pt_cut = 38.;
      //  if(sys_label.find("ptcutdown") != string::npos) pt_cut = 26.;
            
        
            evt_weight = gen_weight;
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            float rap = abs(cm.Rapidity());

            bool pass = gen_weight>0. && m > lq_m_bins[0] && abs(gen_lep_p->Eta()) < 2.5 && abs(gen_lep_m->Eta()) < 2.5 && min(gen_lep_m->Pt(), gen_lep_p->Pt()) > 15.;

          //  && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > pt_cut && min(gen_lep_m->Pt(), gen_lep_p->Pt()) > 15.;
        //bool pass = abs(cm.Rapidity()) < 2.4;
            /*
            if(sys_label == string("RENORMUp")) evt_weight *= mu_R_up;
            else if(sys_label == string ("RENORMDown")) evt_weight *= mu_R_down;
            else if(sys_label == string ("FACUp")) evt_weight *= mu_F_down;
            else if(sys_label == string ("FACDown")) evt_weight *= mu_F_down;
            else if(sys_label == string ("REFACUp")) evt_weight *= mu_RF_down;
            else if(sys_label == string ("REFACDown"))  evt_weight *= mu_RF_down;
            else if(sys_label.find("pdf") != string::npos){
                int n_pdf = 60;
                float comb = 0.;
                for(int j=0; j<n_pdf; j++){
                    float diff = abs(1 - pdf_weights[j]);
                    if(diff >= 1.) diff = 0.9;
                    comb += pow(diff,2);
                }
                comb = sqrt(comb);
                if(sys_label.find("Up") != string::npos) evt_weight *= abs((1+comb));
                if(sys_label.find("Down") != string::npos) evt_weight *= abs((1-comb));
            }
            else if (sys_label.find("ptrw") == string::npos && sys_label.find("ptcut") == string::npos && sys_label != string("")){
                printf("Can't find sys %s \n", sys_label.c_str());
            }

            if(do_ptrw){
                float ptrw = get_ptrw_SF(ptrw_SFs, m, pt, do_ptrw_sys); 
                evt_weight *= ptrw;
            }


            */
           // h_uncut->Fill(cost_st, evt_weight);
            if(pass){


                float gen_cost = cost_st;
                float s = m*m;

                if(evt_weight >0) nEvents++;
                else if(evt_weight<0.)  {
                  nEvents--;
                  printf("gen_weight in templates is negative = %f\n",evt_weight);
                }

                sum_weights+=gen_weight;
                

                float denom = get_reweighting_denom(A0_helper, gen_cost, m, pt, rap);

                float reweight_a = gen_cost/ denom;
                float reweight_s = (1 + gen_cost*gen_cost)/denom;
                float reweight_alpha = (1 - gen_cost*gen_cost)/denom;




                h_raw->Fill(m, rap, gen_cost, evt_weight);

                h_sym->Fill(m, rap, gen_cost, reweight_s * evt_weight *5000); 
                h_sym->Fill(m, rap, -gen_cost, reweight_s * evt_weight *5000); 

                h_asym->Fill(m, rap, gen_cost, reweight_a * evt_weight *5000);
                h_asym->Fill(m, rap, -gen_cost, -reweight_a * evt_weight *5000);

                h_alpha->Fill(m, rap, gen_cost, reweight_alpha * evt_weight *5000); 
                h_alpha->Fill(m, rap, -gen_cost, reweight_alpha * evt_weight *5000); 

                int flag_q=0;
                if((inc_id1 == 1 && inc_id2 == -1)||(inc_id1 == -1 && inc_id2 == 1)) flag_q=1;
                if((inc_id1 == 2 && inc_id2 == -2)||(inc_id1 == -2 && inc_id2 == 2)) flag_q=2;
                if(flag_q!=0){ 

          if(flag_q==1){
            Q_q=Q_d;
            caq=caq_d;
            cvq=cvq_d;
          }

          if(flag_q==2){
            Q_q=Q_u;
            caq=caq_u;
            cvq=cvq_u;
          }

          
            float LQ_denom = get_LQ_denom(gen_cost,s,Q_q,caq,cvq);


              //float reweight_LQpure_norm = (n_conv*LQ_jacobian/(128*M_PI*s));
            float reweight_LQpure_norm = (1/(128*M_PI*s));
              //weight(cost)
            float reweight_LQpure_num1 = ((1 - gen_cost)*(1 - gen_cost));
            float reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1-gen_cost)* ((2*m_LQ*m_LQ/s)+1-gen_cost));
            float reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            float reweight_LQpure_pos;
            reweight_LQpure_pos = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);
              //weight(-cost)
            reweight_LQpure_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1+gen_cost)* ((2*m_LQ*m_LQ/s)+1+gen_cost));
            reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
            float reweight_LQpure_neg;
            reweight_LQpure_neg = (reweight_LQpure_norm*reweight_LQpure_num/LQ_denom);



            float reweight_LQint_norm1 = ((alpha*Q_q)/(16*s));
            float reweight_LQint_norm2_num = ((m_Z0*m_Z0-s)*(cal+cvl)*(caq-cvq)*G_F*m_Z0*m_Z0);
            float reweight_LQint_norm2_denom = (128*1.4142*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
            float reweight_LQint_norm2 = (reweight_LQint_norm2_num/reweight_LQint_norm2_denom);
             // float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2)*n_conv*LQ_jacobian;
            float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2);
               //weight(cost)
            float reweight_LQint_num1 = ((1 - gen_cost)*(1 - gen_cost));
            float reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1-gen_cost);
            float reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            float reweight_LQint_pos;
            reweight_LQint_pos = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);
              //weight(-cost)
            reweight_LQint_num1 = ((1 + gen_cost)*(1 + gen_cost));
            reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1+gen_cost);
            reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
            float reweight_LQint_neg;
            reweight_LQint_neg = (reweight_LQint_norm*reweight_LQint_num/LQ_denom);

            if(flag_q==1){
              h_LQpure_d->Fill(m, rap, gen_cost, reweight_LQpure_pos * evt_weight * 5000 ); 
              //h_LQpure_d->Fill(tm.m, var1, -tm.cost, reweight_LQpure_neg * tm.evt_weight );
              h_LQint_d->Fill(m, rap, gen_cost, reweight_LQint_pos * evt_weight * 5000); 
              //h_LQint_d->Fill(tm.m, var1, -tm.cost, reweight_LQint_neg * tm.evt_weight);
            }
              //uLQ temps
            if(flag_q==2){
              h_LQpure_u->Fill(m, rap, gen_cost, reweight_LQpure_pos * evt_weight * 5000); 
              //h_LQpure_u->Fill(tm.m, var1, -tm.cost, reweight_LQpure_neg * tm.evt_weight);
              h_LQint_u->Fill(m, rap, gen_cost, reweight_LQint_pos * evt_weight * 5000); 
              //h_LQint_u->Fill(tm.m, var1, -tm.cost, reweight_LQint_neg * tm.evt_weight);
            }
            }
	}        
    }
    printf("selected %i events; total evt_weight = %f \n", nEvents,sum_weights);

    //cleanup_template(h_sym);
    //fixup_template_sum(h_sym, h_asym);

    return sum_weights;

}

int make_gen_data_temps(TTree *t_gen, TH3F *h_data, int year = 2016, float sum_weights = 1.){

    printf("Making data generator level templates\n");
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    int inc_id1, inc_id2;
   
    float evt_weight;
  
    Bool_t sig_event(1);

    t_gen->SetBranchAddress("lep_pls", &gen_lep_p);
    t_gen->SetBranchAddress("lep_mns", &gen_lep_m);
    //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
    //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
    //t_gen->SetBranchAddress("m", &m);
    //t_gen->SetBranchAddress("cost", &cost);
    //t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);
   
    //t_gen->SetBranchAddress("pdf_weights", &pdf_weights);
  //  t_gen->SetBranchAddress("inc_id1", &inc_id1);
   // t_gen->SetBranchAddress("inc_id2", &inc_id2);

    A0_helpers A0_helper; 
    setup_A0_helper(&A0_helper, year);
    
    //float pt_cut = 26.;
    float pt_cut = 30.;


    int nEvents=0;
    float sum_weights_data = 0.;

    for(int i=0; i<t_gen->GetEntries();i++){
      t_gen->GetEntry(i);
	cm = *gen_lep_p + *gen_lep_m;
            float m = cm.M();
      if(m > lq_m_bins[0] && gen_weight>0.)
      sum_weights_data+=gen_weight;
      else if(gen_weight<0.) printf("gen_weight is negative = %f\n",gen_weight);
    }
    printf("sum_weights_data = %f\n",sum_weights_data);
    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
    

        
            evt_weight = gen_weight;
            cm = *gen_lep_p + *gen_lep_m;
             m = cm.M();
            float pt = cm.Pt();
            float rap = abs(cm.Rapidity());
            float gen_cost = get_cost(*gen_lep_p, *gen_lep_m, false);

            bool pass = m > lq_m_bins[0]; //&& abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && min(gen_lep_m->Pt(), gen_lep_p->Pt()) > 15.;
         //   && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > pt_cut 
        //bool pass = abs(cm.Rapidity()) < 2.4;
        

            if(pass){

              if(evt_weight >0) nEvents++;
              else  nEvents--;

              h_data->Fill(m, rap, gen_cost, evt_weight*5000*(sum_weights/sum_weights_data));
            }
          }
          printf("selected %i events \n", nEvents);

    //cleanup_template(h_sym);
    //fixup_template_sum(h_sym, h_asym);

    return nEvents;

}

