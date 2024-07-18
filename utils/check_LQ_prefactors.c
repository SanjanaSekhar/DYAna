#include <cmath>
#include <stdio.h>
float alpha = 1/127.9;
float m_Z0 = 91.1875;
float m_W = 80.379;
float sin2_thetaw = 0.231; //sin^2(theta_W) (weak mixing angle)
float G_F = 1.166e-5;
float g_z = 2.4952; //width of Z0

// testing running couplings
float b1 = 41./6.;
float b2 = -19./6.;
// g1 and g2 are SM electroweak couplings
float g1_Z = sqrt(8*G_F*m_Z0*m_Z0*sin2_thetaw/sqrt(2));
float g2_Z = sqrt((4*M_PI*alpha)/sin2_thetaw);



//up quark
float Q_u = 2./3. ;
float I3_u = 1./2. ;

//down quark
float Q_d = -1./3. ;
float I3_d = -1./2. ;


float alpha_run, g1_run, g2_run, G_F_run, sin2_thetaw_run, Q_q, caq, cvq, cvl, cal;

void set_running_couplings(float s, int quark_id){

	//printf("b1 = %f, b2 = %f\n",b1,b2);

	float m_ll = sqrt(s);

	 // testing running couplings

	
	g1_run = pow((1/(g1_Z*g1_Z))+((1/(8*M_PI*M_PI))*b1*log(m_Z0/m_ll)),-0.5);
	g2_run = pow((1/(g2_Z*g2_Z))+((1/(8*M_PI*M_PI))*b2*log(m_Z0/m_ll)),-0.5);


	sin2_thetaw_run = (g1_run*g1_run)/((g1_run*g1_run) + (g2_run*g2_run));
	alpha_run = (g2_run*g2_run)*sin2_thetaw_run/(4*M_PI);
	G_F_run = (sqrt(2)*(g1_run*g1_run))/(8*m_Z0*m_Z0*sin2_thetaw_run);
			//use coupling definitions from Quigg edition 1
	
	//sin2_thetaw_run = sin2_thetaw;
	//alpha_run =alpha;
	//G_F_run = G_F;

	float crl = 2 * sin2_thetaw_run;
	float cll = 2 * sin2_thetaw_run - 1;
	
		
	cvl = crl + cll;
	cal = crl - cll;
//up quark

	float crq_u = -2 *Q_u * sin2_thetaw_run;
	float clq_u = 2*I3_u- 2. *Q_u * sin2_thetaw_run;

	float cvq_u = crq_u + clq_u;
	float caq_u = -(clq_u - crq_u);

//down quark

	float crq_d = -2 *Q_d * sin2_thetaw_run;
	float clq_d = 2*I3_d- 2. *Q_d * sin2_thetaw_run;

	float cvq_d = crq_d + clq_d;
	float caq_d = -(clq_d - crq_d);

	if(quark_id==1){
		Q_q=Q_d;
		caq=caq_d;
		cvq=cvq_d;
	}

	if(quark_id==2){
		Q_q=Q_u;
		caq=caq_u;
		cvq=cvq_u;
	}
}
float get_LQ_denom(float gen_cost,float s,float Q_q, float caq, float cvq){
// if(test_sign)cal = crl - cll;



		float color_factor = 3.;
		float XS1 = (M_PI*pow(alpha_run,2)*pow(Q_q,2)*(pow(gen_cost,2)+1))/(2*color_factor*s);
//pure Z0 term
// float XS2_num = ((((cal*caq*pow(gen_cost,2)+ cal*caq+ 8*gen_cost*cvl*cvq)*caq +(pow(gen_cost,2)+1)*cal*pow(cvq,2))*cal+(pow(caq,2)+pow(cvq,2))*(pow(gen_cost,2)+1)*pow(cvl,2))*pow(G_F,2)*pow(m_Z0,4)*s);
		float XS2_num = (pow(cvl,2)+pow(cal,2))*(pow(caq,2)+pow(cvq,2))*(1+pow(gen_cost,2))*pow(G_F_run,2)*pow(m_Z0,4)*s;
		float XS2_denom = (256*color_factor*M_PI*(pow((m_Z0*m_Z0-s),2) + pow(g_z*m_Z0,2)));
		float XS2 = XS2_num/ XS2_denom;
//Z0 gamma interference
//float XS45_num =  - ((gen_cost*gen_cost+1)*cvl*cvq + 2*cal*caq*gen_cost) * (m_Z0*m_Z0-s) * alpha*G_F*m_Z0*m_Z0*Q_q;
		float XS45_num =  - ((gen_cost*gen_cost+1)*cvl*cvq) * (s- m_Z0*m_Z0) * alpha_run*G_F_run*m_Z0*m_Z0*Q_q;
		float XS45_denom = (8*color_factor*sqrt(2)*(pow((m_Z0*m_Z0-s),2)+pow((g_z*m_Z0),2)));
		float XS45 = XS45_num/XS45_denom;


//use events twice and make denom symmetric only
 //   if(only_sym){
 //     XS2_num = (pow(cvl,2)+pow(cal,2))*(pow(caq,2)+pow(cvq,2))*(1+pow(gen_cost,2))*pow(G_F,2)*pow(m_Z0,4)*s;
 //     XS2 = XS2_num/XS2_denom;
 //     XS45_num =  - ((gen_cost*gen_cost+1)*cvl*cvq) * (m_Z0*m_Z0-s) * alpha*G_F*m_Z0*m_Z0*Q_q;
 //     XS45 = XS45_num/XS45_denom;
 //   }
		//printf("XS1 = %f, XS2 = %f, XS45 = %f\n",XS1, XS2, XS45);
		float LQ_denom = (XS1 + XS2 + XS45); //new LQdenom is basically just LO SM
		return LQ_denom;
	}

float get_LQ_scalar_num(float gen_cost,float s,float Q_q, float caq, float cvq, float m_LQ, bool interference, bool negcos){
//	  if (test_sign) cal = crl - cll;
				 //float reweight_LQpure_norm = (n_conv*LQ_jacobian/(128*M_PI*s));


		float color_factor = 3.;
		float reweight_LQpure_norm = (1/(128*color_factor*M_PI*s));

		float reweight_LQint_norm1 = ((alpha_run*Q_q)/(16*color_factor*s));
		float  reweight_LQint_norm2_num = ((s - m_Z0*m_Z0)*(cal+cvl)*(caq-cvq)*m_Z0*m_Z0*G_F_run);
		
		float reweight_LQint_norm2_denom = (128*color_factor*1.4142*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));
		float reweight_LQint_norm2 = (reweight_LQint_norm2_num/reweight_LQint_norm2_denom);
						 // float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2)*n_conv*LQ_jacobian;
		float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2);
		float reweight_LQ_num, reweight_LQpure_num1, reweight_LQpure_denom1, reweight_LQpure_num, reweight_LQint_num1, reweight_LQint_denom1, reweight_LQint_num;
		
		printf("reweight_LQint_norm1 = %f, reweight_LQint_norm2 = %f\n",reweight_LQint_norm1*1e11, reweight_LQint_norm2*1e11);
		if(!interference and !negcos){
			
							//weight(cost)
			reweight_LQpure_num1 = ((1 - gen_cost)*(1 - gen_cost));
			reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1-gen_cost)* ((2*m_LQ*m_LQ/s)+1-gen_cost));
			reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
			reweight_LQ_num = reweight_LQpure_norm*reweight_LQpure_num;
		}
		else if(!interference and negcos){
							//weight(-cost)
			reweight_LQpure_num1 = ((1 + gen_cost)*(1 + gen_cost));
			reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1+gen_cost)* ((2*m_LQ*m_LQ/s)+1+gen_cost));
			reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
			reweight_LQ_num = reweight_LQpure_norm*reweight_LQpure_num;
		}
		else if(interference and !negcos){
					//weight(cost)
			reweight_LQint_num1 = ((1 - gen_cost)*(1 - gen_cost));
			reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1-gen_cost);
			reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
			reweight_LQ_num = reweight_LQint_norm*reweight_LQint_num;
		}
		else{
							//weight(-cost)
			reweight_LQint_num1 = ((1 + gen_cost)*(1 + gen_cost));
			reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1+gen_cost);
			reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
			reweight_LQ_num = reweight_LQint_norm*reweight_LQint_num;
		}
		return reweight_LQ_num;
	}

	float get_LQ_vec_num(float gen_cost,float s,float Q_q, float caq, float cvq, float m_LQ, bool interference, bool negcos){
//    if (test_sign) cal = crl - cll;
				 //float reweight_LQpure_norm = (n_conv*LQ_jacobian/(128*M_PI*s));
 	// testing running couplings


		float color_factor = 3.;
		float reweight_LQpure_norm = (1/(32*color_factor*M_PI*s));

		float reweight_LQint_norm1 = ((alpha_run*Q_q)/(8*color_factor*s));

		float reweight_LQint_norm2_num = -((s - m_Z0*m_Z0)*(cal-cvl)*(caq-cvq)*G_F_run*m_Z0*m_Z0);
		float reweight_LQint_norm2_denom = (64*color_factor*1.4142*M_PI*((m_Z0*m_Z0-s)*(m_Z0*m_Z0-s)+(g_z*g_z*m_Z0*m_Z0)));

		float reweight_LQint_norm2 = (reweight_LQint_norm2_num/reweight_LQint_norm2_denom);
						 // float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2)*n_conv*LQ_jacobian;
		float reweight_LQint_norm = (reweight_LQint_norm1 + reweight_LQint_norm2);
		float reweight_LQ_num, reweight_LQpure_num1, reweight_LQpure_denom1, reweight_LQpure_num, reweight_LQint_num1, reweight_LQint_denom1, reweight_LQint_num;

		if(!interference and !negcos){

							//weight(cost)
			reweight_LQpure_num1 = ((1 + gen_cost)*(1 + gen_cost));
			reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1-gen_cost)* ((2*m_LQ*m_LQ/s)+1-gen_cost));
			reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
			reweight_LQ_num = reweight_LQpure_norm*reweight_LQpure_num;
		}
		else if(!interference and negcos){
							//weight(-cost)
			reweight_LQpure_num1 = ((1 - gen_cost)*(1 - gen_cost));
			reweight_LQpure_denom1 = (((2*m_LQ*m_LQ/s)+1+gen_cost)* ((2*m_LQ*m_LQ/s)+1+gen_cost));
			reweight_LQpure_num =(reweight_LQpure_num1/reweight_LQpure_denom1);
			reweight_LQ_num = reweight_LQpure_norm*reweight_LQpure_num;
		}
		else if(interference and !negcos){
					//weight(cost)
			reweight_LQint_num1 = ((1 + gen_cost)*(1 + gen_cost));
			reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1-gen_cost);
			reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
			reweight_LQ_num = reweight_LQint_norm*reweight_LQint_num;
		}
		else{
							//weight(-cost)
			reweight_LQint_num1 = ((1 - gen_cost)*(1 - gen_cost));
			reweight_LQint_denom1 =  ((2*m_LQ*m_LQ/s)+1+gen_cost);
			reweight_LQint_num = (reweight_LQint_num1/reweight_LQint_denom1);
			reweight_LQ_num = reweight_LQint_norm*reweight_LQint_num;
		}
		return reweight_LQ_num;
	}

void check_LQ_prefactors(){
    float s = 600*600;
    float LQ_denom, reweight_LQpure_pos, reweight_LQint_pos;
    printf("For scalar up quark case\n");
    for(float c = -1.; c <= 1; c+=0.1){
        set_running_couplings(s,2);
        LQ_denom = get_LQ_denom(c, s, 2./3., -1., 0.384)*1e11;
        reweight_LQpure_pos = get_LQ_scalar_num(c, s, 2./3., -1., 0.384, 2500., false, false)*1e11;
        reweight_LQint_pos = get_LQ_scalar_num(c, s, 2./3., -1., 0.384, 2500., true, false)*1e11;
        printf("LQ_denom = %f, reweight_LQpure_pos = %f, reweight_LQint_pos = %f\n", LQ_denom, reweight_LQpure_pos, reweight_LQint_pos);
    }

    printf("For scalar down quark case\n");
    for(float c = -1.; c <= 1; c+=0.1){
        set_running_couplings(s,1);
        LQ_denom = get_LQ_denom(c, s, -1./3., 1., -0.692)*1e11;
        reweight_LQpure_pos = get_LQ_scalar_num(c, s, -1./3., 1., -0.692, 2500., false, false)*1e11;
        reweight_LQint_pos = get_LQ_scalar_num(c, s, -1./3., 1., -0.692, 2500., true, false)*1e11;
        printf("LQ_denom = %f, reweight_LQpure_pos = %f, reweight_LQint_pos = %f\n",LQ_denom, reweight_LQpure_pos, reweight_LQint_pos);
    }
}
