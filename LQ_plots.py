import pylab
import numpy as np
from math import *
from scipy.integrate import quad

#define constants

alpha = 1/127.9
m_Z0 = 91.1875
sin2_thetaw = 0.231 #sin^2(theta_W) (weak mixing angle)
G_F = 1.166e-5
g_z = 2.4952 #width of Z0

#checking ssh

def constants(flag):
    #use coupling definitions from Quigg edition 1
    #up quark
    Q_q = 2./3. 
    I_3 = 1./2.
    if flag== 1:
        #down quark
        Q_q = -1./3.
        I_3 = -1./2.

    crq = -2 *Q_q * sin2_thetaw
    clq = 2*I_3- 2. *Q_q * sin2_thetaw

    crl = 2 * sin2_thetaw
    cll = 2 * sin2_thetaw - 1

    cvq = crq + clq
    caq = -(clq - crq)

    cvl = crl + cll
    cal = crl - cll

    return cvl,cal,cvq,caq,Q_q




def LQ_cost(flag, cost, s, m_lq, prnt = False):
        cvl,cal,cvq,caq,Q_q = constants(flag)
        SM_part = SM_cost(flag,cost, s)

        #pure LQ term
        XS3 = ((cost-1)**2*s*y_lq**4)/(128*pi*(2*m_lq**2+ s*(1-cost))**2)
        #LQ gamma interference
        XS67 = (alpha*Q_q*y_lq**2*(cost-1)**2)/(16*(2*m_lq**2+s*(1-cost)))
        #LQ Z0 interference
        XS89_num = (G_F*m_Z0**2*s*y_lq**2*(cal+cvl)*(caq-cvq)*(cost-1)**2*(s-m_Z0**2))
        XS89_denom = (128*sqrt(2)*pi*(2*m_lq**2+s*(1-cost))*((m_Z0**2-s)**2+g_z**2*m_Z0**2))
        XS89 = XS89_num / XS89_denom
        #if(prnt): print("LQ terms: ", SM_part, XS3, XS67, XS89)
        
        return SM_part + XS3 + XS67 + XS89
        #return XS3 + XS67 + XS89

def LQ_vec_cost(flag, cost, s, prnt = False):
        cvl,cal,cvq,caq,Q_q = constants(flag)
        SM_part = SM_cost(flag,cost, s)

        #pure LQ term
        XS3 = ((cost+1)**2*s*g_lq**4)/(32*pi*(2*m_lq**2+ s*(1-cost))**2)
        #LQ gamma interference
        XS67 = (alpha*Q_q*g_lq**2*(cost+1)**2)/(8*(2*m_lq**2+s*(1-cost)))
        #LQ Z0 interference
        XS89_num = -(G_F*m_Z0**2*s*g_lq**2*(cal-cvl)*(caq-cvq)*(cost+1)**2*(s-m_Z0**2))
        XS89_denom = (64*sqrt(2)*pi*(2*m_lq**2+s*(1-cost))*((m_Z0**2-s)**2+g_z**2*m_Z0**2))
        XS89 = XS89_num / XS89_denom
        #if(prnt): print("LQ terms: ", SM_part, XS3, XS67, XS89)
        
        return SM_part + XS3 + XS67 + XS89


def SM_cost(flag, cost, s, prnt = False):
        cvl,cal,cvq,caq,Q_q = constants(flag)
        #pure gamma term
        #cost = -cost # there is some weird negative sign issue 
        XS1 = (pi*alpha**2*Q_q**2*(cost**2+1))/(2*s)
        #pure Z0 term
        XS2_num = ((((cal*caq*cost**2+ cal*caq+ 8*cost*cvl*cvq)*caq +(cost**2+1)*cal*cvq**2)*cal+(caq**2+cvq**2)*(cost**2+1)*cvl**2)*G_F**2*m_Z0**4*s)
        XS2_denom = (256*pi*((m_Z0**2-s)**2 + g_z**2*m_Z0**2))
        XS2 = XS2_num/ XS2_denom
        #Z0 gamma interference
        XS45_num =  - ((cost**2+1)*cvl*cvq + 2*cal*caq*cost) * (s-m_Z0**2) * alpha*G_F*m_Z0**2*Q_q
        XS45_denom = (8*sqrt(2)*((m_Z0**2-s)**2+(g_z*m_Z0)**2))
        XS45 = XS45_num/XS45_denom
        if(prnt): print("SM Terms: ", XS1, XS2, XS45)
        
        return XS1 + XS2 + XS45

'''

#rough form of standard model distribution
def SM(cost):
        norm = 3./8.
        return (1. + cost**2) *norm + 0.6*cost

def tot_xsec(s):
        return quad(lambda x: SM_cost(flag,x, s), -1., 1.)[0]

def AFB(s):
    SM_for = quad(lambda x: SM_cost(flag,x,s), 0., 1.)[0]
    SM_back = quad(lambda x: SM_cost(flag,x,s), -1., 0.)[0]
    return (SM_for - SM_back)/(SM_for + SM_back)




E_range = numpy.linspace(1., 400., 400)
xsecs = []
Afbs = []
for E in E_range:
    xs = tot_xsec(E**2)
    afb = AFB(E**2)
    xsecs.append(xs)
    Afbs.append(afb)

pylab.plot(E_range, xsecs, 'b')
pylab.yscale("log")
pylab.show()
pylab.close()

pylab.plot(E_range, Afbs, 'b')
pylab.show()
pylab.close()
'''
x_axis = np.linspace(-1,1,1000) # 1000 linearly spaced numbers


#SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
#sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#pylab.plot(x_axis,sm_v2,'b',label='SM')
'''
y_lq_list = np.linspace(-1.,1.,60)
y_lqsq_list = [x*x for x in y_lq_list]
m_ll_list = np.linspace(500,5000,300)
cost_list = np.linspace(-1,1,300)
flag = 1

for mLQ in [500,1000,2000,3000]:
    print("======= mLQ = %i ======="%mLQ)
    avg_xsec, avg_xsec_err = [],[]
    for y_lq in y_lq_list:
        f_list = []
	#y_lq = sqrt(abs(y_lq2))
	print('doing y_lq = %.3f'%y_lq)
        for m_ll in m_ll_list:
            for cost in cost_list:
                

                s = m_ll*m_ll
                f_list.append(LQ_cost(flag,cost,s,mLQ))

        avg_xsec.append(np.mean(f_list)/1e-11)
        avg_xsec_err.append(np.std(f_list)/1e-11)
    print(np.amin(avg_xsec), np.amax(avg_xsec))
    pylab.plot(y_lq_list, avg_xsec, label = 'mLQ=%i GeV'%mLQ)

pylab.ylim(0.5,35)
pylab.title('LQ angular distribution averaged over s and cost')
pylab.ylabel('(LQ xsec averaged over s and cost)/1e-11')
pylab.xlabel('yLQ (mu-d)')
pylab.legend()
pylab.savefig("xsec_vs_yLQ_dm.png")

'''

for m_ll in [500,600,700,800,900,1000,1500,2000]:

    flag = 1
    s = m_ll*m_ll

    SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
    sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
   # sm = SM(x_axis)
    #pylab.plot(x_axis,sm_v2,label='SM')
#   pylab.plot(x_axis,sm,'.',label=r'SM (rough),$m_{ll}=$'+str(m_ll))

    m_lq = 2000.
   

    y_lq = .6
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,label=r'$m_{ll}=$'+str(m_ll))

   #  g_lq = 1.6
   #  LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
   #  F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # # print F
   #  pylab.plot(x_axis,F,'black',label=r'$g_{ue}=1.6$')


   #  g_lq = 3.0
   #  LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
   #  F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # # print F
   #  pylab.plot(x_axis,F,'g',label=r'$g_{ue}=3.0$')


    

pylab.plot(x_axis,sm_v2,label='SM')
pylab.title(r'2 TeV $S_{ed}$, $y_{ed}$ = '+str(y_lq) )
pylab.xlabel(r'$cos \theta$')
pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
pylab.ylim([0., 1.4])
pylab.legend()
pylab.savefig("Sed_mll.png")
pylab.close()

for m_ll in [500]:

    flag = 2
    s = m_ll*m_ll

    SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
    sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
   # sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#   pylab.plot(x_axis,sm,'.',label='SM (rough)')

    m_lq = 2000.
   

    g_lq = .6
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$g_{ue}=0.6$')

    g_lq = 1.6
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$g_{ue}=1.6$')


    g_lq = 3.0
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'g',label=r'$g_{ue}=3.0$')


    


    pylab.title(r"ElectroUp, $m_{ee}$ = "+str(m_ll)+ r" GeV, $m_{LQ}$ ="+(" %.1f TeV" % (m_lq/1000)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.savefig("vLQ_ElectroUp_gLQ.png")
    pylab.close()

    SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
    sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    #sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#pylab.plot(x_axis,sm,'.',label='SM (rough)')

    g_lq = 2
    m_lq = 1000.
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$m_{LQ}=1$ TeV')


    m_lq = 2000.
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$m_{LQ}=2$ TeV')


    m_lq = 3000.
    LQ_norm = quad(lambda x: LQ_vec_cost(flag,x,s), -1., 1.)[0]
    F = LQ_vec_cost(flag,x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'g',label=r'$m_{LQ}=3$ TeV')


    pylab.title(r"ElectroUp, $m_{ee}$ = "+str(m_ll)+ r" GeV, $g_{ue}$ = "+("%.1f " % (g_lq)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.savefig("vLQ_ElectroUp_mLQ.png")
    pylab.close()

    flag = 1

    SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
    sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    #sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#   pylab.plot(x_axis,sm,'.',label='SM (rough)')

    m_lq = 1000.
    y_lq = .8
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$y_{de}=0.8$')

    y_lq = 1.6
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$y_{de}=1.6$')


    y_lq = 2.4
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'g',label=r'$y_{de}=2.4$')




    pylab.title(r"ElectroDown, $m_{ee}$ = "+str(m_ll)+ r" GeV, $m_{LQ}$ ="+(" %.1f TeV" % (m_lq/1000)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.savefig("sLQ_ElectroUp_gLQ.png")
    pylab.close()

    SM_norm = quad(lambda x: SM_cost(flag,x,s), -1., 1.)[0]
    sm_v2 = SM_cost(flag,x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    #sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#pylab.plot(x_axis,sm,'.',label='SM (rough)')

    y_lq = 2
    m_lq = 1000.
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$m_{LQ}=2$ TeV')


    m_lq = 2000.
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$m_{LQ}=3$ TeV')


    m_lq = 3000.
    LQ_norm = quad(lambda x: LQ_cost(flag,x,s), -1., 1.)[0]
    F = LQ_cost(flag,x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'g',label=r'$m_{LQ}=4$ TeV')


    pylab.title(r"ElectroDown, $m_{ee}$ = "+str(m_ll)+ r" GeV, $y_{de}$ = "+("%.1f " % (y_lq)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.savefig("sLQ_ElectroUp_mLQ.png")
    pylab.close()


