import pylab
import numpy
from math import *
from scipy.integrate import quad

#define constants

alpha = 1/127.9
m_Z0 = 91.1875
sin2_thetaw = 0.231 #sin^2(theta_W) (weak mixing angle)
G_F = 1.166e-5
g_z = 2.4952 #width of Z0

#use coupling definitions from Quigg edition 1

#up quark
Q_q = 2./3. 
I_3 = 1./2.

#down quark
#Q_q = -1./3.
#I_3 = -1./2.

crq = -2 *Q_q * sin2_thetaw
clq = 2*I_3- 2. *Q_q * sin2_thetaw

crl = 2 * sin2_thetaw
cll = 2 * sin2_thetaw - 1

cvq = crq + clq
caq = -(clq - crq)

cvl = crl + cll
cal = cll - crl

m_lq = 1000
y_lq = 1


def LQ_cost(cost,s, prnt = False):
        SM_part = SM_cost(cost, s)

        #pure LQ term
        XS3 = ((cost-1)**2*s*y_lq**4)/(128*pi*(2*m_lq**2+ s*(1-cost))**2)
        #LQ gamma interference
        XS67 = (alpha*Q_q*y_lq**2*(cost-1)**2)/(16*(2*m_lq**2+s*(1-cost)))
        #LQ Z0 interference
        XS89_num = (G_F*m_Z0**2*s*y_lq**2*(cal+cvl)*(caq-cvq)*(cost-1)**2*(m_Z0**2-s))
        XS89_denom = (128*sqrt(2)*pi*(2*m_lq**2+s*(1-cost))*((m_Z0**2-s)**2+g_z**2*m_Z0**2))
        XS89 = XS89_num / XS89_denom
        #if(prnt): print("LQ terms: ", SM_part, XS3, XS67, XS89)
        
        return SM_part + XS3 + XS67 + XS89
        #return XS3 + XS67 + XS89


def SM_cost(cost, s, prnt = False):
        #pure gamma term
        #cost = -cost # there is some weird negative sign issue 
        XS1 = (pi*alpha**2*Q_q**2*(cost**2+1))/(2*s)
        #pure Z0 term
        XS2_num = ((((cal*caq*cost**2+ cal*caq+ 8*cost*cvl*cvq)*caq +(cost**2+1)*cal*cvq**2)*cal+(caq**2+cvq**2)*(cost**2+1)*cvl**2)*G_F**2*m_Z0**4*s)
        XS2_denom = (256*pi*((m_Z0**2-s)**2 + g_z**2*m_Z0**2))
        XS2 = XS2_num/ XS2_denom
        #Z0 gamma interference
        XS45_num =  - ((cost**2+1)*cvl*cvq + 2*cal*caq*cost) * (m_Z0**2-s) * alpha*G_F*m_Z0**2*Q_q
        XS45_denom = (8*sqrt(2)*((m_Z0**2-s)**2+(g_z*m_Z0)**2))
        XS45 = XS45_num/XS45_denom
        if(prnt): print("SM Terms: ", XS1, XS2, XS45)
        
        return XS1 + XS2 + XS45



#rough form of standard model distribution
def SM(cost):
        norm = 3./8.
        return (1. + cost**2) *norm + 0.6*cost

def tot_xsec(s):
        return quad(lambda x: SM_cost(x, s), -1., 1.)[0]

def AFB(s):
    SM_for = quad(lambda x: SM_cost(x,s), 0., 1.)[0]
    SM_back = quad(lambda x: SM_cost(x,s), -1., 0.)[0]
    return (SM_for - SM_back)/(SM_for + SM_back)

x_axis = numpy.linspace(-1,1,1000) # 1000 linearly spaced numbers

'''
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

for m_ll in [1000]:


    s = m_ll*m_ll

    SM_norm = quad(lambda x: SM_cost(x,s), -1., 1.)[0]
    sm_v2 = SM_cost(x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#   pylab.plot(x_axis,sm,'.',label='SM (rough)')

    m_lq = 2000.
    y_lq = .8
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$y_{ue}=0.8$')

    y_lq = 1.8
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$y_{ue}=1.6$')


    y_lq = 2.6
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'g',label=r'$y_{ue}=2.4$')


    y_lq = 3.2
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'y',label=r'$y_{ue}=3.2$')



    pylab.title(r"ElectroUp, $m_{ee}$ = "+str(m_ll)+ r" GeV, $m_{LQ}$ ="+(" %.1f TeV" % (m_lq/1000)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.show()
    pylab.close()

    SM_norm = quad(lambda x: SM_cost(x,s), -1., 1.)[0]
    sm_v2 = SM_cost(x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#pylab.plot(x_axis,sm,'.',label='SM (rough)')

    y_lq = 2
    m_lq = 2000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$m_{LQ}=2$ TeV')


    m_lq = 3000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$m_{LQ}=3$ TeV')


    m_lq = 4000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'g',label=r'$m_{LQ}=4$ TeV')

    m_lq = 5000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'y',label=r'$m_{LQ}=5$ TeV')

    pylab.title(r"ElectroUp, $m_{ee}$ = "+str(m_ll)+ r" GeV, $y_{ue}$ = "+("%.1f " % (y_lq)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.show()
    pylab.close()

    #down quark
    Q_q = -1./3.
    I_3 = -1./2.


    SM_norm = quad(lambda x: SM_cost(x,s), -1., 1.)[0]
    sm_v2 = SM_cost(x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#   pylab.plot(x_axis,sm,'.',label='SM (rough)')

    m_lq = 2000.
    y_lq = .8
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$y_{ue}=0.8$')

    y_lq = 1.8
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$y_{ue}=1.6$')


    y_lq = 2.6
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'g',label=r'$y_{ue}=2.4$')


    y_lq = 3.2
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'y',label=r'$y_{ue}=3.2$')



    pylab.title(r"ElectroDown, $m_{ee}$ = "+str(m_ll)+ r" GeV, $m_{LQ}$ ="+(" %.1f TeV" % (m_lq/1000)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.show()
    pylab.close()

    SM_norm = quad(lambda x: SM_cost(x,s), -1., 1.)[0]
    sm_v2 = SM_cost(x_axis,s)/SM_norm
#afb = AFB(s)
#print("SM AFB is %.3f" % afb)
    sm = SM(x_axis)
    pylab.plot(x_axis,sm_v2,'b',label='SM')
#pylab.plot(x_axis,sm,'.',label='SM (rough)')

    y_lq = 2
    m_lq = 2000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'r',label=r'$m_{LQ}=2$ TeV')


    m_lq = 3000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
   # print F
    pylab.plot(x_axis,F,'black',label=r'$m_{LQ}=3$ TeV')


    m_lq = 4000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'g',label=r'$m_{LQ}=4$ TeV')

    m_lq = 5000.
    LQ_norm = quad(lambda x: LQ_cost(x,s), -1., 1.)[0]
    F = LQ_cost(x_axis,s)/LQ_norm
    #print F
    pylab.plot(x_axis,F,'y',label=r'$m_{LQ}=5$ TeV')

    pylab.title(r"ElectroDown, $m_{ee}$ = "+str(m_ll)+ r" GeV, $y_{ue}$ = "+("%.1f " % (y_lq)) )
    pylab.xlabel(r'$cos \theta$')
    pylab.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
    pylab.ylim([0., 1.4])
    pylab.legend()
    pylab.show()
    pylab.close()


