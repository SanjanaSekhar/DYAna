import pylab
import numpy as np
from math import *
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
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



def LQ_vec_cost(flag, cost, s, prnt = False):
        cvl,cal,cvq,caq,Q_q = constants(flag)
        SM_part = SM_cost(cost, s)

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

def LQ_cost(flag, cost, s, m_lq, y_lq, prnt = False):

        cvl,cal,cvq,caq,Q_q = constants(flag=1)

        A1 = (pi * alpha**2 * Q_q**2)/(6 * s)
        A2 = (G_F**2 * m_Z0**4 * s * (cal**2 + cvl**2) * (caq**2 + cvq**2))/(768 * pi * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A3 = -(alpha * Q_q * G_F * m_Z0**2 * (s-m_Z0**2) * cvl * cvq)/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A = A1 + A2 + A3

        B1 = (8 * G_F**2 * m_Z0**4 * s * cal * caq * cvl * cvq)/(768 * pi * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B2 = -(2 * alpha * Q_q * G_F * m_Z0**2 * cal * caq * (s - m_Z0**2))/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B = B1 + B2

        xsec1 = A * (1 + cost**2) + B * cost
        
        cvl,cal,cvq,caq,Q_q = constants(flag=2)

        A1 = (pi * alpha**2 * Q_q**2)/(6 * s)
        A2 = (G_F**2 * m_Z0**4 * s * (cal**2 + cvl**2) * (caq**2 + cvq**2))/(768 * pi * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A3 = -(alpha * Q_q * G_F * m_Z0**2 * (s-m_Z0**2) * cvl * cvq)/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A = A1 + A2 + A3

        B1 = (8 * G_F**2 * m_Z0**4 * s * cal * caq * cvl * cvq)/(768 * pi * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B2 = -(2 * alpha * Q_q * G_F * m_Z0**2 * cal * caq * (s - m_Z0**2))/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B = B1 + B2

        xsec2 = A * (1 + cost**2) + B * cost 
        SM_xsec = (1./3.) * xsec1 + (2./3.) * xsec2

        cvl,cal,cvq,caq,Q_q = constants(flag)

        # A1 = (pi * alpha**2 * Q_q**2)/(6 * s)
        # A2 = (G_F**2 * m_Z0**4 * s * (cal**2 + cvl**2) * (caq**2 + cvq**2))/(768 * pi * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        # A3 = -(alpha * Q_q * G_F * m_Z0**2 * (s-m_Z0**2) * cvl * cvq)/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        # A = A1 + A2 + A3

        # B1 = (8 * G_F**2 * m_Z0**4 * s * cal * caq * cvl * cvq)/(768 * pi * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        # B2 = -(2 * alpha * Q_q * G_F * m_Z0**2 * cal * caq * (s - m_Z0**2))/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        # B = B1 + B2

        C = 1/(384 * pi * s)
        D1 = (alpha * Q_q)/(48 * s)
        D2 = (G_F * m_Z0**2 * (cal + cvl) * (caq - cvq) * (s - m_Z0**2))/(384 * pi * 2**0.5 * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        D = D1 + D2
        if flag==1: C, D = 0.1875*C, 0.1875*D # d quark PDF is 1/4th of 75%
        else: C, D = 0.5625*C , 0.5625*D # u quark PDF is 3/4th of 75%


        return SM_xsec + C * y_lq**4 * ((1 - cost)/(1 - cost + (2 * m_lq**2)/s))**2 + D * y_lq**2 * (1 - cost)**2/(1 - cost + (2 * m_lq**2)/s)


def SM_cost(cost, s, prnt = False):
        
        '''
        xsec = A(1+cost^2)+B(cost)
        '''
        cvl,cal,cvq,caq,Q_q = constants(flag=1)

        A1 = (pi * alpha**2 * Q_q**2)/(6 * s)
        A2 = (G_F**2 * m_Z0**4 * s * (cal**2 + cvl**2) * (caq**2 + cvq**2))/(768 * pi * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A3 = -(alpha * Q_q * G_F * m_Z0**2 * (s-m_Z0**2) * cvl * cvq)/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A = A1 + A2 + A3

        B1 = (8 * G_F**2 * m_Z0**4 * s * cal * caq * cvl * cvq)/(768 * pi * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B2 = -(2 * alpha * Q_q * G_F * m_Z0**2 * cal * caq * (s - m_Z0**2))/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B = B1 + B2

        xsec1 = A * (1 + cost**2) + B * cost
        
        cvl,cal,cvq,caq,Q_q = constants(flag=2)

        A1 = (pi * alpha**2 * Q_q**2)/(6 * s)
        A2 = (G_F**2 * m_Z0**4 * s * (cal**2 + cvl**2) * (caq**2 + cvq**2))/(768 * pi * ((s - m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A3 = -(alpha * Q_q * G_F * m_Z0**2 * (s-m_Z0**2) * cvl * cvq)/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        A = A1 + A2 + A3

        B1 = (8 * G_F**2 * m_Z0**4 * s * cal * caq * cvl * cvq)/(768 * pi * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B2 = -(2 * alpha * Q_q * G_F * m_Z0**2 * cal * caq * (s - m_Z0**2))/(24 * 2**0.5 * ((s-m_Z0**2)**2 + m_Z0**2 * g_z**2))
        B = B1 + B2

        xsec2 = A * (1 + cost**2) + B * cost 
        xsec = (1./3.) * xsec1 + (2./3.) * xsec2

        return xsec


'''
generate an up type LQ with mass 2 TeV for mll = 1300 GeV with ylq = 1
'''

def SM_cost_fit(cost,A,B):
        return (A * (1 + cost**2) + B * cost)

mll = 1100.
s = mll * mll
m_lq_list = [1000., 1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.]


#y_lq_list = [0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,.85,0.9,1.]
#,0.8,0.85,0.9,0.95,1.]
#y_lq_list = [0.]
cost_list = np.linspace(-1,1,10000)



SM_xsec = SM_cost(0.3, s)
SM_norm = quad(lambda x: SM_cost(x,s), -1., 1.)[0]
print(" SM xsec at cost 0.3 = ", SM_xsec, "SM_norm = ",SM_norm)
y_lq_list = np.linspace(0.0,2.,100)
chan = 'mumu'
        
for flag in [1]:
        results = np.array([0.,0.,0.,0.])
        for m_lq in m_lq_list:
                
                # if m_lq < 2500:
                #         if chan=='ee': y_lq_list = np.linspace(0.05,0.7,40)
                #         else: y_lq_list = np.linspace(0.1,0.8,40)
                # elif m_lq >= 2500 and m_lq < 4000:
                #         if chan=='ee': y_lq_list = np.linspace(0.3,1.2,40)
                #         else: y_lq_list = np.linspace(0.5,1.3,40)
                # elif m_lq >= 4000 and m_lq < 6000:
                #         if chan=='ee': y_lq_list = np.linspace(0.5,1.6,40)
                #         else: y_lq_list = np.linspace(0.8,2.,40)
                for y_lq in y_lq_list:
                        #print("LQ mass = ", m_lq," y_lq = ",y_lq)
                        LQ_norm = quad(lambda x: LQ_cost(flag, x, s, m_lq, y_lq), -1., 1.)[0]
                        LQ_xsec = []
                        for cost in cost_list:
                                
                                LQ_xsec.append(LQ_cost(flag, cost, s, m_lq, y_lq)/LQ_norm)

                        AB,pcov = curve_fit(SM_cost_fit, cost_list, LQ_xsec)
                        A,B = AB[0], AB[1]
                        #print(A,B)
                        # plt.plot(cost_list,LQ_xsec,'b-',label='LQ signal - m_LQ = %.1f TeV, y_LQ = %.3f'%(m_lq/1000., y_lq))
                        # plt.plot(cost_list, SM_cost_fit(cost_list, A, B), 'r-', label='fit: A = %f, B = %f'%(A,B))
                        # plt.xlabel(r'$cos\theta$')
                        # plt.ylabel(r'$(1/\sigma)d\sigma/dcos\theta$')
                        # plt.title(r'Fitting $S_{\ell u}$ signal xsec with SM DY only terms')
                        # plt.legend()
                        # plt.savefig('SMP_reinterpret_fits/m%i_yLQ%.3f_SMtoLQfit.png'%(m_lq,y_lq))
                        # plt.close()

                        results = np.vstack((results,[m_lq, y_lq, A, B]))

        np.savetxt("smp_reinterpret_%s.txt"%chan,results)

        mlq_l, ylq_l, B_l = results[:,0], results[:,1], results[:,3]
        afb_limit=[]
        for mlq in m_lq_list:
                print("chan = ",chan," flag = ", flag, " mLQ = ", mlq)
                ylq = ylq_l[mlq_l==mlq]
                B = B_l[mlq_l==mlq]
                #print(ylq,B)
                
                # plt.axline((ylq[0],0.5318),(ylq[-1],0.5318),'b-',label=r'$\pm 2\sigma$')
                # plt.axline((ylq[0],0.5849),(ylq[-1],0.55849),'g-',label=r'$\pm 1\sigma$')
                # plt.axline((ylq[0],0.638),(ylq[-1],0.638),'r-',label=r'$A_{FB}$ in mass bin > 1 TeV from SMP-21-002')
                # plt.axline((ylq[0],0.6911),(ylq[-1],0.6911),'g-')
                # plt.axline((ylq[0],0.7442),(ylq[-1],0.7442),'b-')
                # plt.axhline(y=0.5318,color='b',label=r'$\pm 2\sigma$')
                # plt.axhline(y=0.5849,color='g',label=r'$\pm 1\sigma$')
                if chan=="ee":
                        plt.axhline(y=0.6091,color='r',label=r'$A_{FB}$ when $y_{e%s} = 0$'%("u" if flag!=1 else "d"))
                        # plt.axhline(y=0.6911,color='g')
                        # plt.axhline(y=0.7442,color='b')
                        plt.fill_between(np.linspace(-0.5,2.),0.4547,0.5319,color='yellow',label=r'$\pm 2\sigma$ (from SMP-21-002)')
                        plt.fill_between(np.linspace(-0.5,2.),0.6863,0.7635,color='yellow')
                        plt.fill_between(np.linspace(-0.5,2.),0.6091,0.6863,color='chartreuse',label=r'$\pm 1\sigma$ (from SMP-21-002)')
                        plt.fill_between(np.linspace(-0.5,2.),0.5319,0.6091,color='chartreuse')
                        plt.plot(ylq,B,color='black',label=r'$A_{FB}$ extracted from analytical fit to LQ signal')
                        s = [0.4547]#,0.694,0.77,0.846]
                        for sig in s:
                                diff, ylq_diff = B[B - sig < 0.0001], ylq[B - sig < 0.0001]
                                diff, ylq_diff = diff[sig - diff > 0.0001], ylq_diff[sig - diff > 0.0001]
                                afb_limit.append(ylq_diff[0])

                        plt.plot(ylq_diff[0],diff[0],'go') 
                        plt.xlabel(r"$y_{e%s}$"%("u" if flag!=1 else "d"))
                        plt.ylabel(r"$A_{FB}$")
                        plt.title(r'$S_{e%s}$ mass = %.1f TeV, $m_{ee}$ = %.2f TeV'%("u" if flag==2 else "d",mlq/1000.,mll/1000))
                        plt.ylim(0,1.2)
                        plt.xlim(np.amin(ylq),np.amax(ylq))
                        plt.legend()
                        plt.savefig('SMP_reinterpret_fits/m%i_yLQAFB_e%s.png'%(mlq,"u" if flag==2 else "d"))
                        plt.close()

                        
                else:
                        plt.axhline(y=0.6091,color='r',label=r'$A_{FB}$ when $y_{\mu %s} = 0$'%("u" if flag!=1 else "d"))
                        # plt.axhline(y=0.6911,color='g')
                        # plt.axhline(y=0.7442,color='b')
                        plt.fill_between(np.linspace(-0.5,2.),0.4681,0.5386,color='yellow',label=r'$\pm 2\sigma$ (from SMP-21-002)')
                        plt.fill_between(np.linspace(-0.5,2.),0.6796,0.75,color='yellow')
                        plt.fill_between(np.linspace(-0.5,2.),0.6091,0.6796,color='chartreuse',label=r'$\pm 1\sigma$ (from SMP-21-002)')
                        plt.fill_between(np.linspace(-0.5,2.),0.5386,0.6091,color='chartreuse')
                        plt.plot(ylq,B,color='black',label=r'$A_{FB}$ extracted from analytical fit to LQ signal')
                        s = [0.4681]#,0.595,0.665,0.735]
                        for sig in s:
                                diff, ylq_diff = B[B - sig < 0.0001], ylq[B - sig < 0.0001]
                                diff, ylq_diff = diff[sig - diff > 0.0001], ylq_diff[sig - diff > 0.0001]
                                afb_limit.append(ylq_diff[0])
                        plt.plot(ylq_diff[0],diff[0],'go') 
                        plt.xlabel(r"$y_{\mu %s}$"%("u" if flag!=1 else "d"))
                        plt.ylabel(r"A_{FB}")
                        plt.title(r'$S_{\mu %s}$ mass = %.1f TeV, $m_{\mu\mu}$ = %.2f TeV'%("u" if flag==2 else "d",mlq/1000.,mll/1000))
                        plt.ylim(0,1.2)
                        plt.xlim(np.amin(ylq),np.amax(ylq))
                        plt.legend()
                        plt.savefig('SMP_reinterpret_fits/m%i_yLQAFB_m%s.png'%(mlq,"u" if flag==2 else "d"))
                        plt.close()

        exp_limits = []
        with open("analyze/combine/AFB_fits/LQ_cards/%s/limit_json/limits_%s_091823.json"%(("d" if flag==1 else "u")+chan[0],("d" if flag==1 else "u")+chan[0]), 'r+') as f:
                data = json.load(f)
                for m in m_lq_list:
                        exp_limits.append(data[str(m)]['exp0'])

        afb_limit = np.array(afb_limit).reshape((9,1))
        plt.plot(m_lq_list, afb_limit, label='Limits recast from SMP-21-002')
        plt.plot(m_lq_list, exp_limits, label='Limits from EXO-22-013')
        #plt.xlim(0.,1.7)
        plt.ylim(0.,1.7)
        plt.xlabel(r"$M_{S_{%s%s}}$"%(chan[0],"u" if flag==2 else "d"))
        plt.ylabel(r"$y_{S_{%s%s}}$"%(chan[0],"u" if flag==2 else "d"))
        plt.title(r"$S_{%s%s}$ Limits"%(chan[0],"u" if flag==2 else "d"))
        plt.legend()
        plt.savefig("SMP_reinterpret_fits/lim_S%s%s.png"%(chan[0],"u" if flag==2 else "d"))
        plt.close()