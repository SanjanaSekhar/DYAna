from hepdata_lib import *
from hepdata_lib.helpers import *
import numpy as np
import json
import pandas as pd

sub = Submission()

sub.read_abstract("hepdata_inputs/abstract.txt")
sub.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/EXO-22-013/")

pf = []

i, j = 0, 0
for vec in [False,True]:
    for chan in ['e','mu']:
        for q in ['u','d']:

            pf.append(Table("Figure %i"%(i+3)))
            pf[i].description = "The observed data in the %s channel and the fitted signal-plus-background templates, shown for the $%s_{%s %s}$ scenario with a candidate LQ mass of 2.5 TeV. Distributions of events are binned in the reconstructed dilepton mass, rapidity, and cosine theta." %(("dielectron" if chan=='e' else "dimuon"), ("V" if vec else "S"), ('\mu' if chan=='mu' else chan),q)
            pf[i].location = "Data from Figure %i"%(i+3)
            pf[i].keywords["observables"] = ["D3N/DM/DYRAP/DTHETA"]
            pf[i].keywords["cmenergies"] = [13000]
            pf[i].keywords["reactions"] = ["P P --> LQ* --> LEPTON+ LEPTON-","%sQ %sQBAR --> LQ* --> %s+ %s-"%(q.upper(),q.upper(),chan.upper(),chan.upper())]
            data = np.loadtxt("hepdata_inputs/%s_%s%s_postfit_table.txt"%(chan+chan,q,("_vec" if vec else "")), skiprows=1)
            # df[['Bin index','Background','Background err','LQ template yield','LQ template yield err','Observed data','Observed data err']]
            bin = Variable("Bin index", is_independent=True, is_binned=True, units="")
            bin.values = [(i,i+1) for i in range(1,61)]
            bkg = Variable("Total Background", is_independent=False, is_binned=False, units="")
            bkg.values = data[:,1]
            bkg_unc = Uncertainty("", is_symmetric=True)
            bkg_unc.values = data[:,2]
            bkg.add_uncertainty(bkg_unc)
            lq = Variable("LQ template yield", is_independent=False, is_binned=False, units="")
            lq.values = data[:,3]
            lq_unc = Uncertainty("", is_symmetric=True)
            lq_unc.values = data[:,4]
            lq.add_uncertainty(lq_unc)
            d = Variable("Observed data", is_independent=False, is_binned=False, units="")
            d.values = data[:,5]
            d_unc = Uncertainty("", is_symmetric=True)
            d_unc.values = data[:,6]
            d.add_uncertainty(d_unc)

            pf[i].add_variable(bin)
            pf[i].add_variable(bkg)
            pf[i].add_variable(lq)
            pf[i].add_variable(d)
            # pf[i].add_image("hepdata_inputs/%s_%s%s_unblinded_LQ_m2500/Postfit_%s%s_m2500.pdf"%(chan+chan,q,("_vec" if vec else ""),q,chan[0]))
            sub.add_table(pf[i])

            with open("hepdata_inputs/limits_%s%s%s_031124_unblinding_y-2001.json"%(q, chan[0], ("_vec" if vec else "")), 'r+') as f:
                data = json.load(f)
            lim_df = pd.DataFrame(data).T
            lim_df.sort_index(ascending=True, inplace=True)

            lim = Table("Figure %i (%s)"%((j+11),("left" if i%2==0 else "right")))
            lim.description = "Upper limits at $95\%%$ CL on the LQ-fermion couplings, $|%s_{%s %s}|$, versus $m_{LQ}$ for %s LQs coupled to %s." %(("g" if vec else "y"), ("e" if chan=='e' else '\mu'), q, ("vector" if vec else "scalar"), ("electrons" if chan=='e' else "muons")) 
            lim.location = "Data from Figure %i (%s)"%((j+11),("left" if i%2==0 else "right"))
            lim.keywords["observables"] = ["CLS"]
            lim.keywords["cmenergies"] = [13000]
            lim.keywords["reactions"] = ["P P --> LQ* --> LEPTON+ LEPTON-","%sQ %sQBAR --> LQ* --> %s+ %s-"%(q.upper(),q.upper(),chan.upper(),chan.upper())]

            bin = Variable("$%s_{%s %s}$ Mass" %(("V" if vec else "S"), ('\mu' if chan=='mu' else chan),q), is_independent=True, is_binned=False, units="GeV")
            bin.values = [float(a) for a in lim_df.index.to_list()]

            exp = Variable("Expected limit on $|%s_{%s %s}|$" %(("g" if vec else "y"), ("e" if chan=='e' else '\mu'), q), is_independent=False, is_binned=False, units="")
            exp.values = lim_df['exp0'].to_list()
            exp_unc = Uncertainty("1 s.d.", is_symmetric=False)
            exp_unc.values = [(a,b) for a,b in zip((lim_df['exp-1']-lim_df['exp0']).to_list(),lim_df['exp+1']-(lim_df['exp0']).to_list())]
            exp.add_uncertainty(exp_unc)
            exp_unc2 = Uncertainty("2 s.d.", is_symmetric=False)
            exp_unc2.values = [(a,b) for a,b in zip((lim_df['exp-2']-lim_df['exp0']).to_list(),lim_df['exp+2']-(lim_df['exp0']).to_list())]
            exp.add_uncertainty(exp_unc2)

            obs = Variable("Observed limit on $|%s_{%s %s}|$" %(("g" if vec else "y"), ("e" if chan=='e' else '\mu'), q), is_independent=False, is_binned=False, units="")
            obs.values = lim_df['obs'].to_list()
            

            lim.add_variable(bin)
            lim.add_variable(exp)
            lim.add_variable(obs)
            
            # pf[i].add_image("hepdata_inputs/%s_%s%s_unblinded_LQ_m2500/Postfit_%s%s_m2500.pdf"%(chan+chan,q,("_vec" if vec else ""),q,chan[0]))
            sub.add_table(lim)

            if i%2 != 0: j+=1
            i+=1

acc = Table("Acceptance x Efficiency fractions for the dielectron channel")
acc.description = "Acceptance x efficiency for the dielectron channel for 2016, 2017 and 2018. Electrons are required to have opposite charge, and one of two electrons in the event must pass the high level trigger with $p_T > 27/35/32$ GeV for each year. We apply an offline selection of $p_T > 40$ GeV for the leading electron, and $p_T > 15$ GeV for the subleading electron. We require that both electrons have $|\eta| < 1.44$ or $1.56 <|\eta| < 2.5$ and $m_{\ell\ell} > 500$ GeV. Additionally, the electrons are subject to the tight working points of cut-based identification and isolation. Distributions of events are binned in the generator level dilepton mass, rapidity, and cosine theta."

acc.keywords["observables"] = ["ACC, EFF"]
acc.keywords["cmenergies"] = [13000]
acc.keywords["reactions"] = ["P P --> LQ* --> E+ E-"]
           
data = np.loadtxt("hepdata_inputs/elel_acc_eff_y16.txt")
bin = Variable("Bin index", is_independent=True, is_binned=True, units="")
bin.values = [(i,i+1) for i in range(1,61)]
acc_16 = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2016)", is_independent=False, is_binned=False, units="")
acc_16.values = data[:,4]
data = np.loadtxt("hepdata_inputs/elel_acc_eff_y17.txt")
acc_17 = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2017)", is_independent=False, is_binned=False, units="")
acc_17.values = data[:,4]
data = np.loadtxt("hepdata_inputs/elel_acc_eff_y18.txt")
acc_18 = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2018)", is_independent=False, is_binned=False, units="")
acc_18.values = data[:,4]

acc.add_variable(bin)
acc.add_variable(acc_16)
acc.add_variable(acc_17)
acc.add_variable(acc_18)

sub.add_table(acc)    

acc_mu = Table("Acceptance x Efficiency fractions for the dimuon channel")
acc_mu.description = "Acceptance x efficiency for the dimuon channel for 2016, 2017 and 2018. Muon are required to have opposite charge, and one of two muons in the event must pass the high level trigger with $p_T > 24/27/24$ GeV for each year. We apply an offline selection of $p_T > 40$ GeV for the leading muon, and $p_T > 15$ GeV for the subleading muon. We require that both muon have $|\eta| < 2.4$ and $m_{\ell\ell} > 500$ GeV. Additionally, the muons are subject to the tight working points of cut-based identification and isolation. Distributions of events are binned in the generator level dilepton mass, rapidity, and cosine theta."

acc_mu.keywords["observables"] = ["ACC, EFF"]
acc_mu.keywords["cmenergies"] = [13000]
acc_mu.keywords["reactions"] = ["P P --> LQ* --> MU+ MU-"]
           
data = np.loadtxt("hepdata_inputs/mumu_acc_eff_y16.txt")
bin_mu = Variable("Bin index", is_independent=True, is_binned=True, units="")
bin_mu.values = [(i,i+1) for i in range(1,61)]
acc_16_mu = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2016)", is_independent=False, is_binned=False, units="")
acc_16_mu.values = data[:,4]
data = np.loadtxt("hepdata_inputs/mumu_acc_eff_y17.txt")
acc_17_mu = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2017)", is_independent=False, is_binned=False, units="")
acc_17_mu.values = data[:,4]
data = np.loadtxt("hepdata_inputs/mumu_acc_eff_y18.txt")
acc_18_mu = Variable("Acc. x eff. fraction after $p_T$ and $\eta$ cuts (2018)", is_independent=False, is_binned=False, units="")
acc_18_mu.values = data[:,4]

acc_mu.add_variable(bin_mu)
acc_mu.add_variable(acc_16_mu)
acc_mu.add_variable(acc_17_mu)
acc_mu.add_variable(acc_18_mu)

sub.add_table(acc_mu)     


outdir="hepdata_outputs/"
sub.create_files(outdir)