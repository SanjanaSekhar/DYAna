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
            bin.values = [(i,i+1) for i in range(1,60)]
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

outdir="hepdata_outputs/"
sub.create_files(outdir)