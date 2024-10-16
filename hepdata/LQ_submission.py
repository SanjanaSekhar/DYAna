from hepdata_lib import *
from hepdata_lib.helpers import *
import numpy as np

sub = Submission()

sub.read_abstract("hepdata_inputs/abstract.txt")
sub.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/EXO-22-013/")

vec = False
i = 0
for chan in ['e','mu']:
    for q in ['u','d']:

        pf.append(Table("Figure %i"%(i+3)))
        pf[i].description = "The observed data in the dielectron channel and the fitted signal-plus-background templates, shown for the $%s_{%s%s}$ scenario with a candidate LQ mass of 2.5 TeV. The black points are the observed data, the stacked histograms represent the backgrounds, and the yellow histogram shows the fitted LQ signal multiplied by 10. Distributions of events are binned in the reconstructed dilepton mass (vertical red dashed lines), rapidity (vertical black dashed lines), and cosine theta" %(("V" if vec else "S"), ('\mu' if chan=='m' else chan),q)
        pf[i].location = "Data from Figure %i"%(i+3)
        pf[i].keywords["observables"] = ["D3N/DM/DYRAP/DTHETA"]
        pf[i].keywords["cmenergies"] = [13000]
        pf[i].keywords["reactions"] = ["P P --> LQ* --> LEPTON+ LEPTON-","%sQ %sQBAR --> LQ* --> %s+ %s-"%(q.upper(),chan.upper())]
        data = np.loadtxt("hepdata_inputs/%s_%s%s_postfit_table.txt"%(chan+chan,q,("_vec" if vec else "")), skiprows=1)
        # df[['Bin index','Background','Background err','LQ template yield','LQ template yield err','Observed data','Observed data err']]
        bin = Variable("Bin index", is_independent=True, is_binned=True, units="")
        bin.values = data[:,0]
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
        pf[i].add_image("hepdata_inputs/%s_%s%s_unblinded_LQ_m2500/Postfit_%s%s_m2500.pdf"%(chan+chan,q,("_vec" if vec else ""),q,chan[0]))
        sub.add_table(pf[i])

        i+=1

outdir="hepdata_outputs/"
sub.create_files(outdir)