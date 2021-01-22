from utils import *
import ROOT
from ROOT import *

par_name = 'METJER18'
fname = 'multidimfit.root'




f= ROOT.TFile.Open(fname)
fit = f.Get('fit_mdf')
#cov = fit.covarianceMatrix()
cov = fit.correlationMatrix()
cor = fit.correlation(par_name)


n_entries = cov.GetNrows()
par_list = fit.floatParsFinal()
par_idx = par_list.index(par_name)
par = par_list.find(par_name)
max_entry = 0
max_idx = -1

for i in range(n_entries):
    if(i == par_idx): continue
    entry = cov(par_idx, i)
    print("%s %.3e " % (par_list[i].GetName(), entry))
    if(entry > max_entry):
        max_entry = entry
        max_idx = i

print("Par %s : %.3e +/- %.3e" % (par_name, par.getVal(), par.getError()))
print("Max Cov %.3e, idx %i, name %s" % (max_entry, max_idx, par_list[max_idx].GetName()))
cor[max_idx].Print();
