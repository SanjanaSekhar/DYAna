from ROOT import *
import os
import sys
from optparse import OptionParser
from optparse import OptionGroup

def get_avg_err(h):
    nbins = h.GetNbinsX()
    avg = 0.
    total = h.Integral()

    for i in range(1,nbins+1):
        err = h.GetBinError(i)
        cont = h.GetBinContent(i)

        #same as weighted average of fractional error 
        avg += err

    avg /= total
    return avg



fin = sys.argv[1]

years = [16,17,18]


dir_names = [ 'Y%i_ee%i_prefit', 'Y%i_mumu%i_prefit']
template_names = ['TotalSig', 'db', 'gam', 'qcd', 'top']

f = TFile.Open(fin, "READ")

for year in years:
    for w_dir  in dir_names:
        n_dir = w_dir % (year, year)
        f.cd(n_dir)
        print("\n\nDir %s:" % n_dir)
        for name  in template_names:
            h = gDirectory.Get(name)
            avg_err = get_avg_err(h)
            print("%s : %.3f " % (name, avg_err))


