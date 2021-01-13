

import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)


#labels = ["gof_postfit_y16", "gof_postfit_y17","gof_postfit_y18", "gof_postfit", "gof_postfit_y16_mumu", "gof_postfit_y16_ee" ]

labels = [ "impacts",
        "bias_test_afb6_a00", "bias_test_afb6_a01", "bias_test_afb0_a00", "bias_test_afb0_a01",
"gof_postfit", 
#"gof_postfit_y16", "gof_postfit_y17","gof_postfit_y18",
#"gof_postfit_y16_ee", "gof_postfit_y16_mumu",
]
date = "Jan12"

odir = "../plots/Misc_plots/fit_checks_jan12/"
print_and_do("mkdir %s" % odir)



for i,label in enumerate(labels):
    out = odir + label + "_" + date
    print_and_do("mkdir %s" % out)
    print_and_do("python doCondor.py -g -o %s -n %s" %(out, label + "_" + date))



