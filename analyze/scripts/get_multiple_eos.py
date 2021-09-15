

import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)


#labels = ["gof_postfit_y16", "gof_postfit_y17","gof_postfit_y18", "gof_postfit", "gof_postfit_y16_mumu", "gof_postfit_y16_ee" ]

labels = [ 
        "gof_prefit",
        "sys_uncs",
        "impacts",  
        "sys_uncs_diff",
        "impacts_diff",  
        #'gof_prefit_mumu',
        #'gof_prefit_ee',
        #'gof_noSym',
        #'gof_mumu16',
        #'gof_mumu17',
        #'gof_mumu18',
        #'gof_ee16',
        #'gof_ee17',
        #'gof_ee18',

#"gof_prefit",
        #"gof_postfit", 
        #"impacts",
        #"expected_impacts",
        #"expected_sys_uncs",
        #"sys_uncs",
        #"bias_test_afb6_a00", "bias_test_afb6_a01", "bias_test_afb0_a00", "bias_test_afb0_a01",
#"gof_postfit_y16", "gof_postfit_y17","gof_postfit_y18",
#"gof_postfit_y16_ee", "gof_postfit_y16_mumu",
]
date = "Sep01"

odir = "../plots/Misc_plots/fit_checks_sep01_unblind/"
#odir = "../plots/Misc_plots/"
print_and_do("mkdir %s" % odir)



for i,label in enumerate(labels):
    out = odir + label + "_" + date
    print_and_do("mkdir %s" % out)
    print_and_do("python doCondor.py -g -o %s -n %s" %(out, label + "_" + date))



