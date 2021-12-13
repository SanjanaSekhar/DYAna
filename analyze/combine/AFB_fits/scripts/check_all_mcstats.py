import ROOT
from ROOT import *
from utils import *
import numpy as np


parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mbin",  default=1, type='int', help="Which mass bin to run on ")
(options, args) = parser.parse_args()

sigma2 = 0.1
if(options.mbin >=5): 
    sigma2 = 0.6




print_and_do("python scripts/check_mcstat_uncs.py --mbin %i --sigma2 %.2f --expected  "  % (options.mbin, sigma2))
print_and_do("python scripts/check_mcstat_uncs.py --mbin %i --sigma2 %.2f --expected  --fullCorr" % (options.mbin, sigma2))
print_and_do("python scripts/check_mcstat_uncs.py --mbin %i --sigma2 %.2f --expected  --noSymMCStats" % (options.mbin, sigma2))

