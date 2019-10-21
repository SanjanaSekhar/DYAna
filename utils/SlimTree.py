import ROOT as r

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product
#from BaconAna.Utils.makeFilelist import *

default_args = []

# Options
parser = OptionParser()
parser.add_option("-c", "--cut", dest='cut_str', default = "m>130.", help="String of cut to apply to tree")
(options, args) = parser.parse_args()

fin_str = sys.argv[1]
fout_str = sys.argv[2]

print("Copying trees from %s to %s. Using cut string %s" % (fin_str, fout_str, options.cut_str)

fin = r.TFile.Open(fin_str, "READ")
fout = r.TFile.Open(fout_str, "RECREATE")

trees = fin.GetListOfKeys()

for t_key in trees:
    t = fin.Get(t_key.GetName())
    print(t)
    t_new = t.CopyTree(options.cut_str)
    fout.cd()
    t_new.Write()
    

