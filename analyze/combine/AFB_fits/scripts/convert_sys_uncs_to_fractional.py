import operator
import ROOT
from utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-i", "--input", default="", help = "Input file name")
parser.add_option("-o", "--output", default="", help = "output file name")
(options, args) = parser.parse_args()

if(options.output == ""):
    options.output = options.input[:-4] + "_fractional" + ".txt"

f_in = open(options.input, 'r')
sys_dict = dict()
tot_unc = 0.

header = f_in.readline()
mbin = int(header.split()[-1])

for l in f_in:
    p  = l.split('&')
    name= p[0]

    unc2 = (float(p[1].split()[0]) * 10**(-3))**2
    print(name, unc2)
    sys_dict[name] = unc2
    tot_unc +=unc2

print("total", tot_unc)

with open(options.output, 'w') as f_out:
    sorted_d = sorted(sys_dict.items(), key=operator.itemgetter(1))
    f_out.write("Fractional Systematic uncertainties for bin %i \n" % mbin)
    for sys_name, val in sorted_d[::-1]:
        frac = (val / tot_unc) * 100.

        if(frac > 0.1):
            f_out.write("%s & %.1f%% \\\\ \n" % (sys_name, frac))
        else:
            f_out.write("%s & $<$0.1\\%% \\\\ \n" % (sys_name))

