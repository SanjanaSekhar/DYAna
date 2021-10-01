import ROOT
from ROOT import *
from contextlib import contextmanager
import os
from optparse import OptionParser
import CMS_lumi, tdrstyle

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("--input", "-i", default = "", help="Input dir")
    parser.add_option("--output", "-o", default = "", help="Input directory")
    parser.add_option("--mbin", "-m", type = 'int', default = -1, help="Mass bin (-1 for all)")
    parser.add_option("--prelim", default = False, action = "store_true", help = "Add 'Preliminary' label to plots")
    (options, args) = parser.parse_args()

if(options.output[-1] != "/"): options.output+="/"

if(options.mbin >0):
    mbins = [options.mbin]
else:
    mbins = range(1,8)

extra = ""
if(options.prelim): extra = " --prelim"

os.system("mkdir %s" % options.output)
for mbin in mbins:
    mbin_dirname = "combined_mbin%i/" % mbin
    root_filename = "combined_fit_shapes_mbin%i.root" % mbin

    output_dir = options.output + mbin_dirname
    os.system("mkdir %s" % output_dir)
    #os.system("python scripts/plot_postfit.py -i %s -o %s --mbin %i" % (options.input + mbin_dirname + root_filename, output_dir, mbin))
    os.system("python scripts/plot_comb_postfit.py -i %s -o %s -m %i %s" % (options.input + mbin_dirname + root_filename,  output_dir, mbin, extra))
    #os.system("python scripts/plot_swaped_axis_postfit.py -i %s -o %s --mbin %i" % (options.input + mbin_dirname + root_filename, output_dir, mbin))


