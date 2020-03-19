import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product


def print_and_do(s):
    print("Exec: " + s)
    os.system(s)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--no_plot",  default=False, action="store_true", help="Don't make postfit plots")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--no_cleanup",  default=False, action="store_true", help="Don't remove root files created by fit")

(options, args) = parser.parse_args()


extra_params = ""

if(options.chan == "ee"):
    print("Chan is ee, will mask mumu channels")
    extra_params += " --setParameters mask_Y16_mumu16=1,mask_Y17_mumu17=1,mask_Y18_mumu18=1" 
elif(options.chan == "mumu"):
    print("Chan is mumu, will mask ee and ee_ss channels")
    #extra_params += " --setParameters mask_Y16_ee16=1,mask_Y17_ee17=1,mask_Y18_ee18=1,mask_Y16_ee16_ss,mask_Y17_ee17_ss,mask_Y18_ee18_ss" 



cleanup = False


#for mbin in range(8):
for mbin in range(1):

    template_card="card_templates/combined_fit_template.txt"
    if(options.no_sys): template_card = "card_templates/combined_fit_template.txt"
    comb_card = "cards/combined_fit_mbin%i.txt" % (mbin)
    for year in (16,17,18):
        card="cards/combined_fit_y%i_mbin%i.txt" % (year, mbin)
        print_and_do("cp %s %s" % (template_card, card))
        print_and_do("""sed -i "s/YR/%i/g" %s""" % (year, card))

    print_and_do("combineCards.py Y16=cards/combined_fit_y16_mbin%i.txt Y17=cards/combined_fit_y17_mbin%i.txt Y18=cards/combined_fit_y18_mbin%i.txt > %s" % (mbin, mbin, mbin, comb_card))



    workspace="workspaces/%s_fit_%i.root" % (options.chan, mbin)
    plotdir="postfit_plots/%s_mbin%i" % (options.chan, mbin)
    print_and_do("[ -e %s ] && rm -r %s" % (plotdir, plotdir))
    print_and_do("mkdir %s" % (plotdir))
    print_and_do("text2workspace.py %s --keyword-value M_BIN=%i -P Analysis.DYAna.my_model:dy_AFB -o %s --channel-masks" % (comb_card, mbin, workspace))
    print_and_do("combine %s -M MultiDimFit  --saveWorkspace --saveFitResult --robustFit 1 %s" %(workspace, extra_params))

    if(not options.no_plot):
        print_and_do("PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfit.root:fit_mdf --skip-prefit --postfit -o %s_fit_shapes_mbin%i.root --sampling --samples 100"
                % (options.chan, mbin))
        print_and_do("python scripts/plot_postfit.py -i %s_fit_shapes_mbin%i.root -o %s -m %i" % (options.chan, mbin, plotdir, mbin))
        print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s" % (workspace, extra_params)) #only to get prefit, probably a better way
        print_and_do("python scripts/my_diffNuisances.py multidimfit.root --multidim --prefit fitDiagnostics.root -p Afb --skipFitB -g %s" % (plotdir))
        print_and_do("mv %s_fit_shapes_mbin%i.root %s" %(options.chan, mbin, plotdir))
        if(not options.no_cleanup): print_and_do("rm fitDiagnostics.root")


    print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
    print_and_do("""echo ".q" >> cmd.txt """)
    print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_fit_results_mbin%i.txt" % (options.chan, mbin))
    if(not options.no_cleanup): print_and_do("rm cmd.txt combine_logger.out higgsCombineTest.MultiDimFit.mH120.root multidimfit.root")

