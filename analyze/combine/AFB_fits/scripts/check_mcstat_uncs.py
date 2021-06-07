import ROOT
from ROOT import *
from utils import *
import numpy as np
from add_group_impact import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject")
parser.add_option("--mbin",  default=1, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="expected_fit_resuts/", help = "output directory")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("--expected",  default=False, action="store_true", help="Compute expected impacts based on toys with AFB=0.6 A0=0.05")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--noSymMCStats", default = False, action="store_true",  help="Don't add constraints to mcStat nuisances")

(options, args) = parser.parse_args()

AFB_shifts = [0.0, 0.016, 0.009, 0.005, 0.002, 0.002, -0.001, -0.001]
AFB_shifts_unc = [0.0, 0.003, 0.002, 0.002, 0.001, 0.002, -0.001, -0.001]

gStyle.SetOptFit(1) 

extra_params = ""

if(options.chan == "ee"):
    print("Chan is ee, will mask mumu channels")
    if(options.year < 0): 
        extra_params += " --setParameters mask_Y16_mumu16=1,mask_Y17_mumu17=1,mask_Y18_mumu18=1" 
    else:
        extra_params += " --setParameters mask_Y%i_mumu%i=1" % (options.year % 2000, options.year % 2000)
elif(options.chan == "mumu"):
    print("Chan is mumu, will mask ee channels")
    if(options.year < 0):
        extra_params += " --setParameters mask_Y16_ee16=1,mask_Y17_ee17=1,mask_Y18_ee18=1" 
    else:
        extra_params += " --setParameters mask_Y%i_ee%i=1" % (options.year % 2000, options.year % 2000, options.year % 2000, options.year % 2000)

extra_params += " --X-rtd MINIMIZER_no_analytic"

bin_start = 1
bin_stop = 8

AFB_err_stat = [0.]*8
AFB_err_sys = [0.]*8
AFB_err_full = [0.]*8

A0_err_stat = [0.]*8
A0_err_sys = [0.]*8
A0_err_full = [0.]*8

fit_name = options.chan + "_expected"
if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)



print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
print_and_do("""echo ".q" >> cmd.txt """)

s = 123456



workspace = "workspaces/%s_fit_expected_unc_%i.root" % (options.chan, options.mbin)

if(not options.reuse_fit):
    make_workspace(workspace, options.mbin, year = options.year, symMCStats = not options.noSymMCStats)
    print_and_do("combine -M MultiDimFit -d %s --saveFitResult --saveWorkspace -n _base --robustFit 1  %s" % (workspace, extra_params))



if(options.expected):
    print("Will inject AFB %.2f A0 %.2f for all toys " %(options.Afb, options.A0))


    print_and_do(("combine -M GenerateOnly -d higgsCombine_base.MultiDimFit.mH120.root --snapshotName MultiDimFit --toysFrequentist"
    " --bypassFrequentistFit --saveToys -t 1 -s %i  --setParameters Afb=%.2f,A0=%.2f")
            % (s, options.Afb, options.A0))

    print_and_do(("combine -M MultiDimFit -d %s -n _nom --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root " +
    "--toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s") %(workspace, s, extra_params))

    extra_params += " --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist -t 1" % s

else:
    print_and_do("cp higgsCombine_base.MultiDimFit.mH120.root higgsCombine_nom.MultiDimFit.mH120.%i.root" % (s))

print_and_do("""combine -M FitDiagnostics -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _nom1 %s""" % (s,extra_params))

print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups autoMCStats -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _NoMCStat %s """ 
        % (s, extra_params))

if(not options.expected): s = -1
sys_unc = compute_sys("nom1", "NoMCStat", s)

print("Mass bin %i, Symmetrized MCStats = %r, expected = %r \n" %(options.mbin, not options.noSymMCStats, options.expected))
print("Uncertainty is %.5f \n" % sys_unc)

#print_and_do("rm higgsCombine* fitDiagnostics* multidimfit*")

