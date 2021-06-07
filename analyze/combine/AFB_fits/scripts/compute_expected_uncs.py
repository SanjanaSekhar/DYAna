import ROOT
from ROOT import *
from utils import *
import numpy as np

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject")
parser.add_option("--mbin",  default=-1, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="expected_fit_resuts/", help = "output directory")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")

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

if(options.mbin >= 0):
    print("Will only do fit for bin %i " % options.mbin)
    bin_start = options.mbin
    bin_stop = bin_start + 1


print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
print_and_do("""echo ".q" >> cmd.txt """)

seed = 12345

for mbin in range(bin_start, bin_stop):


    workspace = "workspaces/%s_fit_expected_unc_%i.root" % (options.chan, mbin)
    make_workspace(workspace, mbin,)



    print("Will inject AFB %.2f A0 %.2f for all toys " %(options.Afb, options.A0))

    print("Sampling toys based on postfit")
    if(not options.reuse_fit):
        print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s" % (workspace, extra_params))


    fitted_afb, fitted_a0 = setSnapshot(Afb_val = -1., mdf = True)
    #full fit
    print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -n FullFit -s %i --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters Afb=%.2f,A0=%.2f" 
        % (seed, options.Afb, options.A0))
    print_and_do(("combine -M MultiDimFit -d %s -n FullFit --saveWorkspace --saveFitResult --toysFile higgsCombineFullFit.GenerateOnly.mH120.%i.root " +
    "--toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s") %(workspace, seed, extra_params))
    #no systematics fit
    #print_and_do(("combine -M MultiDimFit  -d initialFitWorkspace.root -n NoSys -s %i --snapshotName initialFit  --saveFitResult " + 
    #"--toysNoSystematics --freezeParameter allConstrainedNuisances --bypassFrequentistFit -t 1  --setParameters Afb=%.2f,A0=%.2f") % (seed, options.Afb, options.A0))
    print_and_do(("combine -M MultiDimFit -d %s -n NoSys --saveWorkspace --saveFitResult --toysFile higgsCombineFullFit.GenerateOnly.mH120.%i.root " +
    "--toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL --freezeParameter allConstrainedNuisances %s") %(workspace, seed, extra_params))

    #No MC stat uncs fit
    #print_and_do(("combine -M MultiDimFit -d %s -n NoMCStat --saveWorkspace --saveFitResult --toysFile higgsCombineFullFit.GenerateOnly.mH120.%i.root " +
    #"--toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL --freezeNuisanceGroups autoMCStats %s") %(workspace, seed, extra_params))

    f_fit_nosys = TFile.Open("multidimfitNoSys.root")
    f_fit_full = TFile.Open("multidimfitFullFit.root")


    if f_fit_nosys:
        fr = f_fit_nosys.Get('fit_mdf')
        myargs = RooArgSet(fr.floatParsFinal())

        AFB_err_stat[mbin] = myargs.find("Afb").getError()
        A0_err_stat[mbin] = myargs.find("A0").getError()
        f_fit_nosys.Close()
        print_and_do("root -l -b multidimfitNoSys.root < cmd.txt > expected_fit_results/%s_fit_results_mbin%i.txt" % (fit_name + "_nosys", mbin))
    else:
        print("Can't open multidimfit. Looks like the fit failed.")

    if f_fit_full:
        fr = f_fit_full.Get('fit_mdf')
        myargs = RooArgSet(fr.floatParsFinal())

        AFB_err_full[mbin] = myargs.find("Afb").getError()
        A0_err_full[mbin] = myargs.find("A0").getError()
        f_fit_full.Close()
        print_and_do("root -l -b multidimfitFullFit.root < cmd.txt > expected_fit_results/%s_fit_results_mbin%i.txt" % (fit_name, mbin))
    else:
        print("Can't open multidimfit. Looks like the fit failed.")

    print_and_do("rm higgsCombine* multidim* initialFitWorkspace.root")


n_bins = 8
for i in range(1,n_bins):

    AFB_err_sys[i] = (AFB_err_full[i]**2 - AFB_err_stat[i]**2 + AFB_shifts_unc[i]**2)**0.5
    A0_err_sys[i] = (A0_err_full[i]**2 - A0_err_stat[i]**2)**0.5

print("AFB:")
for i in range(1,n_bins):
    #print("%.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(AFB_val[i], AFB_err_stat[i], AFB_err_sys[i]))
    print("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(AFB_err_stat[i], AFB_err_sys[i]))

print("A0:")
for i in range(1,n_bins):
    #print("%.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(A0_val[i], A0_err_stat[i], A0_err_sys[i]))
    print("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(A0_err_stat[i], A0_err_sys[i]))
