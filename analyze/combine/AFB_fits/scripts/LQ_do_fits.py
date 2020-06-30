from LQ_utils import *
import ROOT
from ROOT import *

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--no_plot",  default=False, action="store_true", help="Don't make postfit plots")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--fake_data",  default=False, action="store_true", help="Use fit template without any shape systematics and no fakes")
parser.add_option("--no_cleanup",  default=False, action="store_true", help="Don't remove root files created by fit")
parser.add_option("--mbin", default = -1, type='int', help="Only do fits for this single mass bin, default is all bins")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")

(options, args) = parser.parse_args()


extra_params = ""

if(options.chan == "ee"):
    print("Chan is ee, will mask mumu channels")
    if(options.year < 0): 
        extra_params += " --setParameters mask_Y16_mumu16=1,mask_Y17_mumu17=1,mask_Y18_mumu18=1" 
    else:
        extra_params += " --setParameters mask_Y%i_mumu%i=1" % (options.year % 2000, options.year % 2000)
elif(options.chan == "mumu"):
    print("Chan is mumu, will mask ee and ee_ss channels")
    if(options.year < 0):
        extra_params += " --setParameters mask_Y16_ee16=1,mask_Y17_ee17=1,mask_Y18_ee18=1,mask_Y16_ee16_ss=1,mask_Y17_ee17_ss=1,mask_Y18_ee18_ss=1" 
    else:
        extra_params += " --setParameters mask_Y%i_ee%i=1,mask_Y%i_ee%i_ss=1" % (options.year % 2000, options.year % 2000, options.year % 2000, options.year % 2000)




bin_start = 0
bin_stop = 8
options.chan="ee"
options.no_sys=True
fit_name = options.chan
if(options.no_sys): fit_name +="_nosys"
if(options.fake_data): fit_name +="_fake_data"
if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)
print("\n fit_name = ", fit_name)

if(options.mbin >= 0):
    print("Will only do fit for bin %i " % options.mbin)
    bin_start = options.mbin
    bin_stop = bin_start + 1

for mLQ in [1000,3000]:
#mLQ = 1000.
#for mbin in range(bin_start, bin_stop):
#print(" \n \n Starting fit for bin %i \n\n" % mbin)
    print(" \n \n Starting fit for LQ m = %i\n\n",mLQ)

    workspace="workspaces/%s_LQ.root" % (options.chan)
    make_workspace(workspace, options.no_sys, False, mLQ, year = options.year)
    plotdir="postfit_plots/%s_LQ_m%i" % (fit_name,mLQ)
    print("\n plotdir = ", plotdir)
    print_and_do("[ -e %s ] && rm -r %s" % (plotdir, plotdir))
    print_and_do("mkdir %s" % (plotdir))
    print_and_do("combine %s -M MultiDimFit  --saveWorkspace --saveFitResult --robustFit 1 %s" %(workspace, extra_params))

    if(not options.no_plot):
        print_and_do("PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfit.root:fit_mdf --postfit -o %s_fit_shapes_LQ.root --sampling --samples 100"
                % (fit_name))
        extra_args = ""
        if(options.year > 0): extra_args = " -y %i " % options.year
        print_and_do("python scripts/LQ_plot_postfit.py -i %s_fit_shapes_LQ.root -o %s  %s --mLQ %i" % (fit_name, plotdir, extra_args,mLQ))
        print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s" % (workspace, extra_params)) #only to get prefit, probably a better way
        print_and_do("python scripts/my_diffNuisances.py multidimfit.root --multidim --mLQ %i --prefit fitDiagnostics.root -p Afb --skipFitB -g %s" % (mLQ, plotdir))
        print_and_do("mv %s_fit_shapes_LQ.root %s" %(fit_name, plotdir))
        if(not options.no_cleanup): print_and_do("rm fitDiagnostics.root higgsCombineTest.FitDiagnostics.mH120.root")


    print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
    print_and_do("""echo ".q" >> cmd.txt """)
    print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_fake_nosys_LQ_m%i.txt" % (fit_name,mLQ))
    print_and_do("rm -f cards/sed*")
    if(not options.no_cleanup): print_and_do("rm cmd.txt combine_logger.out higgsCombineTest.MultiDimFit.mH120.root multidimfit.root")
