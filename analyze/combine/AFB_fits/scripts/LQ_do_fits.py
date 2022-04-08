
from LQ_utils import *
import ROOT
from ROOT import *

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--q",  default="combined", type="string", help="What channels to run the fit over (combined, u, or d)")
parser.add_option("--no_plot",  default=False, action="store_true", help="Don't make postfit plots")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--fake_data",  default=False, action="store_true", help="Use fit template without any shape systematics and no fakes")
parser.add_option("--no_cleanup",  default=False, action="store_true", help="Don't remove root files created by fit")
parser.add_option("--mbin", default = -1, type='int', help="Only do fits for this single mass bin, default is all bins")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("-v", "--verbose", default = 0, type='int', help="Turn up verbosity on fits")
parser.add_option("--noSymMCStats", default = False, action="store_true",  help="Don't add constraints to mcStat nuisances")
parser.add_option("--no_LQ",  default=False, action="store_true", help="For sanity check purposes remove LQ temps")
parser.add_option("--gen_level",  default=False, action="store_true", help="gen level fits")
(options, args) = parser.parse_args()



for y in [2017]:
    #for options.chan in ["mumu","ee"]:
    for options.chan in ["ee"]:
        for options.q in ["u"]:

            

            extra_params=""
#            options.chan="mumu"
#            options.q="u"
            options.no_sys=False
            options.fake_data=False
            options.no_LQ=False
            options.year = y
            '''
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
            '''
            if(options.verbose > 0):
        	extra_params +=" --verbose %i" % options.verbose

            #No analytic minimization of MC stats nuisances
            #extra_params += "--X-rtd MINIMIZER_no_analytic"

            
            fit_name = options.chan
            if(options.no_sys): 
                fit_name +="_nosys"
                #extra_params += " --freezeParameters allConstrainedNuisances"
            if(options.noSymMCStats):
                fit_name += "_noSymMC"
            if(options.fake_data): fit_name +="_fake_data"

            if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)
            fit_name+="_"+options.q

            if(options.no_LQ): fit_name+="_noLQ"

            if options.chan=="ee" and options.gen_level : fit_name+="_gen_level_SMdata"

            print("\n fit_name = ", fit_name)


            for mLQ in [1000]:
            #,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]:
            #mLQ = 1000.
            #for mbin in range(bin_start, bin_stop):
            #print(" \n \n Starting fit for bin %i \n\n" % mbin)
                print(" \n \n Starting fit for LQ m = %i\n\n",mLQ)

                workspace="workspaces/%s_LQ.root" % (options.chan)
                make_workspace(workspace, options.gen_level, options.chan, options.q, options.no_LQ, options.no_sys, options.fake_data, mLQ, year = options.year, symMCStats = (options.noSymMCStats))
                plotdir="postfit_plots/%s_LQ_m%i" % (fit_name,mLQ)
                print("\n plotdir = ", plotdir)
                print_and_do("[ -e %s ] && rm -r %s" % (plotdir, plotdir))
                print_and_do("mkdir %s" % (plotdir))
                #print_and_do("combine %s -M MultiDimFit  --saveWorkspace --saveFitResult --robustFit 1 %s --freezeParameters Afb " %(workspace, extra_params))
                print_and_do("combine %s -M MultiDimFit --saveWorkspace --saveFitResult --robustFit 1  %s " %(workspace, extra_params))

                if(not options.no_plot):
                    print_and_do("PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfit.root:fit_mdf --postfit -o %s_fit_shapes_LQ.root --sampling --samples 100"
                            % (fit_name))
                    extra_args = ""
                    if(options.year > 0): extra_args = " -y %i " % options.year
                    print_and_do("python scripts/LQ_plot_postfit.py -i %s_fit_shapes_LQ.root -o %s  %s --mLQ %i --chan %s --q %s --gen_level %s" % (fit_name, plotdir, extra_args,mLQ,options.chan,options.q,options.gen_level))
                    print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s  " % (workspace, extra_params)) #only to get prefit, probably a better way
                    print_and_do("python scripts/my_diffNuisances.py multidimfit.root --multidim --mLQ %i --prefit fitDiagnostics.root -p yLQ --skipFitB -g %s" % (mLQ, plotdir))
                    print_and_do("mv %s_fit_shapes_LQ.root %s" %(fit_name, plotdir))
                    if(not options.no_cleanup): print_and_do("rm fitDiagnostics.root higgsCombineTest.FitDiagnostics.mH120.root")


                print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
                print_and_do("""echo ".q" >> cmd.txt """)
                print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
                print_and_do("rm -f cards/sed*")
                if(not options.no_cleanup): print_and_do("rm cmd.txt combine_logger.out higgsCombineTest.MultiDimFit.mH120.root multidimfit.root")

