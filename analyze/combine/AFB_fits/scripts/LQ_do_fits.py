
from LQ_utils import *
import ROOT
from ROOT import *
#import matplotlib.pyplot as plt
from array import array


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
parser.add_option("--noSymMCStats", default = True, action="store_true",  help="Don't add constraints to mcStat nuisances")
parser.add_option("--no_LQ",  default=False, action="store_true", help="For sanity check purposes remove LQ temps")
parser.add_option("--gen_level",  default=False, action="store_true", help="gen level fits")
(options, args) = parser.parse_args()



#for y in [2016,2017,2018]:
for y in [-1]:
    for options.chan in ["mumu"]:
    #for options.chan in ["mumu"]:
        for options.q in ["u","d"]:
            mLQ_list = [2500]
            #mLQ_list = [3500,4000,4500,5000]
	    is_vec = False
	    statuncs = False
	    #options.gen_level = False
            extra_params=""
            
#            options.q="u"
            print("options.gen_level = ",options.gen_level);
	    options.no_sys = False
            if not options.gen_level and not options.no_sys: options.fake_data = True
            options.no_LQ = False
            options.year = y
            if(options.verbose > 0):
        	extra_params +=" --verbose %i" % options.verbose

            
            fit_name = options.chan
            if(options.no_sys): 
                fit_name +="_nosys"
            if(not options.noSymMCStats):
                fit_name += "_SymMC"

            if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)
            fit_name+="_"+options.q

            if(options.no_LQ): fit_name+="_noLQ"

            if options.chan=="ee" and options.gen_level : fit_name+="_gen_level_SMdata_nlosys"
		
	    if is_vec: fit_name+="_vec"
	    if statuncs: fit_name += "_statuncs"
            fit_name+="_unblinded_bonly"
	    print("\n fit_name = ", fit_name)
	    #if y > -1: extra_args = "--combined False"
	    
            for mLQ in mLQ_list:
            #for mLQ in [1000]:
            #,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]:
            #mLQ = 1000.
            #for mbin in range(bin_start, bin_stop):
            #print(" \n \n Starting fit for bin %i \n\n" % mbin)
                
                print(" \n \n Starting fit for LQ m = %i\n\n",mLQ)
		
                workspace="workspaces/%s_LQ.root" % (options.chan)
                make_workspace(workspace, options.gen_level, options.chan, options.q, is_vec, options.no_LQ, options.no_sys, options.fake_data, mLQ, year = options.year,noSymMCStats = options.noSymMCStats)
                plotdir="postfit_plots/%s_LQ_m%i" % (fit_name,mLQ)
                print("\n plotdir = ", plotdir)
                #if not os.path.isdir(plotdir) or not os.listdir(plotdir):

		#print("Folder %s NOT FOUND"% plotdir)
		#make_workspace(workspace, options.gen_level, options.chan, options.q, is_vec, options.no_LQ, options.no_sys, options.fake_data, mLQ, year = options.year,noSymMCStats = options.noSymMCStats)
		print_and_do("rm -r %s" % (plotdir))
                print_and_do("mkdir %s" % (plotdir))
                if not statuncs:
		   print_and_do("combine %s -M MultiDimFit   --saveWorkspace --saveFitResult --robustFit 1 --trackErrors yLQ2 %s  -n .%s_%s%s_bonly_%i -s 3456 --setParameters yLQ2=0 --freezeParameters yLQ2" %(workspace, extra_params,options.chan,options.q,("_vec" if is_vec else ""),options.year))
                else:
		   print_and_do("combine %s -M MultiDimFit   --saveWorkspace --saveFitResult --robustFit 1  %s  -n .snapshot -s 3456" %(workspace, extra_params))
		   print_and_do("combine  -M MultiDimFit higgsCombine.snapshot.MultiDimFit.mH120.3456.root  --saveWorkspace --saveFitResult --robustFit 1  --freezeParameters allConstrainedNuisances --snapshotName MultiDimFit -s 3456")
                
		# higgsCombine.mumu_u_vec_2016.MultiDimFit.mH120.root
                if(not statuncs):
                    print_and_do("PostFitShapesFromWorkspace -w higgsCombine.%s_%s%s_%i.MultiDimFit.mH120.3456.root -f multidimfit.%s_%s%s_bonly_%i.root:fit_mdf --postfit -o %s_fit_shapes_LQ.root --sampling --samples 100"
                            % (options.chan,options.q,("_vec" if is_vec else ""),options.year,options.chan,options.q,("_vec" if is_vec else ""),options.year,fit_name))
                    extra_args = ""
                    if(options.year > 0): extra_args = " -y %i " % options.year
                    
                    print_and_do("python scripts/LQ_plot_postfit.py -i %s_fit_shapes_LQ.root -o %s  %s --mLQ %i --chan %s --q %s  %s" % (fit_name, plotdir, extra_args,mLQ,options.chan,options.q,("--vec True" if is_vec else "")))
                    #print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s  --robustFit 1 " % (workspace, extra_params)) #only to get prefit, probably a better way
                    print_and_do("python scripts/my_diffNuisances.py multidimfit.%s_%s%s_%i.root --multidim --mLQ %i --prefit fitDiagnosticsTest.root -p yLQ2 --skipFitB -g %s" % (options.chan,options.q,("_vec" if is_vec else ""),options.year,mLQ, plotdir))
                    print_and_do("mv %s_fit_shapes_LQ.root %s" %(fit_name, plotdir))
                    #if(not options.no_cleanup): print_and_do("rm fitDiagnosticsTest.root higgsCombineTest.FitDiagnostics.mH120.root")

		
                #print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
		print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
		print_and_do(""" cat cmd.txt """)
                print_and_do("""echo ".q" >> cmd.txt """)
                #print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
                print_and_do("root -l -b multidimfit.%s_%s%s_%i.root < cmd.txt > fit_results/%s_m%i.txt" % (options.chan,options.q,("_vec" if is_vec else ""),options.year,fit_name,mLQ))
		
		if(statuncs): print_and_do("root -l -b multidimfitTest.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
		
		print_and_do(""" echo "auto a=fit_mdf->floatParsFinal();" > cmd.txt """)
		print_and_do(""" echo "auto A0=(RooRealVar *) a.at(396);" >> cmd.txt """)
		print_and_do(""" echo "std::cout  << A0->getValV() << ' '  << A0->getErrorHi() << ' ' << A0->getErrorLo() << std::endl;" >> cmd.txt """)
		
                print_and_do(""" echo "auto Afb=(RooRealVar *) a.at(397);" >> cmd.txt """)
                print_and_do(""" echo "std::cout  << Afb->getValV() << ' '  << Afb->getErrorHi() << ' ' << Afb->getErrorLo() << std::endl;" >> cmd.txt """)
		if options.chan == "ee" and options.q != 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(398);" >> cmd.txt """)
		if options.chan == 'ee' and options.q == 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(192);" >> cmd.txt """) 
		if options.chan == "mumu" and options.q != 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(402);" >> cmd.txt """)
		if options.chan == "mumu" and options.q == 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(186);" >> cmd.txt """)
                print_and_do(""" echo "std::cout  << yLQ->getValV() << ' '  << yLQ->getErrorHi() << ' ' << yLQ->getErrorLo() << std::endl;" >> cmd.txt """)
		print_and_do("root -l -b multidimfit.%s_%s%s_%i.root < cmd.txt > %s/results_%s_m%i.txt" % (options.chan,options.q,("_vec" if is_vec else ""),options.year,plotdir,fit_name,mLQ))
                
                
		print_and_do("rm -f cards/sed*")
