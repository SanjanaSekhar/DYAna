
from LQ_utils import *
import ROOT
from ROOT import *
import matplotlib.pyplot as plt
from array import array
plt.ioff()


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



for y in [-1]:
    for options.chan in ["mumu"]:
    	#for options.chan in ["ee"]:
        for options.q in ["u"]:
            #mLQ_list = [500,1000,2000,3000]
            mLQ_list = [2000]
	    is_vec = False
	    statuncs = False
	    #options.gen_level = False
            extra_params=""
#            options.chan="mumu"
#            options.q="u"
            print("options.gen_level = ",options.gen_level);
	    options.no_sys = False
            if not options.gen_level and not options.no_sys: options.fake_data = True
            options.no_LQ = False
            options.year = y
            likelihood_scan = False
	    if likelihood_scan: 
		poi = 'nlo_sys'
		ending = "freezeA0A4_%s"%poi
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
            #extra_params += "  --freezeParameters A4,A0 "
	    #extra_params += " --cminApproxPreFitTolerance 1.0 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerStrategy 0 "
	    if statuncs: extra_params += " --freezeParameters allConstrainedNuisances"
            
            fit_name = options.chan
            if(options.no_sys): 
                fit_name +="_nosys"
                #extra_params += " --freezeParameters allConstrainedNuisances"
            if(not options.noSymMCStats):
                fit_name += "_SymMC"
            if(options.fake_data): fit_name +="_fake_data"

            if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)
            fit_name+="_"+options.q

            if(options.no_LQ): fit_name+="_noLQ"

            if options.chan=="ee" and options.gen_level : fit_name+="_gen_level_SMdata_nlosys"
		
	    if is_vec: fit_name+="_vec"
	    if statuncs: fit_name += "_statuncs"
            fit_name+=""
	    print("\n fit_name = ", fit_name)
	    
	    
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
                
		print_and_do("[ -e %s ] && rm -r %s" % (plotdir, plotdir))
                print_and_do("mkdir %s" % (plotdir))
                print_and_do("combine %s -M MultiDimFit   --saveWorkspace --saveFitResult --robustFit 1 --trackErrors yLQ2 %s  --cminDefaultMinimizerStrategy 0" %(workspace, extra_params))
                #print_and_do("combine %s -M MultiDimFit --saveWorkspace --saveFitResult --robustFit 1  %s " %(workspace, extra_params))
                if likelihood_scan: print_and_do("combine %s -M MultiDimFit --algo grid --points 200 --squareDistPoiStep  --autoRange 2 -P %s --floatOtherPOIs 1 --saveWorkspace --saveFitResult --robustFit 1  %s " %(workspace, poi,  extra_params))

                if(not options.no_plot):
                    print_and_do("PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfitTest.root:fit_mdf --postfit -o %s_fit_shapes_LQ.root --sampling --samples 100"
                            % (fit_name))
                    extra_args = ""
                    if(options.year > 0): extra_args = " -y %i " % options.year
                    print_and_do("python scripts/LQ_plot_postfit.py -i %s_fit_shapes_LQ.root -o %s  %s --mLQ %i --chan %s --q %s --vec %s" % (fit_name, plotdir, extra_args,mLQ,options.chan,options.q,is_vec))
                    #print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s  --robustFit 1 " % (workspace, extra_params)) #only to get prefit, probably a better way
                    #print_and_do("python scripts/my_diffNuisances.py multidimfitTest.root --multidim --mLQ %i --prefit fitDiagnosticsTest.root -p yLQ2 --skipFitB -g %s" % (mLQ, plotdir))
                    print_and_do("mv %s_fit_shapes_LQ.root %s" %(fit_name, plotdir))
                    #if(not options.no_cleanup): print_and_do("rm fitDiagnosticsTest.root higgsCombineTest.FitDiagnostics.mH120.root")

		
                #print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
		print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
		print_and_do(""" cat cmd.txt """)
                print_and_do("""echo ".q" >> cmd.txt """)
                #print_and_do("root -l -b multidimfit.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
                print_and_do("root -l -b multidimfitTest.root < cmd.txt > fit_results/%s_m%i.txt" % (fit_name,mLQ))
		
		print_and_do(""" echo "auto a=fit_mdf->floatParsFinal();" > cmd.txt """)
		print_and_do(""" echo "auto A0=(RooRealVar *) a.at(0);" >> cmd.txt """)
		print_and_do(""" echo "std::cout  << A0->getValV() << ' '  << A0->getErrorHi() << ' ' << A0->getErrorLo() << std::endl;" >> cmd.txt """)
		
                print_and_do(""" echo "auto Afb=(RooRealVar *) a.at(1);" >> cmd.txt """)
                print_and_do(""" echo "std::cout  << Afb->getValV() << ' '  << Afb->getErrorHi() << ' ' << Afb->getErrorLo() << std::endl;" >> cmd.txt """)
		if options.chan == "ee" and options.q != 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(348);" >> cmd.txt """)
		if options.chan == 'ee' and options.q == 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(192);" >> cmd.txt """) 
		if options.chan == "mumu" and options.q != 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(342);" >> cmd.txt """)
		if options.chan == "mumu" and options.q == 's': print_and_do(""" echo "auto yLQ=(RooRealVar *) a.at(186);" >> cmd.txt """)
                print_and_do(""" echo "std::cout  << yLQ->getValV() << ' '  << yLQ->getErrorHi() << ' ' << yLQ->getErrorLo() << std::endl;" >> cmd.txt """)
		print_and_do("root -l -b multidimfitTest.root < cmd.txt > %s/results_%s_m%i.txt" % (plotdir,fit_name,mLQ))
                
                
		print_and_do("rm -f cards/sed*")
                if likelihood_scan: print_and_do("cp higgsCombineTest.MultiDimFit.mH120.root higgsCombineTest.MultiDimFit._%s_%s.root"%(options.chan,options.q))
                #if(not options.no_cleanup): print_and_do("rm cmd.txt combine_logger.out higgsCombineTest.MultiDimFit.mH120.root multidimfit.root")
                 
                if likelihood_scan:

                    deltaNLL, yLQ2_list = [],[]
		    
                    f = ROOT.TFile.Open("higgsCombineTest.MultiDimFit._%s_%s.root"%(options.chan,options.q),"READ")
                    limit_tree = f.Get("limit")
		    poi_value = array('f',[0])
		    limit_tree.SetBranchAddress("%s"%poi, poi_value)

                    for i in range(limit_tree.GetEntries()):
			
                        limit_tree.GetEntry(i)
                        deltaNLL.append(limit_tree.deltaNLL)
                        yLQ2_list.append(poi_value[0])
			print(poi_value)
			#if limit_tree.quantileExpected > 0.49 and limit_tree.quantileExpected < 0.51: print(limit_tree.yLQ,limit_tree.deltaNLL)
			#if limit_tree.quantileExpected > 0.83 and limit_tree.quantileExpected < 0.86: print(limit_tree.yLQ,limit_tree.deltaNLL)	
                    f.Close()
		    #print(yLQ2_list)
		    #print(deltaNLL)
		    idx = np.argsort(np.array(yLQ2_list))
		    yLQ2_list = np.array(yLQ2_list)[idx]
		    deltaNLL = np.array(deltaNLL)[idx]
                    with open('like_scan_%s_%s_m%i_%s.txt'%(options.chan,options.q,mLQ,ending), 'w') as f:
	     		for ylq,dnll in zip(yLQ2_list, deltaNLL):
         		    f.write("%f %f\n" %(ylq,2*dnll)) 
		    #print(np.amax(yLQ2_list),np.amin(yLQ2_list))
	    if likelihood_scan: 
	        for mLQ in mLQ_list:	    
	            respull = []
	            with open('like_scan_%s_%s_m%i_%s.txt'%(options.chan,options.q,mLQ,ending), 'r') as f:
    	                for line in f.readlines():
                            respull.append(line.split(' '))

	            respull = np.asarray(respull, dtype=float)
	            yLQ2_list = respull[:,0].tolist()
	            deltaNLL = respull[:,1].tolist()
	            plt.ylim(0,10)		    
	            plt.plot(yLQ2_list,deltaNLL,label='mLQ=%s GeV'%mLQ)
            	plt.xlabel("%s"%poi)
            	plt.ylabel("-2deltaLL")
            	plt.legend()
            	plt.title("Likelihood Scan: channel %s %s"%(options.chan,options.q))
            	plt.savefig("like_scan_%s_%s_%s.jpg"%(options.chan,options.q,ending))
            	plt.close()
		    
