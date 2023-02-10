
import ROOT
from ROOT import *
from LQ_utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=1000, type='int', help="Which mass bin to run on ")
parser.add_option("--nThreads",  default=10, type='int', help="Number of threads ")
parser.add_option("--expected",  default=False, action="store_true", help="Compute expected impacts based on toys with AFB=0.6 A0=0.05")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject if expected")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject if expected")
parser.add_option("-o", "--odir", default="impacts/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--diff", default = False, action="store_true",  help="Measure difference between electron and muon AFB's")

(options, args) = parser.parse_args()

for chan in ["mumu","ee"]:
    for q in ["s"]:

        options.mLQ = 2000
        fake_data = True
        no_sys = False
        gen_level = False
        no_LQ = False
        year = -1
	is_vec = False	
        extra_params = "--X-rtd MINIMIZER_no_analytic"        
        #extra_params = ""
	ending="020723_LQscaled10"
	if is_vec: ending += "_vec"
        if chan=="ee":
        #all_sys =   ["METJEC", "BTAGCOR","BTAGUNCOR", "BTAGLIGHT" , 
	    
            all_sys =   ["elScaleSyst", "elScaleStat","elScaleGain", "elSmear", "Pu",
                        "elHLTBARPTHIGH", "elIDBARPTHIGH", "elRECOBARPTHIGH", "elHLTENDPTHIGH", "elIDENDPTHIGH", "elRECOENDPTHIGH",
                        "elHLTBARPTLOW", "elIDBARPTLOW", "elRECOBARPTLOW", "elHLTENDPTLOW", "elIDENDPTLOW", "elRECOENDPTLOW",
                        #"ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b",
                        "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                        "elfakesrw1b", "elfakesrw2b", "elfakesrw3b", "elfakesrw4b",
                        "RENORM", "FAC", "REFAC", "alphaS","nlo_sys", 
                        "dy_xsec","db_xsec","top_xsec","gam_xsec" ,"elFakes",
                        "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC", "lumiGS" ,"lumi", 
                        ]

            correlate_all = ["elScaleSyst", "elSmear", "Pu", "elRECOBARPTHIGH", "elRECOENDPTHIGH", "elRECOBARPTLOW", "elRECOENDPTLOW",
                             "elIDBARPTHIGH", "elIDENDPTHIGH", "elIDBARPTLOW", "elIDENDPTLOW", 
                             "nlo_sys","dy_xsec","db_xsec"  ,"top_xsec","gam_xsec",
                             "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC" , "lumiGS",
                             ] 

            correlate_1718 = [#"ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b", 
                                "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                                "RENORM", "FAC", "REFAC", "alphaS"
			     ]
	    '''
	    all_sys = ["dy_xsec"]
	    correlate_all = ["dy_xsec"]
	    correlate_1718 = []
	    '''
        else:
	    
            all_sys =   [ "muPref","muRC", "Pu",
                        "muHLTBAR", "muIDBAR", "muISOBAR",  "muHLTEND", "muIDEND", "muISOEND",  "muIDSYS", "muISOSYS",  
                        #"ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b",
                        "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                        "mufakesrw1b", "mufakesrw2b", "mufakesrw3b", "mufakesrw4b",
                        "RENORM", "FAC", "REFAC", "alphaS", "nlo_sys",
                        "dy_xsec","db_xsec","top_xsec","gam_xsec" ,"muFakes", 
                        "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC", "lumiGS" ,"lumi", 
                        ]

            correlate_all = [ "Pu", "muIDSYS", "muISOSYS",
                             "nlo_sys","dy_xsec","db_xsec"  ,"top_xsec","gam_xsec",
                             "lumiXY" ,"lumiLS" ,"lumiDB" ,"lumiBC" , "lumiGS",
                             ] 

            correlate_1718 = ["muPref",
			#	"ptrw1b", "ptrw2b", "ptrw3b", "ptrw4b", "ptrw5b", "ptrw6b", "ptrw7b", 
                                "emucostrw1b", "emucostrw2b", "emucostrw3b", "emucostrw4b",
                                "RENORM", "FAC", "REFAC", "alphaS"
			     ]
	    '''
	    all_sys = ["dy_xsec"]
            correlate_all = ["dy_xsec"]
            correlate_1718 = []
	    '''
        pars16 = []
        pars17 = []
        pars18 = []
        pars_comb = []

        for i in range(1,61):
            pars_comb.append("pdf" + str(i))

        for par in all_sys:
            if(par not in correlate_all): 
                pars16.append(par + "16")
                if(par not in correlate_1718):
                    pars17.append(par + "17")
                    pars18.append(par + "18")
                else:
                    pars_comb.append(par + "1718")

            else:
                pars_comb.append(par)



	if chan=="ee":        pars16.append("prefire16")
	if chan=="ee":        pars17.append("prefire17")
        
	#pars18.append("METHEM18")


        par_str = ""

        for par in (pars16 + pars17 + pars18 + pars_comb):
            par_str += par +","

        #remove last comma
        par_str=par_str[:-1]

        ws_label = "%s_%s_impacts_mLQ%i" % (chan, q, options.mLQ)

        #if chan=="ee" and q=="u": yLQ_str = 'y_ue'
        #if chan=="ee" and q=="d": yLQ_str = 'y_ud'
        #if chan=="mumu" and q=="u": yLQ_str = 'y_um'
        #if chan=="mumu" and q=="d": yLQ_str = 'y_dm'

        yLQ2_str = 'yLQ2'
        A0_str = 'A0'
        Afb_str = 'Afb'
        if(options.diff):
            A0_str = 'dA0'
            Afb_str = 'dAfb'

        if(options.expected):
            ws_label = "%s_%s_expected_impacts_mbin%i" % (chan, q, options.mbin)




        workspace = "workspaces/%s.root" % (ws_label)
        #make_workspace(workspace, options.mbin, diff = options.diff)
        make_workspace(workspace, gen_level, chan, q, is_vec,  no_LQ , no_sys, fake_data, options.mLQ, year, False, False)

        print("Num pars = %i " % (len(pars16) + len(pars17) + len(pars18) + len(pars_comb)))
        print(par_str)

        s = 123456
        if(options.expected):
            print("Will inject AFB %.2f A0 %.2f for all toys " %(options.Afb, options.A0))

            #print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -n ToyGen -s %i --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters Afb=%.2f,A0=%.2f" 
            #    % (seed, options.Afb, options.A0))

            if(not options.reuse_fit):
                print_and_do("combine -M MultiDimFit -d %s --saveFitResult --saveWorkspace --robustFit 1 %s" % (workspace, extra_params))

            print_and_do("combine -M GenerateOnly -d higgsCombineTest.MultiDimFit.mH120.root -s %i --snapshotName MultiDimFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters Afb=%.2f,A0=%.2f" 
                    % (s, options.Afb, options.A0))


            extra_params += " --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist -t 1" % s


        print_and_do("combineTool.py -M Impacts -m 125 -d %s --doInitialFit --robustFit 1 %s " % (workspace, extra_params))

        if(options.expected):
            print_and_do("cp higgsCombine_initialFit_Test.MultiDimFit.mH125.%i.root higgsCombine_initialFit_Test.MultiDimFit.mH125.root" % s)

        print_and_do("combineTool.py -M Impacts -m 125 -d %s --doFits --named %s --parallel %i %s" % (workspace, par_str, options.nThreads, extra_params))

        if(options.expected):
            print("Renaming sys fits")
            for par in (pars16 + pars17 + pars18 + pars_comb):
                os.system("cp higgsCombine_paramFit_Test_%s.MultiDimFit.mH125.%i.root higgsCombine_paramFit_Test_%s.MultiDimFit.mH125.root" % (par, s, par))

        print_and_do("combineTool.py -M Impacts -m 125 -d %s -o %s/%s.json --named %s" % (workspace, options.odir, ws_label, par_str))
        print_and_do("python scripts/my_plotImpacts.py -i %s/%s.json -o %s/%s_plot_yLQ2_%s --POI %s --blind" % (options.odir, ws_label, options.odir, ws_label,ending, yLQ2_str))
        #print_and_do("python scripts/my_plotImpacts.py -i %s/%s.json -o %s/%s_plot_afb --POI %s --blind" % (options.odir, ws_label, options.odir, ws_label, Afb_str))
        #print_and_do("python scripts/my_plotImpacts.py -i %s/%s.json -o %s/%s_plot_a0 --POI %s --blind" % (options.odir, ws_label, options.odir, ws_label, A0_str))
        #print_and_do("rm higgsCombine_*")


