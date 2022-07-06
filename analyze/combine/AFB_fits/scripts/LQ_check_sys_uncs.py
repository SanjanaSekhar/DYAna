import operator
import ROOT
from ROOT import *
from LQ_utils import *
from add_group_impact import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--expected",  default=False, action="store_true", help="Compute expected impacts based on toys with AFB=0.6 A0=0.05")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject if expected")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject if expected")
parser.add_option("-o", "--odir", default="sys_uncs/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--diff", default=False, action="store_true", help="Diff")

(options, args) = parser.parse_args()

for chan in ["mumu"]:
    for q in ["u"]:
        options.mLQ = 1500
        fake_data = True
        no_sys = False
        gen_level = False
        no_LQ = False
        year = -1
	is_vec = False
        #extra_params = "--X-rtd MINIMIZER_no_analytic"
	ending = "062422"
        s = 123456
        extra_params = " -s %i" % s

        if chan == "ee":
            individual_pars = ["nlo_sys", "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "elFakesYR",  "Pu", "prefireYR"]
            group_pars =[  "RFscalesYRC", "emucostrwsYRC", "ptrwsYRC", "pdfs", "lumisYR","elScalesYR", "elHLTsYR", "elIDs", "elRECOs",  
                            "elfakesrwsYR", "autoMCStats"] 
            #"BTAGSYR","muPrefYRC","METJECYR",
        else:
            individual_pars = ["nlo_sys", "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "muFakesYR", "Pu", "muPrefYRC",  "muRCYR", ]
            group_pars =[  "RFscalesYRC", "emucostrwsYRC", "ptrwsYRC", "pdfs", "lumisYR","muIDsYR", "muHLTsYR", 
                            "mufakesrwsYR",  "autoMCStats"] 

        sys_name_conv = dict()
	sys_name_conv['nlo_sys'] = "LQ LO reweighting"
        sys_name_conv['dy_xsec'] = "DY Cross Section"
        sys_name_conv['db_xsec'] = "Diboson Cross Section"
        sys_name_conv['top_xsec'] = "$\\ttbar$ Cross Section"
        sys_name_conv['gam_xsec'] = "$\\PGg\\PGg$ Cross Section"
        sys_name_conv['RFscalesYRC'] = "$\\alpha_s$ + Renormalization/Factorization Scales"
        sys_name_conv['emucostrwsYRC'] = "$e\\mu$ Shape Corrections"
        sys_name_conv['ptrwsYRC'] = "DY $p_{T}$ Correction"
        sys_name_conv['pdfs'] = "PDFs"
        sys_name_conv['lumisYR'] = "Luminosity"
        sys_name_conv['autoMCStats'] = "MC and MisID Backgrounds Statistical Uncertainty"
        sys_name_conv['Pu'] = "Pileup"

        if chan == "ee": 
            sys_name_conv['elFakesYR'] = "Electron MisID Normalization"
            sys_name_conv['elScalesYR'] = "Electron Momentum Scale"
            sys_name_conv['elHLTsYR'] = "Electron Trigger"
            sys_name_conv['elIDs'] = "Electron Identification/Isolation"
            sys_name_conv['elRECOs'] = "Electron Reconstruction"
            sys_name_conv['elfakesrwsYR'] = "Electron MisID Shape"
            sys_name_conv['prefireYR'] = "Electron Trigger Prefire Correction"

        if chan== "mumu": 
            sys_name_conv['muFakesYR'] = "Muon MisID Normalization"
            sys_name_conv['muRCYR'] = "Muon Momentum Scale"
            sys_name_conv['muIDsYR'] = "Muon Identification/Isolation"
            sys_name_conv['muHLTsYR'] = "Muon Trigger"
            sys_name_conv['muPrefYRC'] = "Muon Trigger Prefire Correction"
            sys_name_conv['mufakesrwsYR'] = "Muon MisID Shape"
            sys_name_conv['muPrefYRC'] = "Muon Trigger Prefire Correction"


        #sys_name_conv['METJECYR'] = "MET Uncertainties"
        #sys_name_conv['BTAGSYR'] = "b-tagging Uncertainty"





        def par_to_freezestr(par):
            if("YRC" in par):
                par16 = par.replace("YRC", "16")
                par1718 = par.replace("YRC", "1718")
                return par16 + "," + par1718
            elif("YR" in par):
    
		if  "prefire" not in par:
                    par16 = par.replace("YR", "16")
                    par17 = par.replace("YR", "17")
                    par18 = par.replace("YR", "18")
                    return par16 + "," + par17 + "," + par18
                else: #prefire only sys with no 18
                    if chan=="ee":   return "prefire16,prefire17"
            else:
                return par

                





        workspace = "workspaces/%s_%s_sys_uncs.root" % (chan, q)


        make_workspace(workspace, gen_level, chan, q, is_vec, no_LQ , no_sys, fake_data, options.mLQ, year,False)
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

        print_and_do("""combine -M FitDiagnostics -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _nom1 %s -v 2""" % (s,extra_params))




        d = dict()
        n = 0
        for indi_par in individual_pars:
            n+=1
            freeze_str = par_to_freezestr(indi_par)
            print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _%s %s -v 2""" 
                    % (freeze_str,s, indi_par, extra_params))
            sys_unc = compute_sys("nom1", indi_par, -1)
            #sys_unc = compute_sys("nom", indi_par, s)
            d[indi_par] = sys_unc
            #if(n>2): break

        for group_par in group_pars:
            n+=1
            freeze_str = par_to_freezestr(group_par)
            print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _%s %s -v 2""" % 
                    (freeze_str, s, group_par, extra_params))
            sys_unc = compute_sys("nom1", group_par, -1)
            #sys_unc = compute_sys("nom", indi_par, s)
            d[group_par] =sys_unc
            #if(n>4): break

        print(d)
        os.system("mkdir %s \n" % options.odir)
        with open("%s/%s_%s_m%s_sys_uncs_%s.txt" % (options.odir, chan, q, options.mLQ,ending), 'w') as f_out:
            sorted_d = sorted(d.items(), key=operator.itemgetter(1))
            f_out.write("Systematic uncertainties (values x1000) for bin %i \n" % options.mbin)
            for sys_name, val in sorted_d[::-1]:
                if (sys_name in sys_name_conv.keys()):
                        out_name = sys_name_conv[sys_name]
                else:
                    out_name = sys_name
                f_out.write("%s & %.2f \\\\ \n" % (out_name, val*1000.))


        #print_and_do("rm higgsCombine* fitDiagnostics* multidimfit*")

