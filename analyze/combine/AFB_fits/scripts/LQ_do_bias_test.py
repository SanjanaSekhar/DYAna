
import ROOT
from ROOT import *
from LQ_utils import *
import numpy as np
from itertools import *

gROOT.SetBatch(True)

gStyle.SetOptFit(1) 




parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--yLQ",  default=0.0, type='float', help="yLQ value to inject")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--q",  default="u", help=" channel u,d,c,s ")
parser.add_option("-o", "--odir", default="signal_injection/condor/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--mLQ",  default=2000, type='int', help="mLQ")
parser.add_option("--is_vec", default=False, action="store_true", help="is vec?")
parser.add_option("--ending",  default="102022", help="png ext")
parser.add_option("--job",  default=1, type='int',help="job index")
parser.add_option("--plot",  default=False, help="make residuals and pulls")
(options, args) = parser.parse_args()



#options.nToys = 200
is_vec = options.is_vec
mLQ = options.mLQ
gen_level = False
no_LQ = False
fake_data = True
no_sys = False
year = -1
ending = options.ending
if is_vec: ending = "vec_"+ending

#if options.chan=="mumu" and options.q=="d": is_vec = True
yLQ2 = options.yLQ**2
yLQ = options.yLQ
print(options.chan,options.q)

if not options.plot:
    workspace = "workspaces/%s_%s_%i_fit_bias_tests.root" % (options.chan, options.q, options.job)
    make_workspace(workspace, gen_level, options.chan, options.q, is_vec, no_LQ , no_sys, fake_data, mLQ, year, True, False)

    #extra_params = "--X-rtd MINIMIZER_no_analytic"
    extra_params = ""

    print("Will inject options.yLQ %.2f for options.chan %s%s for all toys " %(options.yLQ,options.chan,options.q))

    res_yLQ = []
    pull_yLQ = []


    if(not options.prefit):
        print("Sampling toys based on postfit")
        if(not options.reuse_fit):
            print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s -s %i --freezeParameters A4,A0 --robustHesse=1" % (workspace, extra_params,123457+options.job))

    for i in range(options.nToys):

        if(not options.prefit):
            fitted_yLQ2 = setSnapshot(yLQ2_val = -1., mdf = True, s = 123457+options.job)
            if(options.no_sys):
                print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i  --snapshotName initialFit --toysNoSystematics --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ=%.2f --freezeParameters A4,A0 " % (i,yLQ))

            else:
                print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ=%.2f --freezeParameters A4,A0 " % (i,yLQ))
        else:

            print_and_do("combine -M GenerateOnly -d %s -s %i  --saveToys -t 1 --toysFrequentist --setParameters yLQ=%.2f --freezeParameters A4,A0 --robustHesse=1"% (workspace, i, yLQ))

        print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s -n _%i --freezeParameters A4,A0 --robustHesse=1" %(workspace,  i, extra_params, i))
        f_fit = TFile.Open("multidimfit_%i.root"%i)
        if f_fit:
            fr = f_fit.Get('fit_mdf')
            myargs = RooArgSet(fr.floatParsFinal())

            yLQ_fit = myargs.find("yLQ").getVal()
            yLQ_err = myargs.find("yLQ").getError()
            #A0 = myargs.find("A0").getVal()
            #A0_err = myargs.find("A0").getError()
            print("yLQ %.3f err %.3f " % (yLQ_fit, yLQ_err))
            res_yLQ.append(yLQ_fit - yLQ)
            if(yLQ_err > 0.): pull_yLQ.append((yLQ_fit-yLQ)/ yLQ_err)
            

            #h_res_afb.Fill(Afb - options.Afb)
            #h_res_a0.Fill(A0 - options.A0)
            #if(Afb_err > 0.): h_pull_afb.Fill((Afb- options.Afb)/ Afb_err)
            #if(A0_err > 0.):  h_pull_a0.Fill((A0- options.A0)/ A0_err)
            f_fit.Close()
            #print_and_do("rm higgsCombineTest.GenerateOnly.mH120.%i.root" % i)
        else:
            print("Can't open multidimfit. Looks like the fit failed.")



    #f_ws = TFile.Open("higgsCombineTest.FitDiagnostics.mH120.123456.root")
    #f_toys = TFile.Open("higgsCombineTest.GenerateOnly.mH120.123456.root")
    #ws = f_ws.Get('w')
    #toy = f_toys.Get('toys/toy_1')
    #getattr(ws, 'import')(toy)
    #ws.writeToFile("toy_ws.root")
    #print_and_do("PostFitShapesFromWorkspace -w toy_ws.root --dataset model_sData  -f fitDiagnostics.root:fit_s -o toy_shapes.root --sampling --samples 100")
    ##print_and_do("python scripts/plot_postfit.py -i toy_ws.root -o test/ -m %i" % (mbin))
    with open('%srespull_%s_%s_%i_yLQ%.1f_%s.txt'%(options.odir,options.chan,options.q,options.job,options.yLQ,ending), 'w') as f:
        for res,pull in zip(res_yLQ, pull_yLQ):
            f.write("%f %f\n" %(res,pull))

else:

    respull = []

    for job_idx in [0,1,2,3,4,5,6,7,8,9]:
        print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/sasekhar/Condor_outputs/bias_test_yLQ%.1f_%s_%s%s_m%s_%i_%s/respull_%s_%s_%i_yLQ%.1f_%s.txt %s%s" % (options.yLQ, options.chan, options.q, ("_vec" if is_vec else ""), options.mLQ, job_idx, ending[-6:],options.chan,options.q,job_idx,options.yLQ,ending, options.odir,options.mLQ))
        with open('%s%s/respull_%s_%s_%i_yLQ%.1f_%s.txt'%(options.odir,options.mLQ,options.chan,options.q,job_idx,options.yLQ,ending), 'r') as f:
            for line in f.readlines():
        	respull.append(line.split(' '))

    respull = np.asarray(respull, dtype=float)
    res_yLQ2 = respull[:,0].tolist()
    pull_yLQ2 = respull[:,1].tolist()

    print(res_yLQ2)

    n_bins = 20
    h_pull_yLQ2 = TH1F("h_pull_yLQ2", "", n_bins, -3.5, 3.5)

    res_yLQ2_range = max(3.5*np.std(res_yLQ2), 0.15)

    h_res_yLQ2 = TH1F("h_res_yLQ2", "", n_bins, -res_yLQ2_range, res_yLQ2_range)



    def fill_h(arr, h):
        #print(arr)
        for arg in arr:
            h.Fill(arg)

    fill_h(pull_yLQ2, h_pull_yLQ2)
    #fill_h(pull_a0, h_pull_a0)
    fill_h(res_yLQ2, h_res_yLQ2)
    #fill_h(res_a0, h_res_a0)

    if options.chan=="mumu": chan_label = "#mu"
    else : chan_label = "e"

    c1 = TCanvas("c1", "", 900, 900)
    h_pull_yLQ2.Fit("gaus")
    fit_yLQ2= h_pull_yLQ2.GetFunction("gaus")
    if(fit_yLQ2): fit_yLQ2.SetLineColor(kBlue)
    h_pull_yLQ2.Draw()

    if is_vec: 
        h_pull_yLQ2.SetTitle("Signal Inject Test : Inject g_{%s %s} = %.1f (M_{LQ} = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
        h_pull_yLQ2.GetXaxis().SetTitle("Pull g_{%s %s}"%(chan_label,options.q))
        c1.Print("%s%s/bias_test_pull_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.mLQ,options.yLQ, options.chan, options.q,ending))
    else: 
        h_pull_yLQ2.SetTitle("Signal Inject Test : Inject y_{%s %s} = %.1f (M_{LQ} = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
        h_pull_yLQ2.GetXaxis().SetTitle("Pull y_{%s %s}"%(chan_label,options.q))
        c1.Print("%s%s/bias_test_pull_yLQ%.1f_%s_%s_%s.png" %(options.odir, options.mLQ,options.yLQ, options.chan, options.q,ending))


    # c2 = TCanvas("c1", "", 900, 900)
    # h_pull_a0.Fit("gaus")
    # fit_a0= h_pull_a0.GetFunction("gaus")
    # if(fit_a0): fit_a0.SetLineColor(kBlue)
    # h_pull_a0.Draw()
    # h_pull_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
    # h_pull_a0.GetXaxis().SetTitle("Pull A0")
    # c2.Print("%sbias_test_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))


    c3 = TCanvas("c3", "", 900, 900)
    h_res_yLQ2.Fit("gaus")
    fit_yLQ2= h_res_yLQ2.GetFunction("gaus")
    if(fit_yLQ2): fit_yLQ2.SetLineColor(kBlue)
    h_res_yLQ2.Draw()
    if is_vec: 
        h_res_yLQ2.SetTitle("Signal Inject Test : Inject g_{%s %s} = %.1f (M_{LQ} = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
        h_res_yLQ2.GetXaxis().SetTitle("#Delta g_{%s %s}"%(chan_label,options.q))
        c3.Print("%s%s/bias_test_res_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.mLQ, options.yLQ, options.chan, options.q,ending))
    else: 
        h_res_yLQ2.SetTitle("Signal Inject Test : Inject y_{%s %s} = %.1f (M_{LQ} = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
        h_res_yLQ2.GetXaxis().SetTitle("#Delta y_{%s %s}"%(chan_label,options.q))
        c3.Print("%s%s/bias_test_res_yLQ%.1f_%s_%s_%s.png" %(options.odir,options.mLQ, options.yLQ, options.chan, options.q,ending))



    # c4 = TCanvas("c1", "", 900, 900)
    # h_res_a0.Fit("gaus")
    # fit_a0= h_res_a0.GetFunction("gaus")
    # if(fit_a0): fit_a0.SetLineColor(kBlue)
    # h_res_a0.Draw()
    # h_res_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
    # h_res_a0.GetXaxis().SetTitle("#Delta A0")
    # c4.Print("%sbias_test_res_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))

