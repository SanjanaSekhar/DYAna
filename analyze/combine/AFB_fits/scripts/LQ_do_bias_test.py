
import ROOT
from ROOT import *
from LQ_utils import *
import numpy as np
from itertools import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--yLQ",  default=0.0, type='float', help="yLQ value to inject")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--q",  default="u", help=" channel u,d,c,s ")
parser.add_option("-o", "--odir", default="signal_injection/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")

(options, args) = parser.parse_args()

gStyle.SetOptFit(1) 


options.nToys = 200
is_vec = False
mLQ = 1500
gen_level = False
no_LQ = False
fake_data = True
no_sys = False
year = -1
ending = "062422"


if options.chan=="mumu" and options.q=="d": is_vec = True
yLQ2 = options.yLQ**2
print(options.chan,options.q)
workspace = "workspaces/%s_%s_fit_bias_tests.root" % (options.chan, options.q)
make_workspace(workspace, gen_level, options.chan, options.q, is_vec, no_LQ , no_sys, fake_data, mLQ, year,False)

#extra_params = "--X-rtd MINIMIZER_no_analytic"
extra_params = ""

print("Will inject options.yLQ %.2f for options.chan %s%s for all toys " %(options.yLQ,options.chan,options.q))

res_yLQ2 = []
pull_yLQ2 = []


if(not options.prefit):
    print("Sampling toys based on postfit")
    if(not options.reuse_fit):
        print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s" % (workspace, extra_params))

for i in range(options.nToys):

    if(not options.prefit):
        fitted_yLQ2 = setSnapshot(yLQ2_val = -1., mdf = True)
        if(options.no_sys):
            print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root  --snapshotName initialFit --toysNoSystematics --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ2=%.2f" 
                % (yLQ2))

        else:
            print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root  --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters yLQ2=%.2f" 
            % (yLQ2))
    else:

        print_and_do("combine -M GenerateOnly -d %s  --saveToys -t 1 --toysFrequentist --setParameters yLQ2=%.2f"% (workspace,  yLQ2))

    print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.root --toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s" 
            %(workspace,  extra_params))
    f_fit = TFile.Open("multidimfit.root")
    if f_fit:
        fr = f_fit.Get('fit_mdf')
        myargs = RooArgSet(fr.floatParsFinal())

        yLQ2_fit = myargs.find("yLQ2").getVal()
        yLQ2_err = myargs.find("yLQ2").getError()
        #A0 = myargs.find("A0").getVal()
        #A0_err = myargs.find("A0").getError()
        print("yLQ2 %.3f err %.3f " % (yLQ2_fit, yLQ2_err))
        res_yLQ2.append(yLQ2_fit - yLQ2)
        if(yLQ2_err > 0.): pull_yLQ2.append((yLQ2_fit-yLQ2)/ yLQ2_err)
        

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
with open('%srespull_%s_%s_yLQ%.1f_%s.txt'%(options.odir,options.chan,options.q,options.yLQ,ending), 'w') as f:
    for res,pull in zip(res_yLQ2,pull_yLQ2):
        f.write("%f %f\n" %(res,pull))

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

if options.chan=="mumu": chan_label = "\mu"
else : chan_label = "e"

c1 = TCanvas("c1", "", 900, 900)
h_pull_yLQ2.Fit("gaus")
fit_yLQ2= h_pull_yLQ2.GetFunction("gaus")
if(fit_yLQ2): fit_yLQ2.SetLineColor(kBlue)
h_pull_yLQ2.Draw()

if is_vec: 
    h_pull_yLQ2.SetTitle(r"Signal Inject Test : Inject $g_%s %s$ = %.1f ($M_{LQ}$ = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
    h_pull_yLQ2.GetXaxis().SetTitle(r"Pull $g_{%s %s}^2$"%(chan_label,options.q))
    c1.Print("%sbias_test_pull_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.yLQ, options.chan, options.q,ending))
else: 
    h_pull_yLQ2.SetTitle(r"Signal Inject Test : Inject $y_%s %s$ = %.1f ($M_{LQ}$ = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
    h_pull_yLQ2.GetXaxis().SetTitle(r"Pull $y_{%s %s}^2$"%(chan_label,options.q))
    c1.Print("%sbias_test_pull_yLQ%.1f_%s_%s_%s.png" %(options.odir, options.yLQ, options.chan, options.q,ending))


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
    h_res_yLQ2.SetTitle(r"Signal Inject Test : Inject $g_%s %s$ = %.1f ($M_{LQ}$ = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
    h_res_yLQ2.GetXaxis().SetTitle(r"#Delta $g_{%s %s}^2$"%(chan_label,options.q))
    c3.Print("%sbias_test_res_yLQ%.1f_%s_%s_vec_%s.png" %(options.odir, options.yLQ, options.chan, options.q,ending))
else: 
    h_res_yLQ2.SetTitle(r"Signal Inject Test : Inject $y_%s %s$ = %.1f ($M_{LQ}$ = %.1f TeV)" % (chan_label,options.q,options.yLQ,mLQ/1000.))
    h_res_yLQ2.GetXaxis().SetTitle(r"#Delta $y_{%s %s}^2$"%(chan_label,options.q))
    c3.Print("%sbias_test_res_yLQ%.1f_%s_%s_%s.png" %(options.odir, options.yLQ, options.chan, options.q,ending))



# c4 = TCanvas("c1", "", 900, 900)
# h_res_a0.Fit("gaus")
# fit_a0= h_res_a0.GetFunction("gaus")
# if(fit_a0): fit_a0.SetLineColor(kBlue)
# h_res_a0.Draw()
# h_res_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
# h_res_a0.GetXaxis().SetTitle("#Delta A0")
# c4.Print("%sbias_test_res_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))

