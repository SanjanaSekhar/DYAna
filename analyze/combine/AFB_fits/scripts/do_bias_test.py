
import ROOT
from ROOT import *
from utils import *

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="", help = "output directory")

(options, args) = parser.parse_args()

chan = "combined"


workspace = "workspaces/%s_fit_bias_tests_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin)

print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1" % (workspace))
fitted_afb, fitted_a0 = setSnapshot(Afb_val = -1., mdf = True)

h_pull_afb = TH1F("h_pull_Afb", "", 20, -3, 3)
h_pull_a0 = TH1F("h_pull_A0", "", 20, -3, 3)

for i in range(options.nToys):
    print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root --snapshotName initialFit --toysFrequentist --bypassFrequentistFit -s %i --saveToys -t 1 --setParameters Afb=%.2f,A0=%.2f" 
            % (i, options.Afb, options.A0))
    #print_and_do("combine -M FitDiagnostics -d %s --saveWorkspace --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -t %i --robustFit 1   --forceRecreateNLL " % (workspace, options.nToys))
    print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root  -t 1 --robustFit 1 --forceRecreateNLL" 
            %(workspace, i))
    f_fit = TFile.Open("multidimfit.root")
    if f_fit:
        fr = f_fit.Get('fit_mdf')
        myargs = RooArgSet(fr.floatParsFinal())

        Afb = myargs.find("Afb").getVal()
        Afb_err = myargs.find("Afb").getError()
        A0 = myargs.find("A0").getVal()
        A0_err = myargs.find("A0").getError()
        print("Afb %.3f err %.3f " % (Afb, Afb_err))
        if(Afb_err > 0.): h_pull_afb.Fill((Afb- options.Afb)/ Afb_err)
        if(A0_err > 0.):  h_pull_a0.Fill((A0- options.A0)/ A0_err)
        f_fit.Close()
        print_and_do("rm higgsCombineTest.GenerateOnly.mH120.%i.root" % i)
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



c1 = TCanvas("c1", "", 900, 900)
h_pull_afb.Fit("gaus")
fit_afb= h_pull_afb.GetFunction("gaus")
fit_afb.SetLineColor(kBlue)
h_pull_afb.Draw()
h_pull_afb.SetTitle("Signal Inject Test Mass bin %i, Inject AFB %.2f" % (options.mbin, options.Afb))
h_pull_afb.GetXaxis().SetTitle("Pull Afb")
c1.Print("%sbias_test_mbin%i_afb%.0f.png" %(options.odir, options.mbin, 100.* options.Afb))


c2 = TCanvas("c1", "", 900, 900)
h_pull_a0.Fit("gaus")
fit_a0= h_pull_a0.GetFunction("gaus")
fit_a0.SetLineColor(kBlue)
h_pull_a0.Draw()
h_pull_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
h_pull_a0.GetXaxis().SetTitle("Pull A0")
c2.Print("%sbias_test_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))

