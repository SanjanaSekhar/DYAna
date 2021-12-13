import ROOT
from ROOT import *
from utils import *
from parse import parse


def parse_log(filename):
    f = open(filename, "r")
    lines = f.readlines()
    data_gof = -1
    toys_gofs = []

    for line in lines:
        if "Best fit test statistic" in line: 
            val = float(parse("Best fit test statistic: {}", line)[0])
            if(data_gof < 0):
                data_gof = val
            else:
                toys_gofs.append(val)
    return data_gof, np.array(toys_gofs)


gROOT.SetBatch(True)
parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-i", "--input",  default="", help="Logfile")
parser.add_option("-o", "--output",  default="", help="Output")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")


(options, args) = parser.parse_args()

chan = "combined"


data_gof, toys_gofs = parse_log(options.input)

t_obs = data_gof


toy_max = np.amax(toys_gofs)
toy_min = np.amin(toys_gofs)

my_max = max(toy_max, t_obs)*1.2
#if is an error in fit can get very large values in toys
my_max =  min(my_max, 4.*t_obs)

my_min = min(toy_min, t_obs)*0.8

h_test = TH1D("h_toys", "Goodness of fit (saturated): Mass bin %i" % (options.mbin), 30, my_min, my_max)
for toy in toys_gofs:
    h_test.Fill(toy)

np_toys = toys_gofs

bin_low = h_test.GetXaxis().FindBin(t_obs)
bin_high = h_test.GetXaxis().FindBin(my_max)
integral  = h_test.Integral()
if(integral < 1.):
    print("Integral toy gof distribution  is 0? Maybe fits failed?" )
    exit(1)


p_val = h_test.Integral(bin_low, bin_high) / integral

np_above = np_toys > t_obs
tot = np_toys.shape[0]
n_above = np_toys[np_above].shape[0]
np_p_val = float(n_above) / tot

print("Data gof is %f p-value (integral) is %.3f based on %.0f toys" %(t_obs, p_val, integral))
print("Numpy version: %i out of %i above t_obs (p-val %.3f) " % (n_above, tot, float(n_above) / tot))

p_val = np_p_val


draw_max = h_test.GetMaximum()
c = TCanvas("c", "", 800, 800)
h_test.Draw("hist")
l = TLine(t_obs, 0., t_obs, draw_max)
l.SetLineColor(kRed)
l.SetLineWidth(2)
l.Draw("same")

latex = TLatex()
latex.SetTextSize(0.025)
latex.SetTextAlign(13)
latex.SetNDC(True)
latex.DrawLatex(0.5, 0.75, "Data gof is %.1f, p-value is %.3f" % (t_obs, p_val))
c.Print(options.output)
