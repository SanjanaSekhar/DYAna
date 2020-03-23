
import ROOT
from ROOT import *
from utils import *

gROOT.SetBatch(True)


def gof_helper(chan, idx=0):
    f2 = TFile.Open("higgsCombineTest.GoodnessOfFit.mH120.root")

    t_data = f2.Get("limit")
    t_obs = t_data[0].limit
    f2.Close()

    f1 = TFile.Open("higgsCombineTest.GoodnessOfFit.mH120.123.root")


    toys = f1.Get("limit")

    toy_max = toys.GetMaximum("limit")
    toy_min = toys.GetMinimum("limit")

    my_max = max(toy_max, t_obs) + 10.
    #if is an error in fit can get very large values in toys
    my_max min(my_max, 2.*t_obs)

    my_min = min(toy_min, t_obs) - 5.

    h_test = TH1D("h_toys", "Goodness of fit: Mass bin %i" % idx, 30, my_min, my_max)
    toys.Draw("limit>>h_toys")
    bin_low = h_test.GetXaxis().FindBin(t_obs)
    bin_high = h_test.GetXaxis().FindBin(my_max)
    integral  = h_test.Integral()
    p_val = h_test.Integral(bin_low, bin_high) / integral

    print("Data gof is %.0f. p-value is %.3f based on %.0f toys" %(t_obs, p_val, integral))

    fout_name = "GoodnessOfFit/gof_%s_bin%i.png" % (chan, idx)

    c = TCanvas("c", "", 800, 800)
    h_test.Draw("hist")
    draw_max = h_test.GetMaximum()
    l = TLine(t_obs, 0., t_obs, draw_max)
    l.SetLineColor(kRed)
    l.SetLineWidth(2)
    l.Draw("same")

    latex = TLatex()
    latex.SetTextSize(0.025)
    latex.SetTextAlign(13)
    latex.SetNDC(True)
    latex.DrawLatex(0.75, 0.6, "Data gof is %.0f, p-value is %.2f" % (t_obs, p_val))
    c.Print(fout_name)
    f1.Close()
