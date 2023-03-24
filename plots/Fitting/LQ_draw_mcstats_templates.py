import ROOT
from ROOT import *
from optparse import OptionParser
import sys

ext = "020923"
mLQ = 2500
odir = "AN_plots/Systematics/UpDown/"
x_start = 0.7;
x_end = 0.9;
y_start = 0.75;
y_end = 0.9;

for year in [2016,2017,2018]:

    fin = "../analyze/combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)

    f = ROOT.TFile.Open(fin, "READ")
    f.cd("LQ")
    keys = ROOT.gDirectory.GetListOfKeys().Clone()
    ee_sys_list = [s for s in keys if 'ee%i_fpl_MCStatBin'%(year-2000) in s]
    mumu_sys_list = [s for s in keys if 'mumu%i_fpl_MCStatBin'%(year-2000) in s]
    for sys in ee_sys_list:
        if 'Up' in sys:
            h_up = ROOT.gDirectory.Get(sys)
            h_down = ROOT.gDirectory.Get(sys[:-2]+"Down")
            h = ROOT.gDirectory.Get("ee%i_fpl"%(year-2000))
            c1 = TCanvas("c1","c1",200,300,700,500)
            c1.cd()
            h.SetLineWidth(2)
            h_up.SetLineWidth(2)
            h_down.SetLineWidth(2)
            h.SetLineColor(kBlue)
            h_up.SetLineColor(kRed)
            h_down.SetLineColor(kGreen+3)
            h.SetTitle(sys[:-2])
            h.Draw("hist")
            h_up.Draw("hist same")
            h_down.Draw("hist same")
            leg1 = TLegend(x_start, y_start, x_end, y_end);
            leg1.AddEntry(h, "Nominal Template", "l");
            leg1.AddEntry(h_up, "Sys Up Template", "l");
            leg1.AddEntry(h_down, "Sys Down Template", "l");
            c1.Print("%s/ee_fpl_MCStatBin.pdf"%odir)

    for sys in mumu_sys_list:
        if 'Up' in sys:
            h_up = ROOT.gDirectory.Get(sys)
            h_down = ROOT.gDirectory.Get(sys[:-2]+"Down")
            h = ROOT.gDirectory.Get("mumu%i_fpl"%(year-2000))
            c1 = TCanvas("c1","c1",200,300,700,500)
            c1.cd()
            h.SetLineWidth(2)
            h_up.SetLineWidth(2)
            h_down.SetLineWidth(2)
            h.SetLineColor(kBlue)
            h_up.SetLineColor(kRed)
            h_down.SetLineColor(kGreen+3)
            h.SetTitle(sys[:-2])
            h.Draw("hist")
            h_up.Draw("hist same")
            h_down.Draw("hist same")
            leg1 = TLegend(x_start, y_start, x_end, y_end);
            leg1.AddEntry(h, "Nominal Template", "l");
            leg1.AddEntry(h_up, "Sys Up Template", "l");
            leg1.AddEntry(h_down, "Sys Down Template", "l");
            c1.Print("%s/mumu_fpl_MCStatBin.pdf"%odir)


