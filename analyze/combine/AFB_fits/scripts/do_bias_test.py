import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product

gROOT.SetBatch(True)

nToys = 1
chan = "combined"

inject_AFB = 0.6
inject_A0 = 0.05
mbin = 0

comb_card = "cards/%s_fit_mbin%i.txt" % (chan, mbin)
workspace = "workspaces/%s_fit_bias_tests_%i.root" % (chan, mbin)
os.system("text2workspace.py %s --keyword-value M_BIN=%i -P Analysis.DYAna.my_model:dy_AFB -o %s --channel-masks" % (comb_card, mbin, workspace))
os.system("combine -M GenerateOnly -d %s --toysNoSystematics -t %i --saveToys --setParameters Afb=%.2f,A0=%.2f" % (workspace, nToys, inject_AFB, inject_A0))
os.system("combine -M FitDiagnostics -d %s --saveWorkspace --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root  -t %i --skipBOnlyFit " %(workspace, nToys))

#os.system("combine -M FitDiagnostics %s --saveWorkspace --skipBOnlyFit -t %i --toysNoSystematics --robustFit 1 --setParameters Afb=%.2f,A0=%.2f -v 2" % (workspace, nToys, inject_AFB, inject_A0))
exit(1)


f_fit = TFile.Open("fitDiagnostics.root")
tree_fit = f_fit.Get("tree_fit_sb")

h_pull_afb = TH1F("h_pull_afb", "", 20, -4, 4)
h_pull_a0 = TH1F("h_pull_a0", "", 20, -4, 4)

tree_fit.Draw("(Afb - %.2f)/AfbErr>>h_pull_afb" % inject_AFB)
tree_fit.Draw("(A0 - %.2f)/A0Err>>h_pull_a0" % inject_A0)

h_pull_afb.Fit("gaus")
h_pull_a0.Fit("gaus")
