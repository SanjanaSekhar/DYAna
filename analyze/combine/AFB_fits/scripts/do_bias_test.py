import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from numpy import arange
from itertools import product

gROOT.SetBatch(True)

nToys = 3
chan = "combined"

inject_AFB = 0.6
inject_A0 = 0.05
mbin = 0

template_card="card_templates/combined_fit_template.txt"
#comb_card="cards/combined_fit_mbin%i.txt" % mbin
comb_card = "cards/%s_fit_mbin%i.txt" % (chan, mbin)
for year in (16,17,18):
    card="cards/combined_fit_y%i_mbin%i.txt" % (year, mbin)
    os.system("cp %s %s" % (template_card, card))
    os.system("""sed -i "s/YR/%i/g" %s""" % (year, card))

os.system("combineCards.py Y16=cards/combined_fit_y16_mbin%i.txt Y17=cards/combined_fit_y17_mbin%i.txt Y18=cards/combined_fit_y18_mbin%i.txt > %s" % (mbin, mbin, mbin, comb_card))



workspace = "workspaces/%s_fit_bias_tests_%i.root" % (chan, mbin)
os.system("text2workspace.py %s --keyword-value M_BIN=%i -P Analysis.DYAna.my_model:dy_AFB -o %s --channel-masks" % (comb_card, mbin, workspace))
os.system("combine -M GenerateOnly --toysFrequentist -d %s -t %i --saveToys --setParameters Afb=%.2f,A0=%.2f" % (workspace, nToys, inject_AFB, inject_A0))
#os.system("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root  -t %i --robustFit 1" %(workspace, nToys))
os.system("combine -M FitDiagnostics -d %s --saveWorkspace --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -t %i --robustFit 1" % (workspace, nToys))

#f_ws = TFile.Open("higgsCombineTest.FitDiagnostics.mH120.123456.root")
#f_toys = TFile.Open("higgsCombineTest.GenerateOnly.mH120.123456.root")
#ws = f_ws.Get('w')
#toy = f_toys.Get('toys/toy_1')
#getattr(ws, 'import')(toy)
#ws.writeToFile("toy_ws.root")
#os.system("PostFitShapesFromWorkspace -w toy_ws.root --dataset model_sData  -f fitDiagnostics.root:fit_s -o toy_shapes.root --sampling --samples 100")
##os.system("python scripts/plot_postfit.py -i toy_ws.root -o test/ -m %i" % (mbin))

#exit(1)


#f_fit = TFile.Open("multidimfit.root")
#tree_fit = f_fit.Get("fit_mdf")
f_fit = TFile.Open("fitDiagnostics.root")
tree_fit = f_fit.Get("tree_fit_sb")

h_pull_afb = TH1F("h_pull_afb", "", 20, -4, 4)
#h_pull_a0 = TH1F("h_pull_a0", "", 20, -4, 4)

tree_fit.Draw("(Afb - %.2f)/AfbErr>>h_pull_afb" % inject_AFB)
#tree_fit.Draw("(A0 - %.2f)/A0Err>>h_pull_a0" % inject_A0)

h_pull_afb.Fit("gaus")
#h_pull_a0.Fit("gaus")
