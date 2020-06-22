
import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from itertools import product
import numpy as np


def gof_helper(chan, mbin=0, odir = "GoodnessOfFit/", teststat = 'saturated'):
    ROOT.ROOT.EnableImplicitMT()
    f2 = TFile.Open("higgsCombineTest.GoodnessOfFit.mH120.root")

    t_data = f2.Get("limit")
    #array = np.squeeze(t_data.AsMatrix(columns=['limit']))
    array = [ev.limit for ev in t_data]
    print(array)
    t_obs = array[0]
    f2.Close()

    f1 = TFile.Open("higgsCombineTest.GoodnessOfFit.mH120.123456.root")


    toys = f1.Get("limit")

    toy_max = toys.GetMaximum("limit")
    toy_min = toys.GetMinimum("limit")

    my_max = max(toy_max, t_obs)*1.2
    #if is an error in fit can get very large values in toys
    my_max =  min(my_max, 4.*t_obs)

    my_min = min(toy_min, t_obs)*0.8

    h_test = TH1D("h_toys", "Goodness of fit (%s): Mass bin %i" % (teststat, mbin), 30, my_min, my_max)
    np_toys = np.array([])
    for toy in toys:
        np_toys = np.append(np_toys, toy.limit)
        h_test.Fill(toy.limit)

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


    fout_name = "%sgof_%s_bin%i.png" % (odir, chan, mbin)
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
    latex.DrawLatex(0.5, 0.75, "Data gof is %.3f, p-value is %.2f" % (t_obs, p_val))
    c.Print(fout_name)
    f1.Close()


def print_and_do(s):
    print("Exec: " + s)
    os.system(s)

def setSnapshot(mdf = False, Afb_val = 0.6, A0_val= 0.05, d=''):
    fit_name = 'fit_s'
    workspace = d+'higgsCombineTest.FitDiagnostics.mH120.root'
    fit_file = d+'fitDiagnostics.root'
    if(mdf):
        fit_name = 'fit_mdf'
        workspace = d+'higgsCombineTest.MultiDimFit.mH120.root'
        fit_file = d+'multidimfit.root'
    w_f = TFile.Open(workspace)
    w = w_f.Get('w')
    fr_f = TFile.Open(fit_file)
    fr = fr_f.Get(fit_name)
    myargs = RooArgSet(fr.floatParsFinal())
    fitted_afb = myargs.find("Afb").getVal()
    fitted_a0 = myargs.find("A0").getVal()
    if(Afb_val > -0.9):

        myargs.find("Afb").setVal(Afb_val)
        #myargs.find("Afb").setError(0.)
        myargs.find("A0").setVal(A0_val)
        #myargs.find("A0").setError(0.)
    # end new lines
    importPars = w.saveSnapshot('initialFit',myargs)
    fout = TFile('initialFitWorkspace.root',"recreate")
    fout.WriteTObject(w,'w')
    fout.Close()
    return fitted_afb, fitted_a0

def make_workspace(workspace, mbin, no_sys = False, fake_data = False, year = -1):
    print("Making workspace %s mbin %i" % (workspace, mbin))
    template_card="card_templates/combined_fit_template.txt"
    if(no_sys): template_card = "card_templates/combined_fit_template_nosys.txt"
    if(fake_data): template_card = "card_templates/combined_fit_template_fake_data.txt"
    #comb_card="cards/combined_fit_mbin%i.txt" % mbin
    comb_card = "cards/combined_fit_mbin%i.txt" % (mbin)

    if(year > 0): years = [year % 2000]
    else: years = [16,17,18]

    for yr in years:
        card="cards/combined_fit_y%i_mbin%i.txt" % (yr, mbin)
        print_and_do("cp %s %s" % (template_card, card))
        print_and_do("""sed -i "s/YR/%i/g" %s""" % (yr, card))
        if(yr == 16 or yr == 17): print_and_do("""sed -i "s/#prefire/prefire/g" %s""" % (card))
        if(yr == 18): print_and_do("""sed -i "s/#METHEM/METHEM/g" %s""" % (card))


    if(year < 0 ):
        print_and_do("combineCards.py Y16=cards/combined_fit_y16_mbin%i.txt Y17=cards/combined_fit_y17_mbin%i.txt Y18=cards/combined_fit_y18_mbin%i.txt > %s" % (mbin, mbin, mbin, comb_card))
    else:
        print_and_do("combineCards.py Y%i=cards/combined_fit_y%i_mbin%i.txt > %s" % (yr,yr, mbin, comb_card))

    print_and_do("text2workspace.py %s --keyword-value M_BIN=%i -P Analysis.DYAna.my_model:dy_AFB -o %s --channel-masks" % (comb_card, mbin, workspace))
