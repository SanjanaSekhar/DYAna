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

def setSnapshot(mdf = False, yLQ2_val = 0.0,  d=''):
    fit_name = 'fit_s'
    workspace = d+'higgsCombineTest.FitDiagnostics.mH120.root'
    fit_file = d+'fitDiagnosticsTest.root'
    if(mdf):
        fit_name = 'fit_mdf'
        workspace = d+'higgsCombineTest.MultiDimFit.mH120.root'
        fit_file = d+'multidimfitTest.root'
    w_f = TFile.Open(workspace)
    w = w_f.Get('w')
    fr_f = TFile.Open(fit_file)
    fr = fr_f.Get(fit_name)
    myargs = RooArgSet(fr.floatParsFinal())
    fitted_yLQ2 = myargs.find("yLQ2").getVal()
    #fitted_a0 = myargs.find("A0").getVal()
    if(yLQ2_val > -0.9):

        myargs.find("yLQ2").setVal(yLQ2_val)
        #myargs.find("Afb").setError(0.)
        #myargs.find("A0").setVal(A0_val)
        #myargs.find("A0").setError(0.)
    # end new lines
    importPars = w.saveSnapshot('initialFit',myargs, True)
    fout = TFile('initialFitWorkspace.root',"recreate")
    fout.WriteTObject(w,'w')
    fout.Close()
    return fitted_yLQ2

def do_lumi(card, year):
        l_vals = dict()
        l_vals['LYRV'] = [1.022, 1.02, 1.015]
        l_vals['LXYV'] = [1.009, 1.009, 1.02]
        l_vals['LLSV'] = [-1,1.003, 1.002] 
        l_vals['LDBV'] = [1.004, 1.004, -1]
        l_vals['LBCV'] = [-1, 1.003, 1.002]
        l_vals['LGSV'] = [1.004, 1.001, -1]

        idx = year - 16


        for key,l_val in l_vals.items():
            val = l_val[idx]
            f_string = ""
            if(val > 0):
                f_string = "%.3f" % val
            else:
                f_string = "  -  "

            print_and_do("""sed -i "s/%s/%s/g" %s""" % (key, f_string, card))

        print_and_do("""sed -i "s/LUMIYR/LUMI%i/g" %s""" % (year, card))

def make_workspace(workspace, gen_level, chan, q, is_vec = False, no_LQ = False, no_sys = False, fake_data = False, mLQ = 1000, year = -1,symMCStats = False):
    print("\n inside make_workspace()")
    print("Making workspace %s LQ" % (workspace))
    print("nosys =%s"%(no_sys))

    #template_card="card_templates/LQ_combined_fit_template_nosys_fake.txt"
    
    if chan=="ee" and (q=="u" or q=="c"):
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_ue.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_ue.txt"
        if(fake_data) and no_LQ: template_card = "card_templates/LQ_combined_fit_template_fake_ee_noLQ.txt"
    if chan=="ee" and (q=="d" or q=="s"):
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_de.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_de.txt"
        if(fake_data) and no_LQ: template_card = "card_templates/LQ_combined_fit_template_fake_ee_noLQ.txt"
    if chan=="mumu" and (q=="u" or q=="c"):
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_um.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_um.txt"
        if(fake_data) and no_LQ: template_card = "card_templates/LQ_combined_fit_template_fake_mumu_noLQ.txt"
    if chan=="mumu" and (q=="d" or q=="s"):
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_dm.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_dm.txt"
        if(fake_data) and no_LQ: template_card = "card_templates/LQ_combined_fit_template_fake_mumu_noLQ.txt"
    if gen_level and q=="u":
        template_card = "card_templates/LQ_combined_fit_template_genlevel_ue.txt"
    if gen_level and q=="d":
        template_card = "card_templates/LQ_combined_fit_template_genlevel_de.txt"
   
    #comb_card="cards/combined_fit_mbin%i.txt" % mbin
    comb_card = "cards/combined_fit_LQ.txt" 
    if gen_level: comb_card = "cards/combined_fit_genlevel_LQ.txt" 

    if(year > 0): years = [year % 2000]
    else: years = [16,17,18]

    for yr in years:
        if(yr == 16):
            comb_yr = 16
        else:
            #some systematics combined between 17 and 18
            comb_yr = 1718
        card="cards/combined_fit_y%i_LQ.txt" % (yr)
        print_and_do("cp %s %s" % (template_card, card))
        do_lumi(card, yr)
        print_and_do("""sed -i "s/YRC/%i/g" %s""" % (comb_yr, card))
        print_and_do("""sed -i "s/YR/%i/g" %s""" % (yr, card))
        print_and_do("""sed -i "s/MASS/%i/g" %s""" % (mLQ, card))
        if not is_vec: print_and_do("""sed -i "s/QUARK/%s/g" %s""" % (q, card))
        else: print_and_do("""sed -i "s/QUARK/%s_vec/g" %s""" % (q, card))
        if(yr == 16 or yr == 17): print_and_do("""sed -i "s/#prefire/prefire/g" %s""" % (card))
        #if(yr == 18): print_and_do("""sed -i "s/#METHEM/METHEM/g" %s""" % (card))


    if(year < 0 ):
        print_and_do("combineCards.py Y16=cards/combined_fit_y16_LQ.txt Y17=cards/combined_fit_y17_LQ.txt Y18=cards/combined_fit_y18_LQ.txt > %s" % (comb_card))
    else:
        print_and_do("combineCards.py Y%i=cards/combined_fit_y%i_LQ.txt > %s" % (yr,yr,  comb_card))

    if(symMCStats): extra_arg = "--symMCStats"
    else: extra_arg = ""
    print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s --channel-masks %s" % (comb_card, workspace, extra_arg))
    #print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq -o %s --channel-masks %s" % (comb_card, workspace, extra_arg))

def make_gen_level_workspace(workspace, mbin,  year = -1):
    print("Making workspace %s mbin %i" % (workspace, mbin))
    template_card="card_templates/gen_level_fit_template.txt"

    yr = year % 2000

    card="cards/gen_level_fit_y%i_mbin%i.txt" % (yr, mbin)
    print_and_do("cp %s %s" % (template_card, card))

    print_and_do("""sed -i "s/YR/%i/g" %s""" % (yr, card))

    print_and_do("text2workspace.py %s --keyword-value M_BIN=%i -P Analysis.DYAna.my_model:dy_AFB -o %s" % (card, mbin, workspace))
