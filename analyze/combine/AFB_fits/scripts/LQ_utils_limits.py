import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from itertools import product
import numpy as np


def gof_helper(chan, mbin=0, odir = "GoodnessOfFit/"):
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

    my_max = max(toy_max, t_obs) + 10.
    #if is an error in fit can get very large values in toys
    my_max =  min(my_max, 2.*t_obs)

    my_min = min(toy_min, t_obs) - 5.

    h_test = TH1D("h_toys", "Goodness of fit: Mass bin %i" % mbin, 30, my_min, my_max)
    for toy in toys:
        h_test.Fill(toy.limit)

    bin_low = h_test.GetXaxis().FindBin(t_obs)
    bin_high = h_test.GetXaxis().FindBin(my_max)
    integral  = h_test.Integral()
    if(integral < 1.):
        print("Integral toy gof distribution  is 0? Maybe fits failed?" )
        exit(1)


    p_val = h_test.Integral(bin_low, bin_high) / integral

    print("Data gof is %.0f. p-value is %.3f based on %.0f toys" %(t_obs, p_val, integral))

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
    latex.DrawLatex(0.5, 0.75, "Data gof is %.0f, p-value is %.2f" % (t_obs, p_val))
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



def make_workspace(no_sys = False, fake_data = False,  year = -1):
    print("\n inside make_workspace()")
    #print("Making workspace %s LQ" % (workspace))
    print("nosys =%s"%(no_sys))
    #make directory structure: LQ_cards/channel(eu,ed,mu,md)/masses 1000-3500
    #template_card="card_templates/LQ_combined_fit_template_nosys_fake.txt"
    for channel in ['eu','ed','mu','md']:

        if channel=='eu':
            if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_ue.txt"
            if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_ue.txt"
        if channel=='ed':
            if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_de.txt"
            if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_de.txt"
        if channel=='mu':
            if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_um.txt"
            if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_um.txt"
        if channel=='md':
            if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_dm.txt"
            if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_dm.txt"
   
        for mass in [1500,2000,2500,3000,3500]:
        #comb_card="cards/combined_fit_mbin%i.txt" % mbin
            comb_card ="LQ_cards/%s/%i/combined_fit_%s_LQm%i.txt"%(channel,mass,channel,mass) 
            print_and_do("mkdir -p LQ_cards/%s/%i/"%(channel,mass))

            if(year > 0): years = [year % 2000]
            else: years = [16,17,18]

            for yr in years:
                card="cards/combined_fit_y%i_LQ.txt" % (yr)
                print_and_do("cp %s %s" % (template_card, card))
                print_and_do("""sed -i "s/YR/%i/g" %s""" % (yr, card))
                if(yr == 16 or yr == 17): print_and_do("""sed -i "s/#prefire/prefire/g" %s""" % (card))
                if(yr == 18): print_and_do("""sed -i "s/#METHEM/METHEM/g" %s""" % (card))


            if(year < 0 ):
                print_and_do("combineCards.py Y16=cards/combined_fit_y16_LQ.txt Y17=cards/combined_fit_y17_LQ.txt Y18=cards/combined_fit_y18_LQ.txt > %s" % (comb_card))
            else:
                print_and_do("combineCards.py Y%i=cards/combined_fit_y%i_LQ.txt > %s" % (yr,yr,  comb_card))

            print_and_do("""sed -i "s/cards/%s/g" %s""" % ("..", comb_card))
            print("\ncompleted card for channel %s mass %i\n",channel,mass)
   
        print("\n=========making workspace for %s=========\n",channel)
    
        print_and_do("combineTool.py -M T2W  -P LQ_Analysis.DYAna.LQ_my_model:dy_AFB -i LQ_cards/%s/* -o %i/workspace.root "%(channel,mass))
    #print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:dy_AFB -o %s --channel-masks" % (comb_card, workspace))
    #print_and_do("combineTool.py -M Asymptotic -d ed/*/workspace.root --there -n .limit")

