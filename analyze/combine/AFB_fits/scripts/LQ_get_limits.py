from LQ_utils import *

import ROOT
from ROOT import *

import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from itertools import product
import numpy as np
import json
from math import sqrt
from CombineHarvester.CombineTools.plotting import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def print_and_do(s):
    print("Exec: " + s)
    os.system(s)

def plotLimits(channel):
     # Style and pads
    ModTDRStyle()
    canv = ROOT.TCanvas('limit', 'limit')
    pads = OnePad()
     
     # Get limit TGraphs as a dictionary
    graphs = StandardLimitsFromJSONFile('LQ_cards/%s/limit_json/limits_%s.json'%(channel,channel))
    print(graphs)
    del graphs['obs']    
 # Create an empty TH1 from the first TGraph to serve as the pad axis and frame
    axis = CreateAxisHist(graphs.values()[0])
    if is_vec:
	if 'm' in channel: 
		axis.GetXaxis().SetTitle('m_{V_{#mu %s}} (GeV)'%(channel[0]))
		axis.GetYaxis().SetTitle('Limits on g_{#mu %s}'%(channel[0]))
	else: 
		axis.GetXaxis().SetTitle('m_{V_{e %s}} (GeV)'%(channel[0]))
		axis.GetYaxis().SetTitle('Limits on g_{e %s}'%(channel[0]))
    else:
	if 'm' in channel: 
		axis.GetXaxis().SetTitle('m_{S_{#mu %s}} (GeV)'%(channel[0]))
		axis.GetYaxis().SetTitle('Limits on y_{#mu %s}'%(channel[0]))
        else: 
		axis.GetXaxis().SetTitle('m_{S_{e %s}} (GeV)'%(channel[0]))
		axis.GetYaxis().SetTitle('Limits on y_{e %s}'%(channel[0]))
    pads[0].cd()
    axis.Draw('axis')
     
     # Create a legend in the top left
    legend = PositionedLegend(0.3, 0.2, 3, 0.015)
     
     # Set the standard green and yellow colors and draw
    StyleLimitBand(graphs)
    DrawLimitBand(pads[0], graphs, legend=legend)
    legend.Draw()
     
     # Re-draw the frame and tick marks
    pads[0].RedrawAxis()
    pads[0].GetFrame().Draw()
     
     # Adjust the y-axis range such that the maximum graph value sits 25% below
     # the top of the frame. Fix the minimum to zero.
    FixBothRanges(pads[0], 0, 0, GetPadYMax(pads[0]), 0.25)
     
     # Standard CMS logo
    DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.035, 1.2, '', 0.8)
     
    #canv.Print('.pdf')
    if(is_vec): canv.Print('LQ_cards/%s/limit_plots/limits_%s_vec_081922.png'%(channel,channel))
    else: canv.Print('LQ_cards/%s/limit_plots/limits_%s_081922.png'%(channel,channel))

parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
parser.add_option("--vec",  default=False, help="is vec?")
(options, args) = parser.parse_args()

mass = options.mLQ
is_vec = options.vec

extra_params=""
no_sys=False
fake_data=True
year = -1


print("nosys =%s"%(no_sys))
#make directory structure: LQ_cards/channel(eu,ed,mu,md)/masses 1000-3500

for channel in ['ue','de','um','dm']:   
#for channel in ['ce','se','cm','sm']:
#for channel in ['se','sm']:
    if channel=='ue' or channel=='ce':
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_ue.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_ue.txt"
    if channel=='de' or channel=='se':
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_de.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_de.txt"
    if channel=='um' or channel=='cm':
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_um.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_um.txt"
    if channel=='dm' or channel=='sm':
        if(no_sys): template_card = "card_templates/LQ_combined_fit_template_nosys_fake_dm.txt"
        if(fake_data): template_card = "card_templates/LQ_combined_fit_template_fake_dm.txt"

    #for mass in [1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000]:
    
    workspace ="LQ_cards/%s/%i/workspace.root"%(channel,mass)
    comb_card ="LQ_cards/%s/%i/combined_fit_%s_LQm%i.txt"%(channel,mass,channel,mass) 
    print_and_do("mkdir -p LQ_cards/%s/%i/"%(channel,mass))

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
        print_and_do("""sed -i "s/MASS/%i/g" %s""" % (mass, card))
        if not is_vec: print_and_do("""sed -i "s/QUARK/%s/g" %s""" % (channel[0], card))
        else: print_and_do("""sed -i "s/QUARK/%s_vec/g" %s""" % (channel[0], card))
        if(yr == 16 or yr == 17): print_and_do("""sed -i "s/#prefire/prefire/g" %s""" % (card))
       # if(yr == 18): print_and_do("""sed -i "s/#METHEM/METHEM/g" %s""" % (card))


    if(year < 0 ):
        print_and_do("combineCards.py Y16=cards/combined_fit_y16_LQ.txt Y17=cards/combined_fit_y17_LQ.txt Y18=cards/combined_fit_y18_LQ.txt > %s" % (comb_card))
    else:
        print_and_do("combineCards.py Y%i=cards/combined_fit_y%i_LQ.txt > %s" % (yr,yr,  comb_card))

    
    print("\n=========completed card for channel %s mass %i =========\n"%(channel,mass))
    print("\n========= making workspace for %s mass %i =========\n"%(channel,mass))
    print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s --channel-masks" % (comb_card, workspace))
    print("\n========= extracting upper limits for %s mass %i =========\n"%(channel, mass))
    #INCORRECT -> print_and_do("combineTool.py -d %s -M AsymptoticLimits -t -1  -m %i -n .limit --there"%(workspace,mass))
    print_and_do("combineTool.py -d %s -M AsymptoticLimits  -m %i -n .limit --there"%(workspace,mass))



    '''    
	print_and_do("mkdir LQ_cards/%s/limit_json/"%(channel))
	print_and_do("mkdir LQ_cards/%s/limit_plots/"%(channel))
    print("\n========= collecting limits for channel %s and making json =========\n"%(channel))
    print_and_do("combineTool.py -M CollectLimits LQ_cards/%s/*/*.limit.* --use-dirs -o LQ_cards/%s/limit_json/limits.json"%(channel,channel))
    
    with open("LQ_cards/%s/limit_json/limits_%s.json"%(channel,channel), 'r+') as f:
        data = json.load(f)
        for mass in ['1000.0','1500.0','2000.0','2500.0','3000.0','3500.0','4000.0','4500.0','5000.0','5500.0','6000.0','6500.0','7000.0','7500.0','8000.0','8500.0','9000.0']:
            for lim in data[mass]:
                yLQ2 = data[mass][lim]
                data[mass][lim] = sqrt(yLQ2)
        f.seek(0)        # <--- should reset file position to the beginning.
        json.dump(data, f, indent=4)
        f.truncate()     # remove remaining part
    
    print("\n========= making limit plot for channel %s =========\n"%(channel))
    #print_and_do("plotLimits.py LQ_cards/%s/limits_%s.json --auto-style exp"%(channel,channel))
    plotLimits(channel)
    '''

