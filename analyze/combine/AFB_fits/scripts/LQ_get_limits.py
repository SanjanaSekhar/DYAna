from LQ_utils import *

import ROOT
from ROOT import *
from array import array
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from itertools import product
import numpy as np
import json
from math import sqrt
from CombineHarvester.CombineTools.plotting import *
import CMS_lumi, tdrstyle

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)

def print_and_do(s):
    print("Exec: " + s)
    os.system(s)

def plotLimits(channel):
     # Style and pads
    ModTDRStyle()
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadTickX(1)

    canv = ROOT.TCanvas('limit', 'limit')
    pads = OnePad()
     
     # Get limit TGraphs as a dictionary
    #graphs = StandardLimitsFromJSONFile('LQ_cards/%s/limit_json/limits_%s%s_%s.json'%(channel,channel,("_vec" if is_vec else ""),options.ending))
    graphs = StandardLimitsFromJSONFile('LQ_cards/%s/limit_json/limits_%s%s_%s_y%i.json'%(channel,channel, ("_vec" if is_vec else ""),options.ending,year-2000))
    print(graphs)
    #del graphs['obs']    
 # Create an empty TH1 from the first TGraph to serve as the pad axis and frame
    axis = CreateAxisHist(graphs.values()[0])
    #line_sp = ROOT.TLine(1000,1,1755,1)
    line_sp2 = ROOT.TLine(1755,1.,1755,5.3)
    x = array('d',[860,1175,1355,1755])
    y = array('d',[0.4,0.6,0.8,1.0])
    line_sp = ROOT.TGraph(4,x,y)

    line_pp = ROOT.TLine(1435,0.,1435,5.3)
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
    if not is_vec and (channel=='ue' or channel=='de'):
	line_sp.SetLineWidth(2)
	line_sp2.SetLineWidth(2)
	line_pp.SetLineWidth(2)
	line_sp.SetLineColor(kBlue)
	line_sp2.SetLineColor(kBlue)
	line_pp.SetLineColor(kRed+2)
	line_sp.SetLineStyle(9)
	line_sp2.SetLineStyle(9)
	line_pp.SetLineStyle(9)
	line_pp.Draw("same") 
	line_sp.Draw("same")
	line_sp2.Draw("same")
	legend.AddEntry(line_pp,"CMS Limit from arXiv:1811.01197","L")
    	legend.AddEntry(line_sp,"CMS Limit from arXiv:1509.03750","L")
    legend.Draw()
     
     # Re-draw the frame and tick marks
    pads[0].RedrawAxis()
    pads[0].GetFrame().Draw()
     # Adjust the y-axis range such that the maximum graph value sits 25% below
     # the top of the frame. Fix the minimum to zero.
    #FixBothRanges(pads[0], 0, 0, GetPadYMax(pads[0]), 0.25)
    FixBothRanges(pads[0], 0, 0, 4., 0.25)
     
     
# Standard CMS logo
    #DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.035, 1.2, '', 0.8)
    tdrstyle.setTDRStyle()
    CMS_lumi.CMS_lumi(pads[0], year, 11) 
    #canv.Print('.pdf')
    if(is_vec): canv.Print('LQ_cards/%s/limit_plots/limits_%s_vec_y%i_%s.png'%(channel,channel,year-2000, options.ending))
    else: canv.Print('LQ_cards/%s/limit_plots/limits_%s_y%i_%s.png'%(channel,channel,year-2000, options.ending))

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
parser.add_option("--vec",  default=False, help="is vec?")
parser.add_option("-o", "--odir", default="LQ_cards/condor/", help = "output directory")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--q",  default="u", help=" channel u,d,c,s ")
parser.add_option("--ending",  default="102022", help=" date ")
parser.add_option("--inject_yLQ2",  default=0.2, type='float', help="r=X")
parser.add_option("--quantile",  default=0.5, type='float', help="quantile expected")
parser.add_option("--ntoys",  default=10, type='int', help="no of toys")
parser.add_option("--iterations",  default=10, type='int', help="no of iterations")
parser.add_option("--hadd",  default=False, help="hadd")
parser.add_option("--HybridNew",  default=False, help="use HybridNew instead of AsymptoticLimits")
parser.add_option("--year",  default=-1,type='int', help="year")
(options, args) = parser.parse_args()


is_vec = options.vec
channel = options.q+options.chan[0]
mass = options.mLQ
extra_params=""
no_sys=False
fake_data=True
year = options.year


print("nosys =%s"%(no_sys))
#make directory structure: LQ_cards/channel(eu,ed,mu,md)/masses 1000-3500

#for channel in ['ue','de','um','dm']:   
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

#for mass in [1500]:#,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000]:

workspace ="LQ_cards/%s/%i/workspace_%s_%i.root"%(channel,mass,channel,mass)
#workspace = "workspaces/%s_LQ.root"%channel
comb_card ="LQ_cards/%s/%i/combined_fit_%s_LQm%i.txt"%(channel,mass,channel,mass) 
#comb_card ="cards/combined_fit_%s_LQm%i.txt"%(channel,mass)
print_and_do("rm LQ_cards/%s/%i/*"%(channel,mass))
print_and_do("mkdir -p LQ_cards/%s/%i/"%(channel,mass))

if(year > 0): years = [year % 2000]
else: years = [16,17,18]
nlo_sys = 0.3
#else: nlo_sys = 0.3 - 0.04*((mass-1500)/1000)
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
    print_and_do("""sed -i "s/NLO_SYS2/%f/g" %s""" % (1.+(nlo_sys*2), card))
    print_and_do("""sed -i "s/NLO_SYS4/%f/g" %s""" % (1.+(nlo_sys*4), card))
    if not is_vec: print_and_do("""sed -i "s/QUARK/%s/g" %s""" % (channel[0], card))
    else: print_and_do("""sed -i "s/QUARK/%s_vec/g" %s""" % (channel[0], card))
    if(yr == 16 or yr == 17): print_and_do("""sed -i "s/#prefire/prefire/g" %s""" % (card))
   # if(yr == 18): print_and_do("""sed -i "s/#METHEM/METHEM/g" %s""" % (card))


if(year < 0 ):
    print_and_do("combineCards.py Y16=cards/combined_fit_y16_LQ.txt Y17=cards/combined_fit_y17_LQ.txt Y18=cards/combined_fit_y18_LQ.txt > %s" % (comb_card))
else:
    print_and_do("combineCards.py Y%i=cards/combined_fit_y%i_LQ.txt > %s" % (yr,yr,  comb_card))
sigma = 0.6 **0.5
#extra_arg = " --symMCStats --sigma %f"%sigma 
extra_arg = ""
print("\n=========completed card for channel %s mass %i =========\n"%(channel,mass))
#print("\n========= making workspace for %s mass %i =========\n"%(channel,mass))
#print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s  %s" % (comb_card, workspace, extra_arg))

if options.hadd:

    if options.HybridNew:
	print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s  %s" % (comb_card, workspace, extra_arg))
        print_and_do("rm -rf HybridNew_output/%s/%s/"%(channel,mass))
        print_and_do("mkdir -p HybridNew_output/%s/%s/"%(channel,mass))
        print_and_do("rm eosoutput.log")
        for q_string in [0.025,0.160,0.500,0.840,0.975]:
            for point in np.arange(0.2,1.0,0.05):
                print_and_do("xrdfs root://cmseos.fnal.gov/ ls /store/user/sasekhar/Condor_outputs/limits_%s_%s_yLQ2%.2f_q%.3f_%s/ | grep POINT.%.3f | grep %.3f >> eosoutput.log"%(options.chan,options.q,point,q_string,options.ending,point,q_string))
                #print_and_do("eos root://cmseos.fnal.gov file rename  /store/user/sasekhar/Condor_outputs/limits_%s_%s_yLQ2%.2f_q%.3f_%s/higgsCombine.Test.POINT.%.6f.HybridNew.mH%i.*.quant%.3f.root /store/user/sasekhar/Condor_outputs/limits_%s_%s_yLQ2%.2f_q%.3f_%s/higgsCombine.yLQ2%.2f.HybridNew.mH%i.q%.3f.root "%(options.chan,options.q,point,q_string,options.ending,point,mass,q_string,options.chan,options.q,point,q_string,options.ending,point,mass,q_string))
            with open("eosoutput.log","r") as f:
                for line in f: 
                    print_and_do("xrdcp -f root://cmseos.fnal.gov/%s HybridNew_output/%s/%s/"
                    % ( line[:-1], channel, mass ))
            print_and_do("hadd merged_%s_q%.3f_m%i.root HybridNew_output/%s/%s/higgsCombine.Test.POINT.*.HybridNew.mH%i.*.quant%.3f.root"%(channel,q_string,mass,channel,mass,mass,q_string))
            print_and_do("combine %s -M HybridNew --LHCmode LHC-limits --readHybridResults --grid=merged_%s_q%.3f_m%i.root --expectedFromGrid %.3f"%(workspace,channel,q_string,mass,q_string))
    else:
	limits = {}
        for m in range(1000,5500,500):
            # /store/user/ssekhar/Condor_outputs/limits_ee_u_m9000_032823
            print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/limits_%s_%s%s_m%i_y%i_%s/limits_%s_m%i.json LQ_cards/%s/limit_json/limits_%s%s_m%i_y%i.json"
                %(options.chan, options.q, ("_vec" if is_vec else ""), m, year-2000,options.ending, channel, m, channel, channel, ("_vec" if is_vec else ""), m, year-2000))
            
	    with open("LQ_cards/%s/limit_json/limits_%s%s_m%i_y%i.json"%(channel,channel, ("_vec" if is_vec else ""),m, year-2000), 'r+') as f:
                data = json.load(f)
                #for mass in ['1000.0','1500.0','2000.0','2500.0','3000.0','3500.0','4000.0','4500.0','5000.0','5500.0','6000.0','6500.0','7000.0','7500.0','8000.0','8500.0','9000.0']:
                for lim in data[str(m)+".0"]:
                    yLQ2 = data[str(m)+".0"][lim]
                    data[str(m)+".0"][lim] = sqrt(yLQ2)
		    print(yLQ2,sqrt(yLQ2))
	 	limits[str(m)+".0"]=data[str(m)+".0"]
                f.seek(0)        # <--- should reset file position to the beginning.
                json.dump(data, f, indent=4)
                f.truncate()     # remove remaining part
	with open("LQ_cards/%s/limit_json/limits_%s%s_%s_y%i.json"%(channel,channel, ("_vec" if is_vec else ""),options.ending,year-2000), 'w') as f:
	    f.seek(0)
	    json.dump(limits, f, indent=4)
	    

        print("\n========= making limit plot for channel %s =========\n"%(channel))
        # #print_and_do("plotLimits.py LQ_cards/%s/limits_%s.json --auto-style exp"%(channel,channel))
        plotLimits(channel)

else:
    print("\n========= making workspace for %s mass %i =========\n"%(channel,mass))
    print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s  %s" % (comb_card, workspace, extra_arg))
    print("\n========= extracting upper limits for %s mass %i =========\n"%(channel, mass))
    

    if not options.HybridNew: 

        print_and_do("combineTool.py -d %s -M AsymptoticLimits  -m %i -n .limit --there -s -1 "%(workspace,mass)) #INCORRECT -> print_and_do("combineTool.py -d %s -M AsymptoticLimits -t -1  -m %i -n .limit --there"%(workspace,mass))
        print_and_do("mkdir LQ_cards/%s/limit_json/"%(channel))
        print_and_do("mkdir LQ_cards/%s/limit_plots/"%(channel))
        print("\n========= collecting limits for channel %s and making json =========\n"%(channel))
        print_and_do("combineTool.py -M CollectLimits LQ_cards/%s/%i/*.limit.* --use-dirs -o LQ_cards/%s/%i/limits.json"%(channel,mass,channel,mass))
        print_and_do("cp LQ_cards/%s/%i/limits_%s.json %s/limits_%s_m%i.json"%(channel,mass,channel,options.odir,channel,mass))
    
    else:
        print_and_do("combineTool.py %s -M HybridNew -H AsymptoticLimits --LHCmode LHC-limits -m %i --singlePoint %f --clsAcc 0 -s -1  --cminApproxPreFitTolerance 1.0 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerStrategy 0 -T %i -i %i  --X-rtd MINIMIZER_no_analytic --expectedFromGrid=%f --saveHybridResult "
            %(workspace,mass,options.inject_yLQ2,options.ntoys,options.iterations,options.quantile))
        print_and_do("cp *.root %s"%(options.odir))
   




