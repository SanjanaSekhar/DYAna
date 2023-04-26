#for combine review
from LQ_utils import *
import ROOT
from ROOT import *
import matplotlib.pyplot as plt
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from itertools import product
import numpy as np

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=2000, type='int', help="mLQ")
parser.add_option("--vec",  default=False, help="is vec?")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--q",  default="u", help=" channel u,d,c,s ")
parser.add_option("-o", "--odir", default="combine_review/", help = "output directory")
parser.add_option("--hadd",  default=False, help="hadd")
(options, args) = parser.parse_args()

year = -1
mLQ = options.mLQ
chan = options.chan
q = options.q
is_vec = options.vec

# for is_vec in [True]:
# 	for chan in ["ee"]:
# 		for q in ["d"]:

workspace="workspaces/%s_%s_%s_LQ.root" % (chan,q,("vec" if is_vec else ""))

if chan=="ee" and (q=="u" or q=="c"):
	template_card = "card_templates/LQ_combined_fit_template_fake_ue.txt"
if chan=="ee" and (q=="d" or q=="s"):
	template_card = "card_templates/LQ_combined_fit_template_fake_de.txt"
if chan=="mumu" and (q=="u" or q=="c"):
	template_card = "card_templates/LQ_combined_fit_template_fake_um.txt"
if chan=="mumu" and (q=="d" or q=="s"):
	template_card = "card_templates/LQ_combined_fit_template_fake_dm.txt"



comb_card = "cards/EXO-22-013_%s_%s_%s.txt"%(chan,q,("vec" if is_vec else "")) 


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

#sigma = 0.6 **0.5
#extra_arg = "--symMCStats --sigma %f"%sigma	
extra_arg = ""	

if not options.hadd:
	
	print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s %s" % (comb_card, workspace, extra_arg))
	#print_and_do("combine -M FitDiagnostics -d %s -t -1 --setParameters yLQ2=0.0  --forceRecreateNLL -n _t0" %workspace)
	print_and_do("combine -M MultiDimFit -d %s -t -1 --setParameters yLQ2=0.0  --forceRecreateNLL -n _t0 --saveFitResult --robustFit 1" %workspace)	
	print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  -a multidimfit_t0.root -p yLQ2  -g plots_t0.root >> ./fitResults_t0_%s_%s_%s"%(chan,q,("vec" if is_vec else ""))) 

	#print_and_do("combine -M FitDiagnostics -d %s -t -1 --setParameters yLQ2=0.6  --forceRecreateNLL -n _t1" %workspace)
	print_and_do("combine -M MultiDimFit -d %s -t -1 --setParameters yLQ2=0.6  --forceRecreateNLL -n _t1 --saveFitResult --robustFit 1" %workspace)
	print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a multidimfit_t1.root -p yLQ2  -g plots_t1.root >> ./fitResults_t1_%s_%s_%s"%(chan,q,("vec" if is_vec else "")))
	'''
	print_and_do("combineTool.py -M Impacts -d %s -t -1 --setParameters yLQ2=0.0 -m 2000 --doInitialFit --allPars -n t0"%workspace)
	print_and_do("combineTool.py -M Impacts -d %s -t -1 --setParameters yLQ2=0.6 -m 2000 --doInitialFit --allPars -n t1"%(workspace))

	print_and_do("combineTool.py -M Impacts -d %s -o impacts_t0_%s_%s_%s.json -t -1 --setParameters yLQ2=0.0 --doFits -m 2000 -n t0 "%(workspace,chan,q,("vec" if is_vec else "")))

	print_and_do("combineTool.py -M Impacts -d %s -o impacts_t1_%s_%s_%s.json -t -1 --setParameters yLQ2=0.6 --doFits -m 2000 -n t1 "%(workspace,chan,q,("vec" if is_vec else "")))

	print_and_do("combineTool.py -M Impacts -d %s  -m 2000 -n t0 -o impacts_t0_%s_%s_%s.json"%(workspace,chan,q,("vec" if is_vec else "")))
	print_and_do("combineTool.py -M Impacts -d %s  -m 2000 -n t1 -o impacts_t1_%s_%s_%s.json"%(workspace,chan,q,("vec" if is_vec else "")))
	print_and_do("plotImpacts.py -i  impacts_t0_%s_%s_%s.json -o  impacts_t0_%s_%s_%s"%(chan,q,("vec" if is_vec else ""),chan,q,("vec" if is_vec else "")))
	print_and_do("plotImpacts.py -i  impacts_t1_%s_%s_%s.json -o  impacts_t1_%s_%s_%s"%(chan,q,("vec" if is_vec else ""),chan,q,("vec" if is_vec else "")))
	'''
	print_and_do("cp fitResults_t0_%s_%s_%s %s"%(chan,q,("vec" if is_vec else ""),options.odir))
	print_and_do("cp fitResults_t1_%s_%s_%s %s"%(chan,q,("vec" if is_vec else ""),options.odir))
	print_and_do("cp impacts_t0_%s_%s_%s* %s"%(chan,q,("vec" if is_vec else ""),options.odir))
	print_and_do("cp impacts_t1_%s_%s_%s* %s"%(chan,q,("vec" if is_vec else ""),options.odir))

else:

	print_and_do("cp %s ."%(comb_card))
	print_and_do("ValidateDatacards.py %s --jsonFile combine_review/validation_%s_%s_%s.json"%(comb_card[6:],options.chan, options.q, ("vec" if is_vec else "")))
	print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/CR_%s_%s%s/fitResults_t0_%s_%s_%s combine_review/"
                %(options.chan, options.q, ("_vec" if is_vec else ""),options.chan, options.q, ("vec" if is_vec else "")))
	print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/CR_%s_%s%s/fitResults_t1_%s_%s_%s combine_review/"
                %(options.chan, options.q, ("_vec" if is_vec else ""),options.chan, options.q, ("vec" if is_vec else "")))
	print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/CR_%s_%s%s/impacts_t0_%s_%s_%s.pdf combine_review/"
                %(options.chan, options.q, ("_vec" if is_vec else ""),options.chan, options.q, ("vec" if is_vec else "")))
	print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/CR_%s_%s%s/impacts_t1_%s_%s_%s.pdf combine_review/"
                %(options.chan, options.q, ("_vec" if is_vec else ""),options.chan, options.q, ("vec" if is_vec else "")))
