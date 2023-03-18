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


year = -1
mLQ = 2000
for is_vec in [False]:
	for chan in ["ee"]:
		for q in ["u"]:

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

			sigma = 0.6 **0.5
			#extra_arg = "--symMCStats --sigma %f"%sigma	
			extra_arg = ""	
				
			print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s %s" % (comb_card, workspace, extra_arg))
			print_and_do("combine -M FitDiagnostics -d %s  --setParameters yLQ2=0.0,A0=0.03,A4=1.6 --freezeParameters A0,A4  --forceRecreateNLL -n _t0" %workspace)
			# WRONG --> print_and_do("python scripts/my_diffNuisances.py multidimfitTest.root --multidim --mLQ %i --prefit fitDiagnosticsTest.root -p yLQ2  -a fitDiagnostics_t0.root -g plots_t0.root > combine_review/fitResults_t0_%s_%s_%s"%(mLQ, chan,q,("vec" if is_vec else "")))
			print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  -a fitDiagnostics_t0.root -p yLQ2  -g plots_t0.root >> ./fitResults_t0_%s_%s_%s"%(chan,q,("vec" if is_vec else ""))) 
			#print_and_do("combine -M FitDiagnostics -d %s  --setParameters yLQ2=1.0,A0=0.05,A4=1.6  --forceRecreateNLL -n _t1" %workspace)
			# Increase the rMin value if (rMin * Nsig + Nbackground) < 0 for any channel
			#print_and_do("python scripts/my_diffNuisances.py multidimfitTest.root --multidim --mLQ %i --prefit fitDiagnosticsTest.root -p yLQ2  -a fitDiagnostics_t1.root -g plots_t1.root > combine_review/fitResults_t1_%s_%s_%s"%(mLQ, chan,q,("vec" if is_vec else "")))
			#print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics_t1.root -p yLQ2  -g plots_t1.root >> ./fitResults_t1_%s_%s_%s"%(chan,q,("vec" if is_vec else "")))
			#print_and_do("combineTool.py -M Impacts -d %s --setParameters yLQ2=0.0,A0=0.05,A4=1.6 -m 125 --doInitialFit --allPars -n t0"%workspace)
			#print_and_do("combineTool.py -M Impacts -d %s --setParameters yLQ2=1.0,A0=0.05,A4=1.6 -m 125 --doInitialFit --allPars -n t1"%(workspace))

			#print_and_do("combineTool.py -M Impacts -d %s -o impacts_t0_%s_%s_%s.json --setParameters yLQ2=0.0,A0=0.05,A4=1.6 --doFits -m 125 -n t0 "%(workspace,chan,q,("vec" if is_vec else "")))
			
			#print_and_do("combineTool.py -M Impacts -d %s -o impacts_t1_%s_%s_%s.json --setParameters yLQ2=1.0,A0=0.05,A4=1.6 --doFits -m 125 -n t1 "%(workspace,chan,q,("vec" if is_vec else "")))
			
			#print_and_do("combineTool.py -M Impacts -d %s  -m 125 -n t0 -o impacts_t0_%s_%s_%s.json"%(workspace,chan,q,("vec" if is_vec else "")))
			#print_and_do("combineTool.py -M Impacts -d %s  -m 125 -n t1 -o impacts_t1_%s_%s_%s.json"%(workspace,chan,q,("vec" if is_vec else "")))
			#print_and_do("plotImpacts.py -i  impacts_t0_%s_%s_%s.json -o  impacts_t0_%s_%s_%s"%(chan,q,("vec" if is_vec else ""),chan,q,("vec" if is_vec else "")))
			#print_and_do("plotImpacts.py -i  impacts_t1_%s_%s_%s.json -o  impacts_t1_%s_%s_%s"%(chan,q,("vec" if is_vec else ""),chan,q,("vec" if is_vec else "")))
			
			#print_and_do("cd cards")
			#print_and_do("ValidateDatacards.py %s"%comb_card)
