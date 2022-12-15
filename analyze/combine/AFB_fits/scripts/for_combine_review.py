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

for is_vec in [False,True]:
	for chan in ["ee","mumu"]:
    	for q in ["u","d","s"]:

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

		    
		    else: extra_arg = ""
		    print_and_do("text2workspace.py %s -P LQ_Analysis.DYAna.LQ_my_model:lq_ylq_sq -o %s " % (comb_card, workspace))
		    print_and_do("combine -M FitDiagnostics -d %s -t -1 --expectSignal 0 --rMin -10 --forceRecreateNLL -n _t0" %workspace)
		    print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  -a fitDiagnostics_t0.root -g plots_t0.root >> ./fitResults_t0_%s_%s_%s"%(chan,q,("vec" if is_vec else ""))

		    print_and_do("combine -M FitDiagnostics -d %s -t -1 --expectSignal 1  --forceRecreateNLL -n _t1" %workspace)
			# Increase the rMin value if (rMin * Nsig + Nbackground) < 0 for any channel
			print_and_do("python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  -a fitDiagnostics_t1.root -g plots_t1.root >> ./fitResults_t1_%s_%s_%s"%(chan,q,("vec" if is_vec else ""))


