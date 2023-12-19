import ROOT
import subprocess
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
#from AFB_fits.scripts.LQ_utils import *

ext = "020923"

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--mLQ",  default=1000, type='int', help="mLQ")
(options, args) = parser.parse_args()

mLQ = options.mLQ
chan = "ee"

for mLQ in [2500]:
	for year in [2016,2017,2018]:

		fin = "combine/templates/LQm%i_merge_templates%i_%s.root"%(mLQ,year-2000,ext)
		f = ROOT.TFile.Open(fin, "UPDATE")
		f.cd("LQ")


		print("making signal temps: year is %i, chan is %s, mLQ is %i" % (year, chan, mLQ))
		
		hilow_procs = ["elHLTBARPT"]#, "elIDBARPT", "elRECOBARPT", "elHLTENDPT", "elIDENDPT", "elRECOENDPT"]

		procs = ['top']#,'db','fpl','fmn','alpha','tautau','gam','LQpure_u','LQint_u','LQpure_d','LQint_d','LQpure_u_vec','LQint_u_vec','LQpure_d_vec','LQint_d_vec']
		
		for proc in hilow_procs:
			for p in procs:
				hi_up = ROOT.gDirectory.Get("%s%s_%s_%sHIGH%sUp"%(chan,year-2000,p,proc,year-2000))
				hi_down = ROOT.gDirectory.Get("%s%s_%s_%sHIGH%sDown"%(chan,year-2000,p,proc,year-2000))
				lo_up = ROOT.gDirectory.Get("%s%s_%s_%sLOW%sUp"%(chan,year-2000,p,proc,year-2000))
				lo_down = ROOT.gDirectory.Get("%s%s_%s_%sLOW%sDown"%(chan,year-2000,p,proc,year-2000))
				h = ROOT.gDirectory.Get("%s%s_%s"%(chan,year-2000,p))
				
				for idx in range(1,h.GetNbinsX()+1):
					#keys = ROOT.gDirectory.GetListOfKeys().Clone()
					nom = h.GetBinContent(idx)
					nom_e = h.GetBinError(idx)
					new_up = hi_up.GetBinContent(idx) + lo_up.GetBinContent(idx) - nom
					new_down =  nom - hi_down.GetBinContent(idx) - lo_down.GetBinContent(idx)
					
					new_name_up = "%s%s_%s_%s%sUp"%(chan,year-2000,p,proc,year-2000)
					new_name_down = "%s%s_%s_%s%sDown"%(chan,year-2000,p,proc,year-2000)
					print("nom = %f, hi_up = %f, lo_up = %f, hi_down = %f, lo_down = %f\n"%(nom,hi_up.GetBinContent(idx),lo_up.GetBinContent(idx),hi_down.GetBinContent(idx),lo_down.GetBinContent(idx)))
					print("nom = %f, up = %f, down = %f\n"%(nom,new_up,new_down))
					print("New names: %s and %s"%(new_name_up,new_name_down))	
					#if((new_name+"Up" in keys) or (new_name+"Down" in keys) or (new_name2+"Up" in keys) or (new_name2+"Down" in keys)): continue
					#print("Creating %s"%new_name)
					#flag = 0
					#for s in keys:
					#	if s.GetName()+"Up"==new_name or s.GetName()+"Down"==new_name or s.GetName()+"Up"==new_name2 or s.GetName()+"Down"==new_name2: flag = 1
					#if flag == 1: continue
					# UP template
					h_clone = h.Clone(new_name_up)
					h_clone.SetBinContent(idx, new_up)
					
					h_clone.Write()
					
					# DOWN template
					h_clone = h.Clone(new_name_down)
					h_clone.SetBinContent(idx, new_down)
					
					h_clone.Write()
		del f