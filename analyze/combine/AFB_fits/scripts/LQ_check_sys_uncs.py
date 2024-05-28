import operator
import ROOT
from ROOT import *
from LQ_utils import *
from add_group_impact import *
import numpy as np
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--expected",  default=False, action="store_true", help="Compute expected impacts based on toys with AFB=0.6 A0=0.05")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject if expected")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject if expected")
parser.add_option("-o", "--odir", default="sys_uncs/", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--diff", default=False, action="store_true", help="Diff")
parser.add_option("--mLQ",  default=2000, type='int', help="mLQ")
parser.add_option("--nToys",  default=5, type='int', help="no. of toys for expected uncs")
parser.add_option("--vec",  default=False, help="is vec?")
parser.add_option("--chan",  default="ee", help="channel ee or mumu ")
parser.add_option("--q",  default="u", help=" channel u,d")
parser.add_option("--hadd",  default=False, help="hadd")
parser.add_option("--ending", default="041123", help="date")
parser.add_option("-s","--seed",  default=3456, type='int', help="random seed")
(options, args) = parser.parse_args()

chan = options.chan
q = options.q
fake_data = True
no_sys = False
gen_level = False
no_LQ = False
year = -1
is_vec = options.vec
extra_params = ""
ending = options.ending
s = options.seed
if is_vec: ending+="_vec"
extra_params += " -s %i" % s
if options.expected: ending += "_expected"
if options.hadd:

	
	for chan in ["mumu","ee"]:
		for q in ["u","d"]:
			plt.figure(figsize=(9,10))
			for options.mLQ in [1000,1500,2000,2500,3000,3500]:
				for i in range(1,options.nToys+1):
					#print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/sys_%s_%s%s_%s_exp_toy%i/%s_%s_m%s_sys_uncs_%s_toy1.txt sys_uncs/%s_%s_m%s_sys_uncs_%s_toy%i.txt"%(chan, q, ("_vec" if is_vec else ""),options.mLQ,(i+5), chan, q, options.mLQ, ending,  chan, q, options.mLQ, ending, (i+5) ))
					#print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/sys_%s_%s%s_%s_exp/%s_%s_m%s_sys_uncs_%s_toy%i.txt sys_uncs/"%(chan, q, ("_vec" if is_vec else ""),options.mLQ, chan, q, options.mLQ, ending, i))
					
					df1 = pd.read_csv("sys_uncs/%s_%s_m%s_sys_uncs_%s_toy%i.txt"%(chan, q, options.mLQ, ending, i),delimiter="&",names=["Sys name","Contri","%% Contri"],dtype='string')
					df1.drop(index=[0,len(df1)-1], inplace=True)
					df1["%% Contri"] = df1["%% Contri"].str.replace("\\","")
					#df1["Contri"] = pd.to_numeric(df1["Contri"])
					df1.drop(columns="Contri", inplace=True)
					df1["%% Contri"] = pd.to_numeric(df1["%% Contri"])
					if i==1: df_all = df1.copy()
					else: 
						#df_all = pd.concat([df_all,df1])
						df_all = pd.merge(df_all, df1, on="Sys name")
					
				df_all["Mean"] = df_all.mean(axis=1)
				df_all["Std"] = df_all.std(axis=1)
				df_all["Sys name"] = df_all["Sys name"].str.replace("\PGg","\gamma",regex=False)
				print("Average uncs for all mLQ: ",chan, q)
				#print(df_all[['Sys name','Mean']])						
				df_all.to_csv("sys_uncs/%s_%s_m%s_sys_uncs_%s_alltoys.txt"%(chan, q, options.mLQ, ending),sep=' ',index=False)
				plt.errorbar(df_all['Sys name'],df_all['Mean'],yerr=df_all["Std"], linestyle='none', marker='o',label="mLQ=%i GeV"%(options.mLQ))
			plt.legend(fontsize='large')
			plt.title("%% contribution of systematics to %s-%s%s channel"%(chan,q,("-vec" if is_vec else "")))
			plt.xlabel("Systematic",fontsize=10)
			plt.ylabel("%% contribution averaged over 15 toys")
			plt.xticks(rotation=90,fontsize=7)
			plt.tight_layout()
			plt.savefig("sys_uncs/%s_%s_sys_uncs_%s_allmLQ.png"%(chan, q, ending))
			plt.close()
	
	for options.mLQ in [1000,1500,2000,2500,3000,3500]:
		for chan in ["mumu","ee"]:
			plt.figure(figsize=(9,10))
			for q in ["u","d"]:
				for i in range(1,options.nToys+1):
					#print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/sys_%s_%s%s_%s_exp_toy%i/%s_%s_m%s_sys_uncs_%s_toy1.txt sys_uncs/%s_%s_m%s_sys_uncs_%s_toy%i.txt"%(chan, q, ("_vec" if is_vec else ""),options.mLQ,(i+5), chan, q, options.mLQ, ending,  chan, q, options.mLQ, ending, (i+5) ))
					#print_and_do("xrdcp -f root://cmseos.fnal.gov//store/user/ssekhar/Condor_outputs/sys_%s_%s%s_%s_exp/%s_%s_m%s_sys_uncs_%s_toy%i.txt sys_uncs/"%(chan, q, ("_vec" if is_vec else ""),options.mLQ, chan, q, options.mLQ, ending, i))
					
					df1 = pd.read_csv("sys_uncs/%s_%s_m%s_sys_uncs_%s_toy%i.txt"%(chan, q, options.mLQ, ending, i),delimiter="&",names=["Sys name","Contri","%% Contri"],dtype='string')
					df1.drop(index=[0,len(df1)-1], inplace=True)
					df1["%% Contri"] = df1["%% Contri"].str.replace("\\","")
					#df1["Contri"] = pd.to_numeric(df1["Contri"])
					df1.drop(columns="Contri", inplace=True)
					df1["%% Contri"] = pd.to_numeric(df1["%% Contri"])
					if i==1: df_all = df1.copy()
					else: 
						#df_all = pd.concat([df_all,df1])
						df_all = pd.merge(df_all, df1, on="Sys name")
					
				df_all["Mean"] = df_all.mean(axis=1)
				df_all["Std"] = df_all.std(axis=1)
				df_all["Sys name"] = df_all["Sys name"].str.replace("\PGg","\gamma",regex=False)
				print("Average uncs for all channels: ", chan, options.mLQ)
				#print(df_all[['Sys name','Mean','Std']])						
				plt.errorbar(df_all['Sys name'],df_all['Mean'], yerr=df_all["Std"],linestyle='none', marker='o', label="%s-%s%s"%(chan,q,("-vec" if is_vec else "")))
			plt.legend(fontsize='large')
			plt.title("%% contribution of systematics to %s%s channels, mLQ=%i GeV"%(chan,("-vec" if is_vec else ""),options.mLQ))
			plt.xlabel("Systematic")
			plt.ylabel("%% contribution averaged over 15 toys")
			plt.xticks(rotation=90,fontsize=7)
			plt.tight_layout()
			plt.savefig("sys_uncs/%s_sys_uncs_%s_m%i.png"%(chan, ending,options.mLQ))
			plt.close()	
else:
	if chan == "ee":
		#individual_pars, group_pars = ["nlo_sys"],[]		
		individual_pars = ["nlo_sys", "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "elFakesYR",  "Pu", "prefireYR"]
		
		group_pars =["RFscalesYRC", "elScalesYR", "elIDs", "elHLTsYR", "emucostrwsYRC",  "pdfs", "lumisYR", "elRECOs", "elfakesrwsYR", "autoMCStats,MCStatBin16,MCStatBin17,MCStatBin18"] 
		
	else:
		individual_pars = ["nlo_sys", "dy_xsec", "db_xsec",  "top_xsec", "gam_xsec",  "muFakesYR", "Pu", "muPrefYRC",  "muRCYR", ]
		group_pars =["RFscalesYRC","muIDsYR", "emucostrwsYRC", "pdfs", "lumisYR", "muHLTsYR", "mufakesrwsYR", "autoMCStats,MCStatBin16,MCStatBin17,MCStatBin18"] 

	sys_name_conv = dict()
	sys_name_conv['nlo_sys'] = "LQ LO reweighting"
	sys_name_conv['dy_xsec'] = "DY Cross Section"
	sys_name_conv['db_xsec'] = "Diboson Cross Section"
	sys_name_conv['top_xsec'] = "$\\ttbar$ Cross Section"
	sys_name_conv['gam_xsec'] = "$\\PGg\\PGg$ Cross Section"
	sys_name_conv['RFscalesYRC'] = "$\\alpha_s$ + Renormalization/Factorization Scales"
	sys_name_conv['emucostrwsYRC'] = "$e\\mu$ Shape Corrections"
	#sys_name_conv['ptrwsYRC'] = "DY $p_{T}$ Correction"
	sys_name_conv['pdfs'] = "PDFs"
	sys_name_conv['lumisYR'] = "Luminosity"
	sys_name_conv['autoMCStats,MCStatBin16,MCStatBin17,MCStatBin18'] = "MC and MisID Backgrounds Statistical Uncertainty"
	sys_name_conv['Pu'] = "Pileup"

	if chan == "ee": 
		sys_name_conv['elFakesYR'] = "Electron MisID Normalization"
		sys_name_conv['elScalesYR'] = "Electron Momentum Scale"
		sys_name_conv['elHLTsYR'] = "Electron Trigger"
		sys_name_conv['elIDs'] = "Electron Identification/Isolation"
		sys_name_conv['elRECOs'] = "Electron Reconstruction"
		sys_name_conv['elfakesrwsYR'] = "Electron MisID Shape"
		sys_name_conv['prefireYR'] = "Electron Trigger Prefire Correction"

	if chan== "mumu": 
		sys_name_conv['muFakesYR'] = "Muon MisID Normalization"
		sys_name_conv['muRCYR'] = "Muon Momentum Scale"
		sys_name_conv['muIDsYR'] = "Muon Identification/Isolation"
		sys_name_conv['muHLTsYR'] = "Muon Trigger"
		sys_name_conv['muPrefYRC'] = "Muon Trigger Prefire Correction"
		sys_name_conv['mufakesrwsYR'] = "Muon MisID Shape"
		sys_name_conv['muPrefYRC'] = "Muon Trigger Prefire Correction"


	#sys_name_conv['METJECYR'] = "MET Uncertainties"
	#sys_name_conv['BTAGSYR'] = "b-tagging Uncertainty"





	def par_to_freezestr(par):
		if("YRC" in par):
			par16 = par.replace("YRC", "16")
			par1718 = par.replace("YRC", "1718")
			return par16 + "," + par1718
		elif("YR" in par):

			if  "prefire" not in par:
				par16 = par.replace("YR", "16")
				par17 = par.replace("YR", "17")
				par18 = par.replace("YR", "18")
				return par16 + "," + par17 + "," + par18
			else: #prefire only sys with no 18
				if chan=="ee":   return "prefire16,prefire17"
		else:
			return par

			





	workspace = "workspaces/%s_%s_sys_uncs_m%i.root" % (chan, q, options.mLQ)


	make_workspace(workspace, gen_level, chan, q, is_vec, no_LQ , no_sys, fake_data, options.mLQ, year,True, False)
	print_and_do("combine -M MultiDimFit -d %s --saveFitResult --saveWorkspace -n _base --robustFit 1  -s 3456 " % (workspace))
	#print_and_do("combine -M FitDiagnostics -d %s  --saveWorkspace -n _base --robustFit 1  %s"     % (workspace, extra_params))

	if(options.expected):
		print("Will inject A4=1.61 yLQ2=0 for all toys ")

		for i in range(1,options.nToys+1):
			s += 1
			print_and_do(("combine -M GenerateOnly -d higgsCombine_base.MultiDimFit.mH120.3456.root --snapshotName MultiDimFit --toysFrequentist"
			" --bypassFrequentistFit --saveToys -t 1 -s %i  --setParameters A4=1.61,yLQ2=0.")
					% (s))

			print_and_do(("combine -M MultiDimFit -d %s -n _nom --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root " +
			"--toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL -s %i") %(workspace, s, s))

			extra_params = " --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist -t 1 -s %i" % (s,s)

			d = dict()
			n = 0
			for indi_par in individual_pars:
				n+=1
				freeze_str = par_to_freezestr(indi_par)
				#print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit --robustFit 1 -n _%s %s""" % (freeze_str,s, indi_par, extra_params))
					#  % (freeze_str,s, indi_par, extra_params))
				#print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.FitDiagnostics.mH120.%i.root -w w  --robustFit 1 -n _%s %s""" % (freeze_str,s, indi_par, extra_params))
				print_and_do("""combine -M MultiDimFit --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ % (freeze_str,s, indi_par, extra_params))
				sys_unc = compute_sys("nom", indi_par, s)
				#sys_unc = compute_sys("nom", indi_par, s)
				d[indi_par] = sys_unc
				#if(n>2): break

			for group_par in group_pars:
				n+=1
				freeze_str = par_to_freezestr(group_par)
				#print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit --robustFit 1 -n _%s %s""" % 
				#       (freeze_str, s, group_par, extra_params))
				#print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups %s -d higgsCombine_nom.FitDiagnostics.mH120.%i.root -w w  --robustFit 1 -n _%s %s""" %  (freeze_str, s, group_par, extra_params))
				if group_par == 'autoMCStats,MCStatBin16,MCStatBin17,MCStatBin18': 
					print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ %(freeze_str, s, 'mcstats', extra_params))
					sys_unc = compute_sys("nom", "mcstats", s)
				else:
					print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ %(freeze_str, s, group_par, extra_params))
					sys_unc = compute_sys("nom", group_par, s)
				#sys_unc = compute_sys("nom", indi_par, s)
				d[group_par] =sys_unc
				#if(n>4): break

			print(d)
			sum_uncs2 = 0
			os.system("mkdir %s \n" % options.odir)
			with open("%s/%s_%s_m%s_sys_uncs_%s_toy%i.txt" % (options.odir, chan, q, options.mLQ,ending, i), 'w') as f_out:
				sorted_d = sorted(d.items(), key=operator.itemgetter(1))
				f_out.write("Systematic uncertainties for %s-%s%s channel, mLQ = %s, toy = %i\n"%(chan, q, ("-vec" if is_vec else ""), options.mLQ, i))
				
				for sys_name, val in sorted_d[::-1]:
					# if (sys_name in sys_name_conv.keys()):
					# 	out_name = sys_name_conv[sys_name]
					# else:
					# 	out_name = sys_name
					sum_uncs2 += (val*val)

				for sys_name, val in sorted_d[::-1]:
					if (sys_name in sys_name_conv.keys()):
						out_name = sys_name_conv[sys_name]
					else:
						out_name = sys_name
					f_out.write("%s & %.5f & %.2f  \\\\ \n" % (out_name, val, ((val*val)*100)/sum_uncs2))
				
				f_out.write("Total Uncertainty & %.2f  \\\\ \n" % (np.sqrt(sum_uncs2)))

			#print_and_do("rm higgsCombine* fitDiagnostics* multidimfit*")

	else:
		print_and_do("cp higgsCombine_base.MultiDimFit.mH120.%i.root higgsCombine_nom.MultiDimFit.mH120.%i.root" % (s,s))
		#print_and_do("cp higgsCombine_base.FitDiagnostics.mH120.%i.root higgsCombine_nom.FitDiagnostics.mH120.%i.root" % (s,s))
		print_and_do("cp multidimfit_base.root multidimfit_nom.root")

		d = dict()
		n = 0
		for indi_par in individual_pars:
			n+=1
			freeze_str = par_to_freezestr(indi_par)
			#print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit --robustFit 1 -n _%s %s""" % (freeze_str,s, indi_par, extra_params))
				#  % (freeze_str,s, indi_par, extra_params))
			#print_and_do("""combine -M FitDiagnostics --freezeParameters %s -d higgsCombine_nom.FitDiagnostics.mH120.%i.root -w w  --robustFit 1 -n _%s %s""" % (freeze_str,s, indi_par, extra_params))
			print_and_do("""combine -M MultiDimFit --freezeParameters %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ % (freeze_str,s, indi_par, extra_params))
			sys_unc = compute_sys("nom", indi_par, s)
			#sys_unc = compute_sys("nom", indi_par, s)
			d[indi_par] = sys_unc
			#if(n>2): break

		for group_par in group_pars:
			n+=1
			freeze_str = par_to_freezestr(group_par)
			#print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit --robustFit 1 -n _%s %s""" % 
			#       (freeze_str, s, group_par, extra_params))
			#print_and_do("""combine -M FitDiagnostics --freezeNuisanceGroups %s -d higgsCombine_nom.FitDiagnostics.mH120.%i.root -w w  --robustFit 1 -n _%s %s""" %  (freeze_str, s, group_par, extra_params))
			if group_par == 'autoMCStats,MCStatBin16,MCStatBin17,MCStatBin18': 
				print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ %(freeze_str, s, 'mcstats', extra_params))
				sys_unc = compute_sys("nom", "mcstats", s)
			else:
				print_and_do("""combine -M MultiDimFit --freezeNuisanceGroups %s -d higgsCombine_nom.MultiDimFit.mH120.%i.root  --saveFitResult --robustFit 1 -n _%s %s """ %(freeze_str, s, group_par, extra_params))
				sys_unc = compute_sys("nom", group_par, s)
			#sys_unc = compute_sys("nom", indi_par, s)
			d[group_par] =sys_unc
			#if(n>4): break

		print(d)
		sum_uncs2 = 0
		os.system("mkdir %s \n" % options.odir)
		with open("%s/%s_%s_m%s_sys_uncs_%s.txt" % (options.odir, chan, q, options.mLQ,ending), 'w') as f_out:
			sorted_d = sorted(d.items(), key=operator.itemgetter(1))
			f_out.write("Systematic uncertainties for bin %i\n" % options.mbin)
			
			for sys_name, val in sorted_d[::-1]:
				# if (sys_name in sys_name_conv.keys()):
				# 	out_name = sys_name_conv[sys_name]
				# else:
				# 	out_name = sys_name
				sum_uncs2 += (val*val)

			for sys_name, val in sorted_d[::-1]:
				if (sys_name in sys_name_conv.keys()):
					out_name = sys_name_conv[sys_name]
				else:
					out_name = sys_name
				f_out.write("%s & %.5f & %.2f  \\\\ \n" % (out_name, val, ((val*val)*100)/sum_uncs2))
			
			f_out.write("Total Uncertainty & %.2f  \\\\ \n" % (np.sqrt(sum_uncs2)))

			#print_and_do("rm higgsCombine* fitDiagnostics* multidimfit*")
	#print_and_do("""combine -M FitDiagnostics -d higgsCombine_nom.MultiDimFit.mH120.%i.root -w w --snapshotName MultiDimFit -n _nom1 %s """ % (s,extra_params))
	#print_and_do("""combine -M FitDiagnostics -d higgsCombine_nom.FitDiagnostics.mH120.%i.root -w w --snapshotName FitDiagnostics -n _nom1 %s """ % (s,extra_params))
