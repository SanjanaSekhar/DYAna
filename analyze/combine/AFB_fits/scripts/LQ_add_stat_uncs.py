import operator
import ROOT
from ROOT import *
from LQ_utils import *
from add_group_impact import *
gROOT.SetBatch(True)
import pandas as pd
import math

fake_data = True
no_sys = False
gen_level = False
no_LQ = False
year = -1
mLQ = 2500
ending = "040324_unblinding"
s = 3456
extra_params = " -s %i" % s

for is_vec in [False]:
    if is_vec: ending+="_vec"
    for chan in ["mumu","ee"]:
        for q in ["u","d"]:

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

	    
            workspace = "workspaces/%s_%s_sys_uncs_m%i.root" % (chan, q, mLQ)

            make_workspace(workspace, gen_level, chan, q, is_vec, no_LQ , no_sys, fake_data, mLQ, year,True, False)
            print_and_do("combine -M MultiDimFit -d %s --saveFitResult --saveWorkspace -n _base --robustFit 1  -s 3456 " % (workspace))

            print_and_do("cp higgsCombine_base.MultiDimFit.mH120.%i.root higgsCombine_nom.MultiDimFit.mH120.%i.root" % (s,s))
	    print_and_do("cp multidimfit_base.root multidimfit_nom.root")

            print_and_do("""combine -M MultiDimFit --freezeParameters allConstrainedNuisances -d higgsCombine_nom.MultiDimFit.mH120.%i.root --saveWorkspace  --saveFitResult --robustFit 1 -n _%s %s --snapshotName MultiDimFit """ %(s, 'statuncs', extra_params))
	    stat_unc = compute_sys("nom", "statuncs", s)

            df1 = pd.read_csv("%s_%s_m%s_sys_uncs_%s.txt"%(chan, q, mLQ, ending),delimiter="&",names=["Sys name","Contri","%% Contri"],dtype='string')
	    df1.drop(index=[0,len(df1)-1], inplace=True)
	    print(df1)
            df1["%% Contri"] = df1["%% Contri"].str.replace("\\","")
            df1["Contri"] = pd.to_numeric(df1["Contri"])
            df1["%% Contri"] = pd.to_numeric(df1["%% Contri"])

            sum_uncs2 = stat_unc**2 #variable to store total unc, not just sum of sys uncs
	    sum_sys_uncs2 = 0

            for i,val in enumerate(df1["Contri"]):
                print("adding to full uncs**2 -> ",df1.at[i+1,"Sys name"])
                sum_uncs2+=val**2
		sum_sys_uncs2+=val**2

            df1.loc[len(df1.index)+1] = ["Statistical Uncertainty", stat_unc, (stat_unc**2*100)/sum_uncs2]
            for i in range(1,len(df1.index)):
                print("computing \% contri -> ",df1.at[i,"Sys name"])
                df1.at[i,"%% Contri"] = (df1.at[i,"Contri"]**2*100)/(sum_uncs2)
	    
	    #print(df1)
	    print("Background cross sections: ", df1.loc[df1["Sys name"].str.contains("Section")])
	    df_xsec = df1.loc[df1["Sys name"].str.contains("Section")]
            print("MisID: ", df1.loc[df1["Sys name"].str.contains("MisID")])
	    df_misid = df1.loc[df1["Sys name"].str.contains("MisID")]
	    print("Trigger: ", df1.loc[df1["Sys name"].str.contains("Trigger")])
	    df_trigger = df1.loc[df1["Sys name"].str.contains("Trigger")]
	   
	    df1.loc[len(df1.index)+2] = ["Background Cross Sections",df_xsec["Contri"].sum(),df_xsec["%% Contri"].sum()]
	    df1.loc[len(df1.index)+3] = ["Trigger & Prefire",df_trigger["Contri"].sum(),df_trigger["%% Contri"].sum()]
	    df1.loc[len(df1.index)+4] = ["MisID shape and normalization",df_misid["Contri"].sum(),df_misid["%% Contri"].sum()]
	    df1.loc[len(df1.index)+5] = ["Total systematic uncertainty", math.sqrt(sum_sys_uncs2),sum_sys_uncs2*100/(sum_uncs2)] 
	    print(df1)
	    df1.to_csv("%s_%s_m%s_sys_uncs_%s_statuncs.txt"%(chan, q, mLQ, ending),sep=' ',index=False)

            
