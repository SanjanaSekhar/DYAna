import numpy as np
import ROOT
from hepdata_lib import *
from hepdata_lib.helpers import *


def round_to_precision(vals, sigfigs):
    for i,v in enumerate(vals):
        v_new = round(v, sigfigs)
        vals[i] = v_new


def add_AFB(table_AFB, fname, chan, name): 
    np_file = np.load(fname)

    AFB_vals = dict()
    AFB_vals["AFB_err_stat"] = np_file["AFB_err_stat"][1:]
    AFB_vals["AFB_err_sys"] = np_file["AFB_err_sys"][1:]
    #AFB_vals["AFB_val"] = [0.6]*n_mass_bins
    AFB_vals["AFB_val"] = np_file["AFB_val"][1:]

    #This doesn't work b/c uncertainty changes place (from 0.009 to 0.011)
    #round_value_and_uncertainty(AFB_vals, val_key='AFB_val', unc_key='AFB_err_stat', sig_digits_unc = 2)
    #round_value_and_uncertainty(AFB_vals, val_key='AFB_err_sys', unc_key='AFB_err_stat', sig_digits_unc = 2)

    print("val", AFB_vals['AFB_val'])
    print("stat_unc", AFB_vals['AFB_err_stat'])
    print("sys unc", AFB_vals['AFB_err_sys'])

    round_to_precision(AFB_vals['AFB_val'], 3)
    round_to_precision(AFB_vals['AFB_err_stat'], 3)
    round_to_precision(AFB_vals['AFB_err_sys'], 3)

    print("val", AFB_vals['AFB_val'])
    print("stat_unc", AFB_vals['AFB_err_stat'])
    print("sys unc", AFB_vals['AFB_err_sys'])

    AFBs = Variable(name,
            is_independent = False,
            is_binned = False,
            units = "")
    AFBs.values = AFB_vals["AFB_val"]



    AFB_stat_uncs = Uncertainty("Statistical")
    AFB_stat_uncs.values = AFB_vals["AFB_err_stat"]
    AFBs.add_uncertainty(AFB_stat_uncs)

    AFB_sys_uncs = Uncertainty("Systematic")
    AFB_sys_uncs.values = AFB_vals["AFB_err_sys"]
    AFBs.add_uncertainty(AFB_sys_uncs)

    table_AFB.add_variable(AFBs) 

def add_A0(table_A0, fname, chan, name):
    np_file = np.load(fname)
    A0_vals = dict()
    A0_vals["A0_err_stat"] = np_file["A0_err_stat"][1:]
    A0_vals["A0_err_sys"] = np_file["A0_err_sys"][1:]
    
    A0_vals["A0_val"] = np_file["A0_val"][1:]

    #round_value_and_uncertainty(A0_vals, val_key='A0_val', unc_key='A0_err_stat', sig_digits_unc = 2)
    #round_value_and_uncertainty(A0_vals, val_key='A0_err_sys', unc_key='A0_err_stat', sig_digits_unc = 2)

    round_to_precision(A0_vals['A0_val'], 3)
    round_to_precision(A0_vals['A0_err_stat'], 3)
    round_to_precision(A0_vals['A0_err_sys'], 3)


    A0s = Variable(name,
            is_independent = False,
            is_binned = False,
            units = "")
    #A0s.values = [0.05]*n_mass_bins
    A0s.values = A0_vals["A0_val"]


    A0_stat_uncs = Uncertainty("Statistical")
    A0_stat_uncs.values = A0_vals["A0_err_stat"]
    A0s.add_uncertainty(A0_stat_uncs)

    A0_sys_uncs = Uncertainty("Systematic")
    A0_sys_uncs.values = A0_vals["A0_err_sys"]
    A0s.add_uncertainty(A0_sys_uncs)

    table_A0.add_variable(A0s)


sub = Submission()

#main results
table_AFB = Table("Table 2")
table_AFB.keywords["observables"] = ["ASYM"]
table_AFB.keywords["phrases"] = ["Drell Yan", "Asymmetry Measurement", "Angular Coefficient"]
table_AFB.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-"]
table_AFB.description = """Results for the measurement of $A_\mathrm{FB}$ from the maximum likelihood fit to data in different dilepton mass bins in the different channels as well as an inclusive measurement across all mass bins."""
table_AFB.location = "Data from Table 2 located on page 16"



table_A0 = Table("Table 3")
table_A0.keywords["phrases"] = ["Drell Yan", "Angular Coefficient"]
table_A0.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-"]
table_A0.description = """Results for the measurement of $A_0$ from the maximum likelihood fit to data in different dilepton mass bins in the different channels as well as inclusive measurement across all mass bins. To help in the interpretation of these results, we also list the average dilepton $p_{T}$ of the data events in each mass bin."""

table_A0.location = "Data from Table 3 located on page 17"


table_delta = Table("Table 4")
table_delta.keywords["phrases"] = ["Drell Yan", "Angular Coefficient"]
table_delta.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-"]
table_delta.description = """Results for the measurement of $\Delta A_\mathrm{FB}$ and $\Delta A_0$ between the muon and electron channels from the maximum likelihood fit to data in different mass bins as well as an inclusive measurement across all mass bins."""
table_delta.location = "Data from Table 4 located on page 18"



n_mass_bins = 7
mass_bins = Variable("Mass Bins",
        is_independent = True,
        is_binned = True,
        units = "GeV")

mass_bins.values = [(170., 200.), 
                    (200., 250.), 
                    (250., 320.),
                    (320., 510.),
                    (510., 700.), 
                    (700., 1000.), 
                    (1000., 13000.)]

n_mass_bins_ext = 8
mass_bins_ext = Variable("Mass Bins",
        is_independent = True,
        is_binned = True,
        units = "GeV")

mass_bins_ext.values = [(170., 200.), 
                    (200., 250.), 
                    (250., 320.),
                    (320., 510.),
                    (510., 700.), 
                    (700., 1000.), 
                    (1000., 13000.),
                    (170., 13000.),
                    ]


table_AFB.add_variable(mass_bins_ext)
table_A0.add_variable(mass_bins_ext)
channels = ["mumu", "ee", "combined"]
#channels = ["combined"]
file_dir = "../analyze/combine/AFB_fits/fit_results/"

#hardcoded values for avg pts
avg_pts = Variable("Avg Pt", is_independent = False, is_binned = False, units = "GeV")
avg_pts.values = [38., 43., 48., 55., 65., 73., 88., 44.] 
table_A0.add_variable(avg_pts)

for chan in channels:

    fname = file_dir + chan + "_results.npz"
    AFB_name = "AFB (%s)"  % chan
    A0_name = "A0 (%s)"  % chan
    add_AFB(table_AFB, fname, chan,AFB_name)
    add_A0(table_A0, fname, chan, A0_name)




sub.add_table(table_AFB)
sub.add_table(table_A0)



table_delta.add_variable(mass_bins_ext)

chan = "combined_diff"
fname = file_dir + chan + "_results.npz"
dAFB_name = "Delta AFB"
dA0_name = "Delta A0"
add_AFB(table_delta, fname, chan, dAFB_name)
add_A0(table_delta, fname, chan, dA0_name)
sub.add_table(table_delta)

#photon induced fraction
table_phot = Table("Table 5")
table_phot.keywords["phrases"] = ["Drell Yan", "Photon induced dilepton"]
table_phot.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-", "GAMMA GAMMA --> MU+ MU-", "GAMMA GAMMA --> E+ E-"]
table_phot.description = """ The fraction of photon-induced background as compared with the total amount of DY signal plus photon-induced events 
    ($N_{\gamma\gamma}/(N_{\gamma\gamma} + N_\mathrm{DY})$) in different dilepton mass bins.  These numbers are averaged across the different years and channels."""
table_phot.location = "Data from Table 5 located on page 19"
table_phot.add_variable(mass_bins)
phot_var = Variable(r"$\gamma\gamma \rightarrow \ell \ell$ fraction",
        is_independent = False,
        is_binned = False,
        units = "")
phot_var.values =  np.array([0.018, 0.021, 0.025, 0.028, 0.033, 0.037, 0.041])
table_phot.add_variable(phot_var)
sub.add_table(table_phot)


#limit
table_zprime = Table("Figure 5")
table_zprime.keywords["reactions"] = ["P P --> ZPRIME --> MU+ MU-", "P P --> ZPRIME --> E+ E-"]
table_zprime.keywords["observables"] = ["SIG"]
table_zprime.description = "Exclusion limits at 95% CL on the coupling K_L for a Z' in the sequential standard model as a function of the Z' mass."
table_zprime.location = "Data is from Figure 5 located on page 20"

zp_masses_var = Variable("Z' Mass",
        is_independent = True,
        is_binned = False,
        units = "GeV")

kl_limit_var = Variable("Observed Limit on Kappa_L",
        is_independent = False,
        is_binned = False,
        units = "")

exp_limit_var = Variable("Expected Limit on Kappa_L",
        is_independent = False,
        is_binned = False,
        units = "")


f_limit = ROOT.TFile.Open("../setlimits/limit_obs.root", "READ")
obs_limit_graph = f_limit.Get("obs_limit")
exp_limit_graph = f_limit.Get("exp_limit")
one_sig_graph = f_limit.Get("one_sigma_band")
two_sig_graph = f_limit.Get("two_sigma_band")
zp_masses_ = []
kl_limits_ = []
exp_limits_ = []

one_sig_ = []
two_sig_ = []

x=np.array([0.])
y=np.array([0.])
exp=np.array([0.])

for i in range(obs_limit_graph.GetN()):
    obs_limit_graph.GetPoint(i,x,y)
    exp_limit_graph.GetPoint(i,x,exp)
    zp_masses_.append(x[0])
    kl_limits_.append(y[0])
    exp_limits_.append(exp[0])

    one_sig_.append(one_sig_graph.GetErrorY(i))
    two_sig_.append(two_sig_graph.GetErrorY(i))

round_to_precision(exp_limits_, 2)
round_to_precision(one_sig_, 2)
round_to_precision(two_sig_, 2)

one_sig_uncs = Uncertainty("1 s.d.")
one_sig_uncs.values = one_sig_
exp_limit_var.add_uncertainty(one_sig_uncs)


two_sig_uncs = Uncertainty("2 s.d.")
two_sig_uncs.values = two_sig_
exp_limit_var.add_uncertainty(two_sig_uncs)

n_masses = len(zp_masses_)
zp_masses_var.values = np.array(zp_masses_)
exp_limit_var.values = np.array(exp_limits_)
kl_limit_var.values = np.array(kl_limits_)

table_zprime.add_variable(zp_masses_var)
table_zprime.add_variable(exp_limit_var)
table_zprime.add_variable(kl_limit_var)


sub.add_table(table_zprime)
sub.create_files()

