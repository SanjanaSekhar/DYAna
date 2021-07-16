import numpy as np
import ROOT
from hepdata_lib import *
from hepdata_lib.helpers import *

sub = Submission()

#main results
table_AFB = Table("Table 1")
table_AFB.keywords["observables"] = ["ASYM"]
table_AFB.keywords["phrases"] = ["Drell Yan", "Asymmetry Measurement", "Angular Coefficient"]
table_AFB.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-"]
table_AFB.description = """Measurement of the Drell-Yan forward-backward asymmetry as a function of dilepton mass."""
table_AFB.location = "Data from Table 1 located on page 15"



table_A0 = Table("Table 2")
table_A0.keywords["phrases"] = ["Drell Yan", "Angular Coefficient"]
table_A0.keywords["reactions"] = ["P P --> Z0/GAMMA* --> MU+ MU-", "P P --> Z0/GAMMA* --> E+ E-"]
table_A0.description = """Measurement of the Drell-Yan angular coefficient, A0, as a function of dilepton mass."""
table_A0.location = "Data from Table 2 located on page 15"

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

table_AFB.add_variable(mass_bins)
table_A0.add_variable(mass_bins)
channels = ["combined", "ee", "mumu"]
#channels = ["combined"]
file_dir = "../analyze/combine/AFB_fits/fit_results/"


for chan in channels:

    fname = file_dir + chan + "_results.npz"
    np_file = np.load(fname)

    AFB_vals = dict()
    AFB_vals["AFB_err_stat"] = np_file["AFB_err_stat"][1:]
    AFB_vals["AFB_err_sys"] = np_file["AFB_err_sys"][1:]
    AFB_vals["AFB_val"] = [0.6]*n_mass_bins
    #raw_vals["AFB_val"] = ["AFB_val"][1:]

    round_value_and_uncertainty(AFB_vals, val_key='AFB_val', unc_key='AFB_err_stat', sig_digits_unc = 2)
    round_value_and_uncertainty(AFB_vals, val_key='AFB_err_sys', unc_key='AFB_err_stat', sig_digits_unc = 2)

    A0_vals = dict()
    A0_vals["A0_err_stat"] = np_file["A0_err_stat"][1:]
    A0_vals["A0_err_sys"] = np_file["A0_err_sys"][1:]
    A0_vals["A0_val"] = [0.05]*n_mass_bins
    #raw_vals["A0_val"] = ["A0_val"][1:]

    round_value_and_uncertainty(A0_vals, val_key='A0_val', unc_key='A0_err_stat', sig_digits_unc = 2)
    round_value_and_uncertainty(A0_vals, val_key='A0_err_sys', unc_key='A0_err_stat', sig_digits_unc = 2)

    AFBs = Variable("AFB (%s)"  % chan,
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


    A0s = Variable("A0 (%s)"  % chan,
            is_independent = False,
            is_binned = False,
            units = "")
    A0s.values = [0.05]*n_mass_bins
    #A0.values = np_file["A0_val"]

    np_file = np.load(fname)

    A0_stat_uncs = Uncertainty("Statistical")
    A0_stat_uncs.values = A0_vals["A0_err_stat"]
    A0s.add_uncertainty(A0_stat_uncs)

    A0_sys_uncs = Uncertainty("Systematic")
    A0_sys_uncs.values = A0_vals["A0_err_sys"]
    A0s.add_uncertainty(A0_sys_uncs)

    table_A0.add_variable(A0s)



sub.add_table(table_AFB)
sub.add_table(table_A0)

#limit
table_zprime = Table("Figure 12")
table_zprime.keywords["reactions"] = ["P P --> ZPRIME --> MU+ MU-", "P P --> ZPRIME --> E+ E-"]
table_zprime.keywords["observables"] = ["SIG"]
table_zprime.description = "95% upper limits on the coupling parameter of a Zprime in the Sequential Standard Model as a function of mass"
table_zprime.location = "Data is from Figure 6 located on page 18"

zp_masses_var = Variable("Z' Mass",
        is_independent = True,
        is_binned = False,
        units = "GeV")

kl_limit_var = Variable("Limit on Kappa_L",
        is_independent = False,
        is_binned = False,
        units = "")



f_limit = ROOT.TFile.Open("../setlimits/limit.root", "READ")
#obs_limit_graph = f_limit->Get("obs_limit")
obs_limit_graph = f_limit.Get("exp_limit")
zp_masses_ = []
kl_limits_ = []

x=np.array([0.])
y=np.array([0.])
for i in range(obs_limit_graph.GetN()):
    obs_limit_graph.GetPoint(i,x,y)
    zp_masses_.append(x[0])
    kl_limits_.append(y[0])

n_masses = len(zp_masses_)
zp_masses_var.values = np.array(zp_masses_)
kl_limit_var.values = np.array(kl_limits_)

table_zprime.add_variable(zp_masses_var)
table_zprime.add_variable(kl_limit_var)


sub.add_table(table_zprime)
sub.create_files()

