import numpy as np
from hepdata_lib import Submission, Table, Variable, Uncertainty

sub = Submission()

table_AFB = Table("Table 1")
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
                    (1000., 14000.)]

comb_AFBs = Variable("AFB (combined)",
        is_independent = False,
        is_binned = False,
        units = "")
comb_AFBs.values = [0.6]*n_mass_bins

comb_AFB_stat_uncs = Uncertainty("Statistical")
comb_AFB_stat_uncs.values = [0.02]*n_mass_bins
comb_AFBs.add_uncertainty(comb_AFB_stat_uncs)

comb_AFB_sys_uncs = Uncertainty("Systematic")
comb_AFB_sys_uncs.values = [0.02]*n_mass_bins
comb_AFBs.add_uncertainty(comb_AFB_sys_uncs)

table_AFB.add_variable(mass_bins)
table_AFB.add_variable(comb_AFBs)

sub.add_table(table)
sub.create_files()

