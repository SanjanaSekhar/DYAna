from utils import *
parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
(options, args) = parser.parse_args()


AFB_shifts = [0.0, 0.016, 0.009, 0.005, 0.002, 0.002, -0.001, -0.001]
AFB_shifts_unc = [0.0, 0.003, 0.002, 0.002, 0.001, 0.002, -0.001, -0.001]


def get_AFB_A0(filename):
    AFB_val,AFB_err,A0_val,A0_err = 0,0,0,0
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if 'A0' in line:
                res = line.split()
                A0_val, A0_err = float(res[1]), float(res[3])
                #print('A0', A0_val, A0_err)
            if 'Afb' in line:
                res = line.split()
                AFB_val, AFB_err = float(res[1]), float(res[3])
                #print('AFB', AFB_val, AFB_err)
    return AFB_val, A0_val, AFB_err, A0_err


n_bins = 8

chan = options.chan + "_"
fit_dir = 'fit_results/'
ending = 'fit_results_mbin%i.txt'
AFB_val = [0]*n_bins
A0_val = [0]*n_bins
AFB_err_stat = [0]*n_bins
AFB_err_sys = [0]*n_bins
AFB_err_tot = [0]*n_bins
A0_err_stat = [0]*n_bins
A0_err_sys = [0]*n_bins
A0_err_tot = [0]*n_bins

for i in range(1,n_bins):

    fit_name = fit_dir + chan + (ending % i)
    stat_name = fit_dir + chan + 'nosys_' + (ending % i)
    AFB_fit, A0_fit, AFB_err_fit, A0_err_fit = get_AFB_A0(fit_name)
    _,_, AFB_err_stat[i], A0_err_stat[i] = get_AFB_A0(stat_name)

    AFB_err_sys[i] = (AFB_err_fit**2 - AFB_err_stat[i]**2 + AFB_shifts_unc[i]**2)**0.5
    A0_err_sys[i] = (A0_err_fit**2 - A0_err_stat[i]**2)**0.5

    AFB_err_tot[i] = (AFB_err_fit**2 + AFB_shifts_unc[i]**2)**0.5
    A0_err_tot[i] = A0_err_fit

    AFB_val[i] = AFB_fit + AFB_shifts[i]
    A0_val[i] = A0_fit

print("AFB:")
for i in range(1,n_bins):
    #print("%.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(AFB_val[i], AFB_err_stat[i], AFB_err_sys[i]))
    print("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(AFB_err_stat[i], AFB_err_sys[i]))

print("A0:")
for i in range(1,n_bins):
    #print("%.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(A0_val[i], A0_err_stat[i], A0_err_sys[i]))
    print("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(A0_err_stat[i], A0_err_sys[i]))
