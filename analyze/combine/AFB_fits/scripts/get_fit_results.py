from utils import *
import numpy as np
import sys
parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="all", type="string", help="What channels to run the fit over (all, combined, ee, or mumu)")
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

mbin_str =  [ "150-170", "170-200", "200-250", "250-320", "320-510", "510-700", "700-1000", "$\\geq$ 1000"]


AFB_err_stat = np.array([[0.]*n_bins]*3)
AFB_err_sys = np.array([[0.]*n_bins]*3)
AFB_err_tot = np.array([[0.]*n_bins]*3)
A0_err_stat = np.array([[0.]*n_bins]*3)
A0_err_sys = np.array([[0.]*n_bins]*3)
A0_err_tot = np.array([[0.]*n_bins]*3)
AFB_val = np.array([[0.]*n_bins]*3)
A0_val = np.array([[0.]*n_bins]*3)


if(options.chan == "all"):
    chans = ["ee_", "mumu_", "combined_"]
else:
    chans = [options.chan + "_"]




for c_idx, chan in enumerate(chans):

    fit_dir = 'fit_results/'
    ending = 'fit_results_mbin%i.txt'
    output_file = fit_dir + options.chan + "_results"

    for i in range(1,n_bins):

        fit_name = fit_dir + chan + (ending % i)
        stat_name = fit_dir + chan + 'nosys_' + (ending % i)
        AFB_fit, A0_fit, AFB_err_fit, A0_err_fit = get_AFB_A0(fit_name)
        _,_, AFB_err_stat[c_idx][i], A0_err_stat[c_idx][i] = get_AFB_A0(stat_name)

        AFB_err_sys[c_idx][i] = (AFB_err_fit**2 - AFB_err_stat[c_idx][i]**2 + AFB_shifts_unc[i]**2)**0.5
        A0_err_sys[c_idx][i] = (A0_err_fit**2 - A0_err_stat[c_idx][i]**2)**0.5

        if(np.isnan(AFB_err_sys[c_idx][i])): AFB_err_sys[c_idx][i] = 0
        if(np.isnan(A0_err_sys[c_idx][i])): A0_err_sys[c_idx][i] = 0

        AFB_err_tot[c_idx][i] = (AFB_err_fit**2 + AFB_shifts_unc[i]**2)**0.5
        A0_err_tot[c_idx][i] = A0_err_fit

        AFB_val[c_idx][i] = AFB_fit + AFB_shifts[i]
        A0_val[c_idx][i] = A0_fit

        

    np.savez(output_file, AFB_val=AFB_val[c_idx], AFB_err_stat=AFB_err_stat[c_idx], AFB_err_sys=AFB_err_sys[c_idx], 
            A0_val = A0_val[c_idx], A0_err_stat = A0_err_stat[c_idx], A0_err_sys = A0_err_sys[c_idx])

print("Channels are: ")
for chan in chans:
    print(chan)

print ("total uncs")
for j,chan in enumerate(chans):
    print(chan)
    print("AFB")
    for i in range(1, n_bins):
        print( "%.3f $\\pm$ %.3f" %(AFB_val[j][i], AFB_err_tot[j][i]))

    print("A0")
    for i in range(1, n_bins):
        print( "%.3f $\\pm$ %.3f" %(A0_val[j][i], A0_err_tot[j][i]))


print("\nAFB:")
for i in range(1,n_bins):
    sys.stdout.write( "%s    "  % mbin_str[i], )
    for j in range(len(chans)):
        sys.stdout.write( "& %.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)  " %(AFB_val[j][i], AFB_err_stat[j][i], AFB_err_sys[j][i]), )
    sys.stdout.write(" \\\\ \n")
    #sys.stdout.write("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(AFB_err_stat[i], AFB_err_sys[i]))

sys.stdout.write("\nA0:")
for i in range(1,n_bins):
    sys.stdout.write( "%s    "  % mbin_str[i], )
    for j in range(len(chans)):
        sys.stdout.write("& %.3f $\\pm$ %.3f (stat) $\\pm$ %.3f (syst) " %(A0_val[j][i], A0_err_stat[j][i], A0_err_sys[j][i]), )
    sys.stdout.write(" \\\\ \n")
    #sys.stdout.write("0.XXX $\\pm$ %.3f (stat) $\\pm$ %.3f (syst)" %(A0_err_stat[i], A0_err_sys[i]))


