from utils import *
import numpy as np
import sys
from scipy.stats import norm
parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="all", type="string", help="What channels to run the fit over (all, combined, ee, or mumu)")
(options, args) = parser.parse_args()


frac_uncor = [0.42, 0.49, 0.52, 0.36, 0.29, 0.48, 0.2]


channel = 'combined_diff_'


A0_str = 'dA0'
Afb_str = 'dAfb'


def get_AFB_A0(filename):
    AFB_val,AFB_err,A0_val,A0_err = 0,0,0,0
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if A0_str in line:
                res = line.split()
                A0_val, A0_err = float(res[1]), float(res[3])
                #print('A0', A0_val, A0_err)
            if Afb_str in line:
                res = line.split()
                AFB_val, AFB_err = float(res[1]), float(res[3])
                #print('AFB', AFB_val, AFB_err)
    return AFB_val, A0_val, AFB_err, A0_err







n_bins = 7

mbin_str =  [ "150-170", "170-200", "200-250", "250-320", "320-510", "510-700", "700-1000", "$\\geq$ 1000"]


AFB_err_stat = np.array([0.]*n_bins)
AFB_err_sys = np.array([0.]*n_bins)
AFB_err_tot = np.array([0.]*n_bins)
A0_err_stat = np.array([0.]*n_bins)
A0_err_sys = np.array([0.]*n_bins)
A0_err_tot = np.array([0.]*n_bins)
AFB_val = np.array([0.]*n_bins)
A0_val = np.array([0.]*n_bins)

AFB_cov = np.array([[1.]*n_bins]*n_bins)







fit_dir = 'fit_results/'
ending = 'fit_results_mbin%i.txt'
output_file = fit_dir + channel + "_results"

for idx in range(1,n_bins+1):

    fit_name = fit_dir + channel + (ending % idx)
    stat_name = fit_dir + channel + 'nosys_' + (ending % idx)
    i = idx - 1
    AFB_fit, A0_fit, AFB_err_fit, A0_err_fit = get_AFB_A0(fit_name)
    _,_, AFB_err_stat[i], A0_err_stat[i] = get_AFB_A0(stat_name)

    AFB_err_sys[i] = (AFB_err_fit**2 - AFB_err_stat[i]**2 )**0.5
    A0_err_sys[i] = (A0_err_fit**2 - A0_err_stat[i]**2)**0.5

    if(np.isnan(AFB_err_sys[i])): AFB_err_sys[i] = 0
    if(np.isnan(A0_err_sys[i])): A0_err_sys[i] = 0

    AFB_err_tot[i] = AFB_err_fit
    A0_err_tot[i] = A0_err_fit

    AFB_val[i] = AFB_fit 
    A0_val[i] = A0_fit 
        

for i in range(n_bins):
    AFB_cov[i,i] = AFB_err_tot[i]**2
    correlated_err = AFB_err_sys[i] * (1. - frac_uncor[i-1])**0.5
    #correlated_err = AFB_err_sys[i]
    #correlated_err = 0.
    for j in range(n_bins):
        if(j != i): 
            AFB_cov[i,j] *= correlated_err
            AFB_cov[j,i] *= correlated_err

print(AFB_cov)


W = np.linalg.inv(AFB_cov)

print(W)

J = np.array([1.]* n_bins).T
var = 1./np.matmul(J.T, (np.matmul(W, J)))
unc = var**0.5
weighted_avg = var * (np.matmul(J.T, np.matmul(W, AFB_val)))

weights_1d_simple = [1./AFB_err_tot[i]**2 for i in range(n_bins)]
sum_weights = np.sum(weights_1d_simple)
unc_1d_simple = (1./sum_weights)**0.5

weights_1d = [W[i,i] for i in range(n_bins)]
weighted_avg_1d = np.average(AFB_val, weights = weights_1d)


print("\n \n %.4f +/- %.4f \n" % (weighted_avg, np.sqrt(var)))
print("Simple version %.3f \n" % (weighted_avg_1d))


def comp_test_stat(vals, mean, W):
    #vals = np.abs(vals)
    vals = np.expand_dims(vals, axis = -1)
    dim = len(vals.shape)
    vals_T = vals.swapaxes(dim-2, dim -1)

    L_sm = np.matmul(vals_T, np.matmul(W, vals))

    diff = vals - mean
    diff_T = diff.swapaxes(dim-2, dim -1)
    L_np = np.matmul(diff_T, np.matmul(W, diff))

    return (L_sm - L_np).reshape(-1)


test_stat_obs = comp_test_stat(AFB_val, weighted_avg, W)

print(test_stat_obs, abs(test_stat_obs)**0.5)

n_samples = 100000
means = np.array([0.] * n_bins)

samples = np.random.multivariate_normal(mean = means, cov = AFB_cov, size = n_samples)
sample_test_stats = comp_test_stat(samples, weighted_avg, W)
#sample_test_stats = comp_test_stat(np.abs(samples), np.abs(weighted_avg), W)

np_above = sample_test_stats > test_stat_obs
n_above = samples[np_above].shape[0]
np_p_val = float(n_above) / n_samples

print("P-val is %.4f" % np_p_val)
print("Sigma %.3f" % norm.ppf(np_p_val))
print("2.6 Sigma is pval of %.3f" % norm.cdf(2.6))




