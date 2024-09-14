import numpy as np 
import os, sys
mLQ = 2500
final = "\hline"
for vec in [False, True]:
    final = np.vstack((final,"Channel & $%s$ & Statistical unc. & Systematic unc. & Total unc.\\\\"%('\glq^2' if vec else '\ylq^2')))
    for chan in ["ee","mumu"]:
        for q in ["u","d"]:

            name = "%s_%s%s_unblinded"%(chan,q,("_vec" if vec else ""))
            results = np.loadtxt("postfit_plots/%s_LQ_m2500/results_%s_m2500.txt"%(name,name))
            name_stat = "%s_%s%s_statuncs_unblinded"%(chan,q,("_vec" if vec else ""))
            results_stat = np.loadtxt("postfit_plots/%s_LQ_m2500/results_%s_m2500.txt"%(name_stat,name_stat))

            print(name)
            print("Full unc on yLQ2: ",results[-1])
            print("Stat unc on yLQ2: ",results_stat[-1])
            sys_unc = [0,0,0]
            sys_unc[0] = results[2][0]
            sys_unc[1] = np.sqrt(results[2][1]**2 - results_stat[2][1]**2)
            sys_unc[2] = np.sqrt(results[2][2]**2 - results_stat[2][2]**2)
            print("Sys unc on yLQ2: ",sys_unc)

            '''
            \PVmd & 0.01 $\pm$ 0.05  & 1.61 $\pm$ 0.06 & 0.14 $\pm$ 0.05 (stat) $\pm$ 0.09 (syst) \\
            '''
            result_str = "\P%s%s%s & $%0.2f$ & $%0.2f/+%0.2f$ & $-%0.2f/+%0.2f$ & $%0.2f/+%0.2f$\\\\" %(('V' if vec else 'S'), chan[0], q,
                                results[2][0], results_stat[2][1], results_stat[2][2],
                                sys_unc[1], sys_unc[2],
                                results[2][1], results[2][2])
            final = np.vstack((final, result_str))


np.savetxt("unblinded_results.txt",final, delimiter=" ", fmt="%s")
