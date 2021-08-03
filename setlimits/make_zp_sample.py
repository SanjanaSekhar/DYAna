import sys, commands, os
import numpy as np

#Mzp=10000
KL_start = 0.65
KL_stop = 2.05
KL_step = 0.05
eps = 0.001


KLs = np.linspace(KL_start,KL_stop, num = round((KL_stop - KL_start + eps ) / KL_step) + 1)
print(KLs)

for Mzp in range(4100, 4600, 100):
    #for KL_idx in range(3):
        #KL = 0.05 * (KL_idx) + KL_start
    for KL in KLs:
        

        print("Generating events for Mzp = %i KL = %.2f " % (Mzp, KL))

        cur_dir = os.getcwd()
        os.chdir("MG5_aMC_v2_6_2")
        os.system("rm -r Zp_to_mumu_temp")
        os.system("cp scripts/Zp_to_mumu_LO_TEMPLATE.txt scripts/Zp_to_mumu_temp.txt")
        os.system("""sed -i "s/MASS_ZP/%i/g" scripts/Zp_to_mumu_temp.txt""" % Mzp)
        os.system("""sed -i "s/KAPPA_L/%.2f/g" scripts/Zp_to_mumu_temp.txt""" % KL)
        os.system("./bin/mg5_aMC scripts/Zp_to_mumu_temp.txt")
        os.system("gunzip -d Zp_to_mumu_temp/Events/run*/*")
        os.chdir(cur_dir)
        os.system("ls")
        os.system(""" echo ".x record_AFBs.C(%i,%.2f)" > temp """ % (Mzp, KL))
        os.system(""" echo ".q" >> temp """ )
        os.system("root -l < temp")
        os.system("rm temp")
