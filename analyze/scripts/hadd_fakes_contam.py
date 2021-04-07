

import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

years = [2018]
#prefixes = ["MuMu", "ElEl"]
prefixes = ["EMu"]
#years = [2018]
#prefixes = ["ElEl"]
#labels = ["wt", "dy", "ttbar", "diboson", "phot_ind"]
labels = ["wt", "dy", "ttbar", "diboson"]
ending = "_april5"

for year in years:
    for prefix in prefixes:
        f_names = []
        for label in labels:
            f_names.append("output_files/" + str(year) +"/" + prefix + str(year % 2000) + "_" + label + ending + ".root") 

        rm_name = "output_files/" + str(year)+"/" + prefix +str(year %2000) + "_fakes_contam" + "*.root"
        os.system("rm %s" % rm_name)
        fakes_name = "output_files/" + str(year)+"/" + prefix +str(year %2000) + "_fakes_contam" + ending + ".root"
        cmd = "hadd " + fakes_name
        for name in f_names:
            cmd += " " + name
        print(cmd)
        os.system(cmd)
        os.system("./scripts/prune_root_file.sh %s" % fakes_name)

