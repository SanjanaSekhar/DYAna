


import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

years = [2016,2017,2018]
prefixes = ["MuMu", "ElEl"]
labels = ["dy", "wt", "ttbar", "diboson"]
ending = "_mar30"

for year in years:
    for prefix in prefixes:
        f_names = []
        for label in labels:
            os.system("python doCondor.py -e -o output_files/" + str(year) +"/ -n " + prefix + str(year % 2000) + "_" + label + ending) 
            f_names.append("output_files/" + str(year) +"/" + prefix + str(year % 2000) + "_" + label + ending + ".root") 

        fakes_name = "output_files/" + str(year)+"/" + prefix +str(year %2000) + "_fakes_contam" + ending + ".root"
        cmd = "hadd " + fakes_name
        for name in f_names:
            cmd += " " + name
        print(cmd)
        os.system(cmd)
        os.system("./scripts/prune_root_file.sh %s" % fakes_name)

