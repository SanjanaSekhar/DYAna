
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

years = [2016, 2017,2018]
prefixes = ["MuMu", "ElEl"]
#scripts = ["MuMu_reco_background.C", "ElEl_reco_background.C"]
scripts = ["MuMu_reco_mc.C", "ElEl_reco_mc.C"]
njobs = 15

#labels = ["wt", "ttbar", "diboson", "phot_ind"]
#eos_files = ["WT_files.txt", "TTbar_files.txt", "diboson_files.txt", "PhotInd_files.txt"]
labels = ["dy"]
eos_files = ["DY_files.txt"]
ending = "april17"

for i,script in enumerate(scripts):
    prefix = prefixes[i]
    for year in years:
            f_names = []
            for j,label in enumerate(labels):
                script_name = "scripts/script1.sh"
                print_and_do("cp scripts/script_template.sh %s" % script_name)
                print_and_do("sed -i  s/PREFIX/%s/g %s" % (prefix, script_name))
                print_and_do("sed -i  s/EXEC/%s/g %s" % (script, script_name))
                print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
                print_and_do("sed -i  s/EOSINPUT/%s/g %s" % (eos_files[j], script_name))
                print_and_do("chmod +x %s" % script_name)
                print_and_do("python doCondor.py --njobs %i --sub  -s %s -n %s%i_%s_%s"  % (njobs, script_name, prefix, year % 2000, label, ending))
                print_and_do("rm scripts/script1.sh")


