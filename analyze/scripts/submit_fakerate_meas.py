
import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

def print_and_do(s):
    print(s)
    return os.system(s)

years = [2016, 2017,2018]
prefixes = ["FakeRate"]
scripts = ["SingleElectron_data_measure_fake_rate.C", "SingleElectron_mc_contam_fake_rate.C",
          "SingleMuon_data_measure_fake_rate.C", "SingleMuon_mc_contam_fake_rate.C"]
njobs = 15

eos_files = ["SingleElectron_files.txt", "non_QCD_files.txt", "SingleMuon_files.txt", "non_QCD_files.txt" ]
labels = ["el_data", "el_contam", "mu_data", "mu_contam"]
ending = "oct26"

for i,script in enumerate(scripts):
    prefix = prefixes[0]
    for year in years:
            script_name = "scripts/script1.sh"
            print_and_do("cp scripts/script_template.sh %s" % script_name)
            print_and_do("sed -i  s/PREFIX/%s/g %s" % (prefix, script_name))
            print_and_do("sed -i  s/EXEC/%s/g %s" % (script, script_name))
            print_and_do("sed -i  s/YEAR/%i/g %s" % (year, script_name))
            print_and_do("sed -i  s/EOSINPUT/%s/g %s" % (eos_files[i], script_name))
            print_and_do("chmod +x %s" % script_name)
            print_and_do("python doCondor.py --njobs %i --sub  -s %s -n %s%i_%s_%s"  % (njobs, script_name, prefix, year % 2000, labels[i], ending))
            print_and_do("rm %s" % script_name)


