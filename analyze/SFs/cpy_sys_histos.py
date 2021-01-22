from ROOT import *

f_sys_name = "2016/tmp/RunGH_SYS_ID.root"
f_out_name = "2016/tmp/Mu_GH_ID.root"

f_sys = TFile.Open(f_sys_name, "READ")
f_out = TFile.Open(f_out_name, "UPDATE")

f_sys.cd()
keys = gDirectory.GetListOfKeys()
for key in keys:
    key_name = key.GetName()
    if('_syst' not in key_name):
        new_name = key_name + "_syst"
    else:
        new_name = key_name

    f_sys.cd()
    h = gDirectory.Get(key_name)
    h.SetName(new_name)
    f_out.cd()
    h.Write()

f_out.cd()
gDirectory.ls()

