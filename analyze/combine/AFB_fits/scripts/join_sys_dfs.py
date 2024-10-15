import pandas as pd
import math
import numpy as np

mLQ = 2500
ending = "091224"


for is_vec in [False, True]:
    if is_vec: ending+="_vec" 
    df_list=[]
    flag = 0
    for chan in ["mumu","ee"]:
        for q in ["u","d"]:
            if not flag:
                df_list = pd.read_csv("%s_%s_m%s_sys_uncs_%s_statuncs.txt"%(chan, q, mLQ, ending),delimiter=" ", header=0,names=["Sys name","Contri_mumu_u","Percent_mumu_u"],dtype={"Sys name":str,"Contri_mumu_u":np.float64,"Percent_mumu_u":np.float64})
                df_list.drop(["Contri_mumu_u"], axis=1, inplace=True)
                flag = 1
            else:
                df0 = pd.read_csv("%s_%s_m%s_sys_uncs_%s_statuncs.txt"%(chan, q, mLQ, ending),delimiter=" ", header=0,names=["Sys name","Contri","Percent"],dtype={"Sys name":str,"Contri_mumu_u":np.float64,"Percent_mumu_u":np.float64})
                df0.drop(["Contri"], axis=1, inplace=True)
                df_list = df_list.merge(df0, on='Sys name', how="right", suffixes=(None,"_%s_%s%s"%(chan, q, ("_vec" if is_vec else ""))))
    df_list = df_list.round(2)
    print(df_list)
    
    df_list.to_csv("m%s_sysuncs_statuncs_%s.txt"%( mLQ, ending),sep=' ',index=False)




