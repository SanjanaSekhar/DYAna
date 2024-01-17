import numpy as np
import json
from math import sqrt

d1 = "091823"
d2 = "011524"


for is_vec in [True, False]:
    for channel in ['ue','de','um','dm']:

        lim_diff = {}
        with open("LQ_cards/%s/limit_json/limits_%s%s_%s.json"%(channel,channel, ("_vec" if is_vec else ""), d1), 'r+') as f:
            data1 = json.load(f)

        with open("LQ_cards/%s/limit_json/limits_%s%s_%s.json"%(channel,channel, ("_vec" if is_vec else ""), d2), 'r+') as f:
            data2 = json.load(f)

        for m in range(1000,5500,500):
            lim_diff = data1
            avg = 0
            for lim in data1[str(m)+'.0']:
                
                if lim=='obs': 
                    lim_diff[str(m)+'.0']['avg'] = lim_diff[str(m)+'.0']['obs']
                    del lim_diff[str(m)+'.0']['obs']
                    lim_diff[str(m)+'.0']['avg'] = avg

                elif lim!='obs' and lim!='avg':
                    lim_diff[str(m)+'.0'][lim] = (data1[str(m)+'.0'][lim] - data2[str(m)+'.0'][lim])/data1[str(m)+'.0'][lim]
                    avg += lim_diff[str(m)+'.0'][lim]/5.

        with open("LQ_cards/%s/limit_json/compare_limits_%s%s.json"%(channel,channel, ("_vec" if is_vec else "")), 'w') as f:
            f.seek(0)
            json.dump(lim_diff, f, indent=4)
