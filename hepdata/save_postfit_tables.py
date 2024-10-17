import ROOT
from ROOT import *
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

import pandas as pd

m_str, y_str, c_str = [],[],[]
m = [500,700,1000,14000]
y = [0.,0.6,1.,2.4]
c = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]

for i in range(len(m)-1):
    for j in range(len(y)-1):
        for k in range(len(c)-1):
	    if y[j] >= 0.6 and (c[k] == -0.75 or c[k] == 0.75): continue
            m_str.append("%i-%i" % (m[i], m[i+1]))
            y_str.append("%.1f-%0.1f" % (y[j], y[j+1]))
            if j>=1 and (k==0 or k==6): c_str.append("%.2f-%.2f" % (c[k],c[k+2]))
	    else: c_str.append("%.2f-%.2f" % (c[k],c[k+1]))

vec = True
for chan in ["ee","mumu"]:
    for q in ["u","d"]:
        f_in = TFile.Open("analyze/combine/AFB_fits/postfit_plots/%s_%s%s_unblinded_LQ_m2500/%s_%s%s_unblinded_fit_shapes_LQ.root"%(chan,q,"_vec" if vec else "",chan,q,"_vec" if vec else ""))
        for year in [16,17,18]:
            
            dir_ = "Y%i_postfit/"%year
            print ("\n dir_ = ", dir_)
            
            h_totx = f_in.Get(dir_ + "TotalBkg")
            h_datax = f_in.Get(dir_ + "data_obs")
            h_sigx = f_in.Get(dir_ + "TotalSig")

            if year == 16:
                h_tot = h_totx.Clone()
                h_data = h_datax.Clone()
                h_sig = h_sigx.Clone()
            else:
                h_tot.Add(h_totx)
                h_data.Add(h_datax)
                h_sig.Add(h_sigx)
            
        total_proc, total_proc_err = [],[]
        data, data_err = [],[]
        sig, sig_err = [],[]

        for i in range(h_tot.GetNbinsX()):

            total_proc.append(h_tot.GetBinContent(i))
            total_proc_err.append(h_tot.GetBinError(i))

            data.append(h_data.GetBinContent(i))
            data_err.append(h_data.GetBinError(i))

            sig.append(h_sig.GetBinContent(i))
            sig_err.append(h_sig.GetBinError(i))

        f_in.Close()

        df = pd.DataFrame({
            'Bin index': range(1,61),
            'Background': total_proc,
            'Background err': total_proc_err,
            'LQ template yield': sig,
            'LQ template yield err': sig_err,
            'Observed data': data,
            'Observed data err': data_err
        })
	df = df[['Bin index','Background','Background err','LQ template yield','LQ template yield err','Observed data','Observed data err']]
        print(df)
        df.to_csv("hepdata/%s_%s%s_postfit_table.txt" % (chan, q, "_vec" if vec else ""),sep=' ',index=False, header=False)
                


            

        
