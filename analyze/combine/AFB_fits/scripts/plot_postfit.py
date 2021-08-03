import ROOT
from ROOT import *
from contextlib import contextmanager
import os, pickle, subprocess, time,random
import math
from math import sqrt
from array import array
from optparse import OptionParser
import CMS_lumi, tdrstyle


gStyle.SetOptStat(0)
gROOT.SetBatch(1)

def makeCan(name, tag, histlist, bkglist=[],signals=[],totlist = [], colors=[],titles=[],dataName='data',bkgNames=[],signalNames=[],
        logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False,datastyle='pe',year=1):  
    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]

    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1600
        height = 700
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        width = 1600
        height = 1000
        padx = 3
        pady = 2
    else:
        print 'histlist of size ' + str(len(histlist)) + ' not currently supported'
        print histlist
        return 0

    tdrstyle.setTDRStyle()

    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [kRed,kMagenta,kGreen,kCyan,kBlue]
    if len(colors) == 0:   
        colors = default_colors
    stacks = []
    legends = []
    legends_list = []
    mains = []
    subs = []
    pulls = []
    logString = ''
    leg_align_right = True
    CMS_align_right = False

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        # if len(histlist) > 1:
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd()        

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            if logy == True:
                gPad.SetLogy()
            gPad.SetLeftMargin(0.2)
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            hist.GetXaxis().SetTitleOffset(1.5)
            hist.GetYaxis().SetTitleOffset(1.0)
            hist.GetZaxis().SetTitleOffset(1.8)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            hist.Draw('lego')
            if len(bkglist) > 0:
                print 'ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.'
        
        # Otherwise it's a TH1 hopefully
        else:
            titleSize = 0.09
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(kBlack,alpha)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(kBlack,alpha)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)

            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            
            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.SetMaximum(1.13*hist.GetMaximum())
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                    hist.SetTitleOffset(1.1)
                hist.Draw(datastyle)
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.3, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.3))

                else:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                leg_align_right = True
                CMS_align_right = False
                x_max = totlist[hist_index].GetMaximumBin()
                nbins = totlist[hist_index].GetXaxis().GetNbins()
                #if(2 *x_max > nbins):
                #    print("Found max val in bin %i, aligning legend on the left" % x_max)
                #    leg_align_right = False
                #    CMS_align_right = True
                if not logy: 
                    y_end = 0.88
                    y_size = 0.2 + 0.02*(len(bkglist[0])+len(signals))
                    x_size = 0.35
                    if(leg_align_right):
                        x_start = 0.55
                    else:
                        x_start = 0.2

                    legends.append(TLegend(x_start,y_end - y_size,x_start + x_size,y_end))
                else: 
                    legends.append(TLegend(0.2,0.11,0.45,0.2+0.02*(len(bkglist[0])+len(signals))))

                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                legends_list.append([])


                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.17)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.08)

                subs[hist_index].SetLeftMargin(0.17)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0)
                subs[hist_index].SetBottomMargin(0.35)
                mains[hist_index].Draw()
                subs[hist_index].Draw()

                # Build the stack
                for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                    # bkg.Sumw2()
                    bkg.SetLineColor(kBlack)
                    if logy:
                        bkg.SetMinimum(1e-3)

                    if colors[bkg_index] != None:
                        bkg.SetFillColor(colors[bkg_index])
                    else:
                        bkg.SetFillColor(default_colors[bkg_index])

                    stacks[hist_index].Add(bkg)
                    if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                    elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                    else: this_bkg_name = bkgNames[hist_index][bkg_index]
                    legends_list[hist_index].append((bkg,this_bkg_name,'f'))
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()


                # Set y max of all hists to be the same to accomodate the tallest
                max_scaling = 2.0
                histList = [stacks[hist_index],totlist[hist_index],hist]

                yMax = histList[0].GetMaximum()
                maxHist = histList[0]
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                        maxHist = histList[h]
                for h in histList:
                    h.SetMaximum(yMax*max_scaling)
                    if logy == True:
                        h.SetMaximum(yMax*10)
                    else:
                        h.SetMinimum(0.)

                
                mLS = 0.08
                mTS = 0.06
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.GetYaxis().SetTitleOffset(1.4)
                hist.GetXaxis().SetTitleOffset(1.2)
                hist.GetYaxis().SetTitle('Events / bin')
                hist.GetYaxis().SetLabelSize(mLS)
                hist.GetYaxis().SetTitleSize(mTS)
                hist.GetXaxis().SetLabelOffset(999)
                if logy == True:
                    hist.SetMinimum(1e-3)

                hist.SetBinErrorOption(ROOT.TH1.kPoisson)
                hist.Draw(datastyle)
                #print("Drawing %s %s \n" hist.GetName(), datastyle)

                stacks[hist_index].Draw('same hist')
                #print("Drawing %s same hist \n" stacks[hist_index].GetName())

                # Do the signals
                if len(signals) > 0: 
                    signals[hist_index].SetLineColor(kBlue)
                    signals[hist_index].SetLineWidth(2)
                    if logy == True:
                        signals[hist_index].SetMinimum(1e-3)
                    if signalNames == []: this_sig_name = signals[hist_index].GetName().split('_')[0]
                    legends_list[hist_index].append((signals[hist_index],this_sig_name,'L'))
                    signals[hist_index].Draw('hist same')

                totlist[hist_index].SetFillColor(kBlack)
                totlist[hist_index].SetFillStyle(3354)
                totlist[hist_index].SetMarkerStyle(20)
                totlist[hist_index].SetMarkerSize(0.01)

                totlist[hist_index].Draw('e2 same')

                if not dataOff:
                    legends_list[hist_index].append((hist,dataName,datastyle))
                    hist.Draw(datastyle+' same')


            
                legends[hist_index].SetHeader(title, "c")
                legends[hist_index].SetNColumns(2)
                
                for entry in legends_list[hist_index][::-1]:
                    legends[hist_index].AddEntry(entry[0], entry[1], entry[2])



                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()
                gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
		ratio, ratio_sys_unc = makeRatio(hist,totlist[hist_index])
                ratio.Print()
                hist.Print("range")
                chi2 = 0.
                #for i in range(1, pull.GetNbinsX()+1):
                    #chi2 += pull.GetBinContent(i)**2;
                #print("Chi2/nbin for chan %s is %.1f/%i" % (titles[hist_index], chi2, pull.GetNbinsX()))

                LS = .13
                #title size given as fraction of pad width, scale up to have same size as main pad
                TS =  0.06 * 0.7/0.3



                ratio.GetYaxis().SetRangeUser(0.0, 2.0)
                ratio.GetYaxis().SetTitleOffset(0.3)
                ratio.GetXaxis().SetTitleOffset(1.2)
                             
                ratio.GetYaxis().SetLabelSize(LS)
                ratio.GetYaxis().SetTitleSize(TS)
                ratio.GetYaxis().SetNdivisions(306)
                ratio.GetYaxis().SetTitle("Data/Fit")

                ratio.GetXaxis().SetRangeUser(0., hist.GetNbinsX()-.08)
                ratio.GetXaxis().SetLabelOffset(0.05)
                ratio.GetXaxis().SetLabelSize(LS)
                ratio.GetXaxis().SetTitleSize(TS)
                ratio.GetXaxis().SetTitle(xtitle)


                ratio_sys_unc.SetFillColor(ROOT.kBlack)
                totlist[hist_index].SetFillStyle(3354)


                ratio.Draw('AP')
                line = TLine(0, 1.0, hist.GetNbinsX() - 0.08, 1.0)
                line.Draw()
                ratio_sys_unc.Draw("3 same")

                if logy == True:
                    mains[hist_index].SetLogy()

                if(CMS_align_right): CMS_loc = 33
                else: CMS_loc = 11
                CMS_lumi.CMS_lumi(thisPad, year, CMS_loc)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'.pdf','pdf')
        myCan.Print(tag+'/'+name+'.png','png')


def reducedCorrMatrixHist(fit_result,varsOfInterest=[]):
    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if 'Fail_' not in finalPars.at(cm_index).GetName():
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = open('correlation_matrix.txt','w')

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]
        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt.write('%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y]))
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out

def FindCommonString(string_list):
    to_match = ''   # initialize the string we're looking for/building
    for s in string_list[0]:    # for each character in the first string
        passed = True
        for istring in range(1,len(string_list)):   # compare to_match+s against strings in string_list
            string = string_list[istring]
            if to_match not in string:                  # if in the string, add more
                passed = False
            
        if passed == True:
            to_match+=s

    if to_match[-2] == '_':
        return to_match[:-2] 
    else:
        return to_match[:-1]                # if not, return to_match minus final character

    return to_match[:-2]
        
def makeRatio( DATA,BKG):

    nbins = DATA.GetNbinsX()
    x = array('d')
    ratio = array('d')
    y_err_up = array('d')
    y_err_down = array('d')
    sys_err_up = array('d')
    sys_err_down = array('d')
    x_err_low = array('d')
    x_err_high = array('d')
    x_err_zero = array('d', [0.]* nbins)
    y_val1 = array('d', [1.]* nbins)
    
    for ibin in range(1,DATA.GetNbinsX()+1):
        DATAcont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        BKG_errup = BKG.GetBinErrorUp(ibin)
        BKG_errdown = BKG.GetBinErrorLow(ibin)
        DATA_errup = DATA.GetBinErrorUp(ibin)
        DATA_errdown = DATA.GetBinErrorLow(ibin)


        x.append(DATA.GetXaxis().GetBinCenter(ibin))
        x_err_low.append(DATA.GetXaxis().GetBinCenter(ibin) - DATA.GetXaxis().GetBinLowEdge(ibin))
        x_err_high.append(-DATA.GetXaxis().GetBinCenter(ibin) + DATA.GetXaxis().GetBinUpEdge(ibin))

    

        ratio.append(DATAcont/BKGcont)
        y_err_up.append(DATA_errup / BKGcont)
        y_err_down.append(DATA_errdown / BKGcont)
        sys_err_up.append(BKG_errup/BKGcont)
        sys_err_down.append(BKG_errdown/BKGcont)
        
    pull = ROOT.TGraphAsymmErrors(nbins,x,ratio, x_err_zero, x_err_zero, y_err_down, y_err_up)
    sys_unc = ROOT.TGraphAsymmErrors(nbins, x, y_val1, x_err_low, x_err_high, sys_err_down,  sys_err_up)
    return pull, sys_unc

def Make_up_down(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down

if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("--input", "-i", default = "", help="Input file")
    parser.add_option("--output", "-o", default = "", help="Input directory")
    parser.add_option("--mbin", "-m", type = 'int', default = 0, help="Mass bin (for plot label)")
    parser.add_option("--year", "-y", type = 'int', default = -1, help="Year (-1 for all) ")
    parser.add_option("--ss",   default = False, action='store_true',  help="Fit was done with ee_ss region too")
    (options, args) = parser.parse_args()


#fin_ = "combined_fit_shapes_mbin1.root"
#odir = "postfit_plots/combined_fit_mbin1"
#mbin = 1
    if(options.year < 0):
        years = [2016, 2017, 2018]
    else:
        years = [options.year]
    h_names = ["gam", "db", "qcd", "top",  "tautau", "dy"]
    h_ss_names = ["bk", "dy", "qcd"]


    m_bins = [150, 170, 200,  250, 320, 510, 700, 1000, 14000]


    label_color_map = dict()
    label_color_map['dy'] = ("DY Signal", kRed + 1)
    label_color_map['top'] = ("t#bar{t} + tW ", kBlue)
    label_color_map['db'] = ("WW + WZ + ZZ",  kGreen +3)
    label_color_map['tautau'] = ("DY #tau#tau Bkg.", kMagenta + 4)
    label_color_map['gam'] = ("\\gamma\\gamma \\rightarrow \\ell\\ell ", kOrange)
    label_color_map['qcd'] = ("WJets + QCD", kRed - 7)

    datastyle = "pe0x0"

    fracs = dict()
    for name in h_names:
        fracs[name] = 0.

    dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/"]
    if(options.ss): dirs = dirs = ["Y%i_mumu%i_postfit/", "Y%i_ee%i_postfit/", "Y%i_ee%i_ss_postfit/"]
    f_in = TFile.Open(options.input)
    for year in years:
        for idx, dir_name in enumerate(dirs):
            dir_ = dir_name % (year % 2000, year % 2000)
            h_tot = f_in.Get(dir_ + "TotalProcs")
            h_tot = h_tot.Clone("h_tot_c%i_y%i" %(idx, year))
            h_data_orig = f_in.Get(dir_ + "data_obs")
            #copy to new hist so poisson error bars work
            h_data = h_data_orig.Clone("h_data_c%i_y%i" %(idx, year))
            h_data.Reset()
            for b in range(h_data_orig.GetXaxis().GetNbins() + 1):
                h_data.SetBinContent(b, h_data_orig.GetBinContent(b))

            h_tot_sig = f_in.Get(dir_ + "TotalSig")
            h_tot_sig = h_tot_sig.Clone("h_tot_sig_c%i_y%i" %(idx, year))

            mbin_low = m_bins[options.mbin]
            mbin_high = m_bins[options.mbin+1]

            if(idx == 0): title = "Muons %i-%i GeV" % (mbin_low, mbin_high)
            if(idx == 1): title = "Electrons %i-%i GeV" % (mbin_low, mbin_high)
            if(idx == 2): title = "Electrons Samesign %i-%i GeV" % (mbin_low, mbin_high)
            
            if(idx == 2): name_list = h_ss_names
            else: name_list = h_names
            hist_list = []
            color_list = []
            label_list = []

            for name in name_list:
                if(name == "dy"):
                    h = h_tot_sig.Clone("h_%s_c%i_y%i" %(name, idx, year))
                else:
                    h = f_in.Get(dir_ + name)
                    if(h != None):
                        h = h.Clone("h_%s_c%i_y%i" %(name, idx, year))
                if(h != None):
                    #h.Print()
                    hist_list.append(h)
                    label_list.append(label_color_map[name][0])
                    color_list.append(label_color_map[name][1])

                    if("gam" in name):
                        this_frac = h.Integral()/(h_tot_sig.Integral() + h.Integral())
                        print("Chan %i Year %i Name %s frac %.3f \n" % (idx, year, name, this_frac))
                        fracs[name] += this_frac

            makeCan(dir_[:-1], options.output, [h_data], bkglist=[hist_list], totlist=[h_tot], colors = color_list, bkgNames = label_list, 
                    titles = [title], xtitle = "Template Bin", year = year, datastyle=datastyle ) 

    for key in fracs.keys():
        fracs[key] /= (len(years)*len(dirs))

    key_ = "gam"
    print("Average fraction for %s  is  %.3f \n" % (key_, fracs[key_]))


