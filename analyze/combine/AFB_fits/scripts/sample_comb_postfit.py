import ROOT
from utils import *
from ROOT import *
from utils import *
from parse import parse



def generateTripletNames(f):
    # Based on naming scheme, group histogram years together in tuples of length 3
    base_keys = [k.GetName() for k in f.GetListOfKeys() if '16' in k]
    for k in base_keys:
        return [k, k.replace('16','17'), k.replace('16','18')]

def getHistsFromTriplet(f,triplet):
    out = []
    for d in triplet:
        f.cd()
        hname = d + "TotalProcs"
        out.append(f.Get(hname))
    return out
    #return [f.Get(d + "TotalProcs") for d in triplet]

def run2FromTriplet(f,triplet):
    hists = getHistsFromTriplet(f,triplet)
    out = hists[0].Clone(hists[0].GetName()+'_Run2')
    out.Reset()
    for h in hists:
        out.Add(h)
    return out

def sample_comb_postfit(workspace, fitname, output_name, nSample = 100):
    for i in range(nSample):
        print_and_do('PostFitShapesFromWorkspace -w %s -f %s -o sample_%s.root --seed %i --skip-prefit --postfit --samples 1 >& temp_log' % (workspace, fitname, i, i))
    os.system("rm temp_log")


    dirs = [("Y16_mumu16_postfit/", "Y17_mumu17_postfit/", "Y18_mumu18_postfit/"), 
            ("Y16_ee16_postfit/", "Y17_ee17_postfit/", "Y18_ee18_postfit/")]


    out = TFile.Open(output_name,'UPDATE')

    for triplet in dirs:
        print(triplet)
        dummy_file = TFile.Open('sample_0.root','OPEN')
        dummy = run2FromTriplet(dummy_file,triplet)
        dummy.SetDirectory(0)
        maximum_content = int(1.5*dummy.GetMaximum())
        sample_storage_hist = TH2F('result_storage','',
                                   dummy.GetNbinsX(),
                                   dummy.GetXaxis().GetXmin(),
                                   dummy.GetXaxis().GetXmax(),
                                   1000*maximum_content,0,maximum_content) # resolution will be 1/100ths
        sample_storage_hist.SetDirectory(0)
        dummy_file.Close()

        for i in range(nSample):
            print(i)
            f = TFile.Open('sample_%s.root'%i,'OPEN')
            run2 = run2FromTriplet(f, triplet)
            run2.SetDirectory(0)
            run2.Print("range")

            for ibin in range(1,run2.GetNbinsX()+1):
                sample_storage_hist.Fill(run2.GetXaxis().GetBinCenter(ibin),run2.GetBinContent(ibin))
            f.Close()

        total_hist = TH1F(triplet[0].replace('Y16','Run2').replace('16', '').replace("/", ""), '',
                          dummy.GetNbinsX(),
                          dummy.GetXaxis().GetXmin(),
                          dummy.GetXaxis().GetXmax())

        total_hist.SetDirectory(0)

        for ibin in range(1,total_hist.GetNbinsX()+1):
            temp_projy = sample_storage_hist.ProjectionY('temp_projy',ibin,ibin)
            total_hist.SetBinContent(ibin, temp_projy.GetMean())
            total_hist.SetBinError(ibin, temp_projy.GetRMS())

        total_hist.Print("range")

        out.cd()
        total_hist.Write()
    print_and_do("rm sample_*.root")


if(__name__ == "__main__"):
    parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
    parser.add_option("-n", "--nSample",  default=100, type=int, help="How many samples")
    parser.add_option("-w",  "--workspace", default= "", help="Workspace")
    parser.add_option("-f",  "--fitname", default= "", help="fitname")
    parser.add_option("-o",  "--output", default= "", help="output")
    (options, args) = parser.parse_args()

    sample_comb_postfit(options.workspace, options.fitname, options.output, options.nSample)

