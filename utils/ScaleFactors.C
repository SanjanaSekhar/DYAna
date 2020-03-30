#ifndef DYAFB_SFS
#define DYAFB_SFS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "HistUtils.C"





typedef struct {
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
    TH2D *ISO_SF;
    TH2D *ID_SF;
} mu_SFs;


typedef struct {
    TH1D *data_pileup;
    TH1D *pileup_ratio;
} pileup_SFs;

typedef struct {
    TH1D *data_pileup_up;
    TH1D *data_pileup_down;
    TH1D *data_pileup_nom;

    TH1D *ratio_pileup_up;
    TH1D *ratio_pileup_down;
    TH1D *ratio_pileup_nom;
} pileup_systematics;

typedef struct{
    TH2D *ID_SF;
    TH2D *RECO_SF;
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
} el_SFs;

typedef struct{
    TH2F *jet_rate;
    TH2F *el_rate;
} prefire_SFs;


double get_var(Double_t vals[100]){
    float mean(0.), var(0.);
    int n_vars = 100;
    int n_entries = n_vars;
    for(int i =0; i< n_vars; i++){
        //printf("%.2f \n", vals[i]);
        if(std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    //printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vars; i++){
        if(std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    //printf("std %.3f \n\n\n", sqrt(var));
    return var;
}

Double_t get_pileup_SF(Int_t n_int, TH1D *h){
    if(n_int > 99) n_int = 99;

    TAxis* x_ax =  h->GetXaxis();
    int xbin = x_ax->FindBin(n_int);

    Double_t result = h->GetBinContent(xbin);
    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}

Double_t get_Mu_trk_SF(Double_t eta, TGraphAsymmErrors *h, int systematic = 0){
    eta = abs(eta);
    Double_t result = h->Eval(eta);
    
    if(systematic !=0){
        TAxis* x_ax =  h->GetXaxis();
        int xbin = x_ax->FindBin(eta);
        Double_t err;
        if(systematic >0) err = h->GetErrorYhigh(xbin);
        else if(systematic <0) err = h->GetErrorYlow(xbin);
        result += err*systematic;
    }

    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}

void get_pdf_avg_std_dev(Float_t pdf_Weights[100], Float_t *pdf_avg, Float_t *pdf_std_dev){
    Float_t sum = 0;
    for (int i=0; i<100; i++){
        sum+= pdf_Weights[i];
    }
    *pdf_avg = sum/100.;
    Float_t var = 0;
    for (int i=0; i<100; i++){
        var += pow(*pdf_avg - pdf_Weights[i], 2);
    }
    *pdf_std_dev = sqrt(var/99.);
    return;
}

Double_t get_prefire_rate(float pt, float eta, TH2F *map, int systematic = 0){
    //based on
    //https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/PrefireCorr.py
    float min_pt = 20.;
    float max_pt = 499.;
    float min_eta = 2.0;
    float max_eta = 3.0;
    if((pt < min_pt) || (std::abs(eta) < min_eta) || (std::abs(eta) > max_eta)) return 0.;

    int bin = map->FindBin(eta, std::min(pt, max_pt));
    float prefire_prob = map->GetBinContent(bin);

    float stat_err = map->GetBinError(bin);
    float syst_err = 0.2 * prefire_prob;
    if(systematic == -1){
        prefire_prob = std::max(0., prefire_prob - std::pow(stat_err*stat_err + syst_err*syst_err, 0.5));
    }
    else if(systematic == 1){
        prefire_prob = std::min(1., prefire_prob + std::pow(stat_err*stat_err + syst_err*syst_err, 0.5) );
    }
    //printf("pt %.0f eta %.2f, prob %.3f \n", pt, eta, prefire_prob);
    return prefire_prob;
}





Double_t get_mu_SF(Double_t pt, Double_t eta, int year, TH2D *h, int systematic = 0){
    if(year <2016 || year > 2018) printf("Year is not from 2016-2018. This is bad!! \n");
    //stay in range of histogram
    float sys_unc_mult = 1.0;
    if (pt >= 115.){
        sys_unc_mult = 1.5;
        pt = 90.;
    }
    if (pt <= 22.5) pt = 22.5;
    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin(0),ybin(0);
    if(year == 2016){
        //2016 hists are eta and pt
        xbin = x_ax->FindBin(eta);
        ybin = y_ax->FindBin(pt);
    }
    else{
        //2017 & 18 hists are pt and abs(eta)
        xbin = x_ax->FindBin(pt);
        ybin = y_ax->FindBin(abs(eta));
    }


    Double_t result = h->GetBinContent(xbin, ybin);
    if(systematic != 0){
        Double_t err = sys_unc_mult * h->GetBinError(xbin, ybin);
        //printf("SF is %.3f +/- %.3f \n", result, err);
        result += (systematic * err );
    }
    if(result < 0.001){ 
        printf("0 muon SF for Pt %.1f, Eta %1.2f \n", pt, eta);
        result = 1;
    }
    return result;
}

Double_t get_el_SF(Double_t pt, Double_t eta, TH2D *h, int systematic_barrel = 0, int systematic_endcap = 0){
    float sys_unc_mult = 1.0;
    //get SF for eta's in overlap region (already vetoed superclusters in the
    //region)
    if(!goodElEta(eta)) eta = eta * 1.2;

    if( pt <= 25.) pt = 25;
    if (pt >= 450.){
        pt = 350.;
        sys_unc_mult = 1.5;
    }
    TAxis* x_ax =  h->GetXaxis();

    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);
    Double_t result = h->GetBinContent(xbin, ybin);
    int systematic = systematic_barrel;
    if(abs(eta) > 1.4) systematic = systematic_endcap;
    if(systematic != 0){
        Double_t err = sys_unc_mult * h->GetBinError(xbin, ybin);
        err =abs(err);
        result += (systematic * err);
    }
    //printf("eta %.2f sys_b %i sys_e %i, sys %i \n", eta, systematic_barrel, systematic_endcap, systematic);

    if(result < 0.01 || result > 10.){
        printf("%.2f el SF for Pt %.1f, Eta %1.2f \n", result, pt, eta);
        result = 1;
    }
    return result;
}
Double_t get_HLT_SF_1mu(Double_t mu1_pt, Double_t mu1_eta, TH2D *h_SF){
    //get HLT SF for event with just 1 muon
    //stay in range of histogram
    if (mu1_pt >= 350.) mu1_pt = 350.;
    mu1_eta = abs(mu1_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(std::abs(mu1_eta));
    int ybin1_SF = y_ax_SF->FindBin(mu1_pt);


    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);

    Double_t result = SF1;
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan HLT SF 1 mu for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Double_t get_HLT_SF_1el(Double_t el1_pt, Double_t el1_eta, TH2D *h_SF){
    //get HLT SF for event with just 1 elon
    //stay in range of histogram
    if (el1_pt >= 350.) el1_pt = 350.;
    el1_eta = abs(el1_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(el1_eta);
    int ybin1_SF = y_ax_SF->FindBin(el1_pt);


    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);

    Double_t result = SF1;
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", el1_pt, el1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan HLT SF for Pt %.1f, Eta %1.2f \n", el1_pt, el1_eta);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Double_t get_HLT_SF(Double_t lep1_pt, Double_t lep1_eta, Double_t lep2_pt, Double_t lep2_eta, TH2D *h_SF, TH2D *h_MC_EFF, 
       int systematic_barrel = 0, int systematic_endcap = 0){
    float sys1_unc_mult = 1.0;
    float sys2_unc_mult = 1.0;
    //printf("Getting HLT for %.2f %.2f %.2f %.2f \n", lep1_pt, lep1_eta, lep2_pt, lep2_eta);
    //Get HLT SF for event with 2 leptons
    //stay in range of histogram
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int pt_max_bin = y_ax_SF->GetLast();
    double pt_max = y_ax_SF->GetBinUpEdge(pt_max_bin);
    double pt_overflow = y_ax_SF->GetBinCenter(pt_max_bin);
    if (lep1_pt >= pt_max){
        sys1_unc_mult = 1.5;
        lep1_pt = pt_overflow;
    }
    if (lep2_pt >= pt_max){
        sys2_unc_mult = 1.5;
        lep2_pt = pt_overflow;
    }
    lep1_eta = abs(lep1_eta);
    lep2_eta = abs(lep2_eta);
    int xbin1_SF = x_ax_SF->FindBin(std::fabs(lep1_eta));
    int ybin1_SF = y_ax_SF->FindBin(lep1_pt);

    int xbin2_SF = x_ax_SF->FindBin(std::fabs(lep2_eta));
    int ybin2_SF = y_ax_SF->FindBin(lep2_pt);

    Double_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);
    Double_t SF2 = h_SF->GetBinContent(xbin2_SF, ybin2_SF);

    int systematic1 = systematic_barrel;
    int systematic2 = systematic_barrel;
    if(abs(lep1_eta) > 1.4) systematic1 = systematic_endcap;
    if(abs(lep2_eta) > 1.4) systematic2 = systematic_endcap;

    if(systematic1 != 0){
        Double_t SF1_err = h_SF->GetBinError(xbin1_SF, ybin1_SF);
        SF1 += sys1_unc_mult * SF1_err * systematic1;
    }
    if(systematic2 != 0){
        Double_t SF2_err = h_SF->GetBinError(xbin2_SF, ybin2_SF);
        SF2 += sys2_unc_mult * SF2_err * systematic2;
    }
    //printf("eta1 %.2f eta2 %.2f sys_b %i sys_e %i sys1 %i sys2 %i \n", lep1_eta, lep2_eta, systematic_barrel, systematic_endcap, systematic1, systematic2);



    TAxis *x_ax_MC_EFF =  h_MC_EFF->GetXaxis();
    TAxis *y_ax_MC_EFF =  h_MC_EFF->GetYaxis();
    int xbin1_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(lep1_eta));
    int ybin1_MC_EFF = y_ax_MC_EFF->FindBin(lep1_pt);

    int xbin2_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(lep2_eta));
    int ybin2_MC_EFF = y_ax_MC_EFF->FindBin(lep2_pt);

    Double_t MC_EFF1 = h_MC_EFF->GetBinContent(xbin1_MC_EFF, ybin1_MC_EFF);
    Double_t MC_EFF2 = h_MC_EFF->GetBinContent(xbin2_MC_EFF, ybin2_MC_EFF);
    Double_t result = (1 - (1-MC_EFF1*SF1)*(1-MC_EFF2*SF2))/
                      (1 - (1-MC_EFF1)*(1-MC_EFF2));
    //printf("%.1f %.1f SF: %.3f %.3f, MC EFF: %.3f %.3f  comb %.3f \n", lep1_pt, lep2_pt, SF1, SF2, MC_EFF1, MC_EFF2, result);
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", lep1_pt, lep1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan lep HLT SF for Pt1 %.2f, Eta1 %.2f Pt2 %.2f Eta2 %.2f  %.2f %.2f \n", lep1_pt, lep1_eta, lep2_pt, lep2_eta, SF1, SF2);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}





void setup_pu_SFs(pileup_SFs *pu_SF, int year){
    TFile *f7;

    if(year == 2016) f7 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_69200.root");
    else if(year == 2017) f7 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_69200.root");
    else if(year == 2018) f7 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_69200.root");

    TH1D *data_pileup = (TH1D *) f7->Get("pileup")->Clone();
    data_pileup->Scale(1./data_pileup->Integral());
    data_pileup->SetDirectory(0);
    pu_SF->data_pileup = data_pileup;
    pu_SF->pileup_ratio = (TH1D *) data_pileup->Clone("pileup_ratio");
    pu_SF->pileup_ratio->SetDirectory(0);
    f7->Close();

    if(pu_SF->data_pileup == NULL) printf("Something wrong getting Pileup SF!\n\n\n");
}

void setup_prefire_SFs(prefire_SFs *pre_SF, int year){

    if(year == 2016){

        TFile *f_el = TFile::Open("../analyze/SFs/2016/L1prefiring_photonpt_2016BtoH.root");
        TH2F *h1 = (TH2F *) f_el->Get("L1prefiring_photonpt_2016BtoH")->Clone();
        h1->SetDirectory(0);
        pre_SF->el_rate= h1;
        f_el->Close();

        TFile *f_jet = TFile::Open("../analyze/SFs/2016/L1prefiring_jetpt_2016BtoH.root");
        TH2F *h2 = (TH2F *) f_jet->Get("L1prefiring_jetpt_2016BtoH")->Clone();
        h2->SetDirectory(0);
        pre_SF->jet_rate= h2;
        f_jet->Close();
    }
    else if(year == 2017){

        TFile * f_el = TFile::Open("../analyze/SFs/2017/L1prefiring_photonpt_2017BtoF.root");
        TH2F *h1 = (TH2F *) f_el->Get("L1prefiring_photonpt_2017BtoF")->Clone();
        h1->SetDirectory(0);
        pre_SF->el_rate= h1;
        f_el->Close();

        TFile *f_jet = TFile::Open("../analyze/SFs/2017/L1prefiring_jetpt_2017BtoF.root");
        TH2F *h2 = (TH2F *) f_jet->Get("L1prefiring_jetpt_2017BtoF")->Clone();
        h2->SetDirectory(0);
        pre_SF->jet_rate= h2;
        f_jet->Close();
    }

}





void setup_mu_SFs(mu_SFs *era1, mu_SFs *era2, int year){
    TH1::AddDirectory(kFALSE);

    TFile *f1, *f2, *f3, *f4, *f5, *f6;


    if(year == 2016){
        f1 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_HLT.root");
        f1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("../analyze/SFs/2016/Mu_GH_HLT.root");
        f4->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF ->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF ->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        f2->Close();

        f3 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        f3->Close();


        f5 = TFile::Open("../analyze/SFs/2016/Mu_GH_ID.root");
        TH2D *ID_2 = (TH2D *) f5->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_2->SetDirectory(0);
        era2->ID_SF = ID_2;
        f5->Close();

        f6 = TFile::Open("../analyze/SFs/2016/Mu_GH_ISO.root");
        TH2D *ISO_2 = (TH2D *) f6->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_2->SetDirectory(0);
        era2->ISO_SF = ISO_2;
        f6->Close();

    }
    else if(year == 2017){
        printf("hlt\n");
        f1 = TFile::Open("../analyze/SFs/2017/Mu_HLT.root");
        f1->cd("IsoMu27_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();
        era2->HLT_SF = era1->HLT_SF;
        era2->HLT_MC_EFF = era1->HLT_MC_EFF;

        printf("id\n");
        f2 = TFile::Open("../analyze/SFs/2017/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        f2->Close();

        printf("iso\n");
        f3 = TFile::Open("../analyze/SFs/2017/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        f3->Close();


    }
    else if(year == 2018){
        f1 = TFile::Open("../analyze/SFs/2018/Mu_per1_HLT.root");
        f1->cd("IsoMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("../analyze/SFs/2018/Mu_per2_HLT.root");
        f4->cd("IsoMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("../analyze/SFs/2018/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        f2->Close();

        f3 = TFile::Open("../analyze/SFs/2018/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        f3->Close();


    }

    if(era1->HLT_SF == NULL || era1->ID_SF == NULL || era1->ISO_SF == NULL  ||
       era2->HLT_SF == NULL || era2->ID_SF == NULL || era2->ISO_SF == NULL) printf("Something wrong setup muon SF'S !!!! \n\n\n");
}

void setup_el_SF(el_SFs *sf, int year){
    //Setup electron SF's
    TFile *f_id, *f_reco, *f_hlt, *f_hlt2;
    if(year == 2016){
        f_hlt = TFile::Open("../analyze/SFs/2016/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2016/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2016/El_RECO.root");
    }
    else if(year == 2017){
        f_hlt = TFile::Open("../analyze/SFs/2017/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2017/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2017/El_RECO.root");
    }
    else if(year == 2018){
        f_hlt = TFile::Open("../analyze/SFs/2018/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2018/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2018/El_RECO.root");
    }

    TH2D *h1 = (TH2D *) f_id->Get("EGamma_SF2D")->Clone();
    h1->SetDirectory(0);
    sf->ID_SF = h1;
    f_id->Close();

    TH2D *h2 = (TH2D *) f_reco->Get("EGamma_SF2D")->Clone();
    h2->SetDirectory(0);
    sf->RECO_SF = h2;
    f_reco->Close();

    TH2D *h_hltsf1 = (TH2D *) f_hlt->Get("EGamma_SF2D")->Clone();
    h_hltsf1->SetDirectory(0);
    sf->HLT_SF = h_hltsf1;

    TH2D *h_hlt_mceff = (TH2D *) f_hlt->Get("EGamma_MCEff")->Clone();
    h_hlt_mceff->SetDirectory(0);
    sf->HLT_MC_EFF = h_hlt_mceff;

    
    if(sf->HLT_SF == NULL || sf->HLT_MC_EFF == NULL || sf->ID_SF == NULL || sf->RECO_SF == NULL ) printf("Something wrong setup electron SF'S !!!! \n\n\n");
}


void setup_pileup_systematic(pileup_systematics *pu_sys, int year){

    TFile *fin;
    if(year == 2016){
        fin = TFile::Open("../analyze/SFs/2016/pu_SF.root");
    }
    if(year == 2017){
        fin = TFile::Open("../analyze/SFs/2017/pu_SF.root");
    }
    if(year == 2018){
        fin = TFile::Open("../analyze/SFs/2018/pu_SF.root");
    }


    pu_sys->ratio_pileup_nom = (TH1D *) fin->Get("pu_ratio")->Clone();
    pu_sys->ratio_pileup_nom->SetDirectory(0);

    pu_sys->ratio_pileup_up = (TH1D *) fin->Get("pu_ratio_up")->Clone();
    pu_sys->ratio_pileup_up->SetDirectory(0);

    pu_sys->ratio_pileup_down = (TH1D *) fin->Get("pu_ratio_down")->Clone();
    pu_sys->ratio_pileup_down->SetDirectory(0);
}





void setup_pileup_systematic_old(pileup_systematics *pu_sys, int year){


    TFile *f7, *f8, *f9;
    if(year == 2016){
        f7 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_72383.root");
    }
    if(year == 2017){
        f7 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_72383.root");
    }
    if(year == 2018){
        f7 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_72383.root");
    }


    TH1D * pileup_data_nom = (TH1D *) f7->Get("pileup")->Clone();
    pileup_data_nom->Scale(1./pileup_data_nom->Integral());
    pileup_data_nom->SetDirectory(0);

    TH1D * pileup_data_up = (TH1D *) f8->Get("pileup")->Clone();
    pileup_data_up->Scale(1./pileup_data_up->Integral());
    pileup_data_up->SetDirectory(0);

    TH1D * pileup_data_down = (TH1D *) f9->Get("pileup")->Clone();
    pileup_data_down->Scale(1./pileup_data_down->Integral());
    pileup_data_down->SetDirectory(0);

    //Printf(" Ints are %.4e %.4e %.4e \n", pileup_data_nom->Integral(), pileup_data_up->Integral(), pileup_data_down->Integral());
    pu_sys->data_pileup_up = (TH1D *) pileup_data_up->Clone("data_pu_up");
    pu_sys->data_pileup_down = (TH1D *) pileup_data_down->Clone("data_pu_down");
    pu_sys->data_pileup_nom = (TH1D *) pileup_data_nom->Clone("data_pu_nom");
    pu_sys->data_pileup_up->SetDirectory(0);
    pu_sys->data_pileup_down->SetDirectory(0);
    pu_sys->data_pileup_nom->SetDirectory(0);

    pu_sys->ratio_pileup_up = (TH1D *) pileup_data_up->Clone("ratio_pu_up");
    pu_sys->ratio_pileup_down = (TH1D *) pileup_data_down->Clone("ratio_pu_down");
    pu_sys->ratio_pileup_nom = (TH1D *) pileup_data_nom->Clone("ratio_pu_nom");
    pu_sys->ratio_pileup_up->SetDirectory(0);
    pu_sys->ratio_pileup_down->SetDirectory(0);
    pu_sys->ratio_pileup_nom->SetDirectory(0);
    f7->Close();
    f8->Close();
    f9->Close();
}
#endif

