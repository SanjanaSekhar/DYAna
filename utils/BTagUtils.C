#ifndef STAND_ALONE
#include "TH2D.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

typedef struct {
    BTagCalibrationReader b_reader;
    BTagCalibrationReader c_reader;
    BTagCalibrationReader udsg_reader;

} BTag_readers;

typedef struct {
    TH2D *b_eff_dy;
    TH2D *c_eff_dy;
    TH2D *udsg_eff_dy;
    TH2D *b_eff_ttbar;
    TH2D *c_eff_ttbar;
    TH2D *udsg_eff_ttbar;
    TH2D *b_eff_diboson;
    TH2D *c_eff_diboson;
    TH2D *udsg_eff_diboson;
} BTag_effs;


Double_t btag_eff(Double_t pt, Double_t eta,TH2D *mc_eff){
    if (pt >=500) pt = 450;

    TAxis* x_ax =  mc_eff->GetXaxis();
    TAxis *y_ax =  mc_eff->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t eff = mc_eff->GetBinContent(xbin, ybin);
    if(eff < 1e-8){
        eff = 1e-8;
        //printf("Warning: 0 efficiency for pt %.0f, eta %1.1f! \n", pt, eta);
    }
    return eff;
}

Double_t btag_weight_helper(Double_t pt, Double_t eta, Double_t SF, TH2D *mc_eff){
    if (pt >=500) pt = 450;

    TAxis* x_ax =  mc_eff->GetXaxis();
    TAxis *y_ax =  mc_eff->GetYaxis();
    int xbin = x_ax->FindBin(pt);
    int ybin = y_ax->FindBin(std::abs(eta));

    Double_t eff = mc_eff->GetBinContent(xbin, ybin);
    if(eff < 1e-8){
        eff = 1e-8;
        //printf("Warning: 0 efficiency for pt %.0f, eta %1.1f! \n", pt, eta);
    }
    //printf("Efficiency is %f \n", eff);
    Double_t weight = (1-SF*eff)/(1-eff);
    return weight;
}
Double_t get_btag_weight(Double_t pt, Double_t eta, Float_t flavour, BTag_effs btag_effs, BTag_readers b_readers, int systematic = 0, int btag_mc_eff_idx = 0){
    //compute weighting from btagging scale factors
    Double_t weight, bjet_SF;

    char const *sys;

    TH2D *b_eff, *c_eff, *udsg_eff;

    if(btag_mc_eff_idx == 0){ //tbar
        b_eff = btag_effs.b_eff_ttbar;
        c_eff = btag_effs.c_eff_ttbar;
        udsg_eff = btag_effs.udsg_eff_ttbar;
    }
    else if(btag_mc_eff_idx == 1){ //dy
        b_eff = btag_effs.b_eff_dy;
        c_eff = btag_effs.c_eff_dy;
        udsg_eff = btag_effs.udsg_eff_dy;
    }
    else if(btag_mc_eff_idx == 2){ //dy
        b_eff = btag_effs.b_eff_diboson;
        c_eff = btag_effs.c_eff_diboson;
        udsg_eff = btag_effs.udsg_eff_diboson;
    }

    if(std::abs(flavour - 5.) < 0.01){ //bjet

        if (systematic == 0 || systematic == 3 || systematic  == -3) sys = "central";
        if (systematic == 1) sys = "up_correlated";
        if (systematic == -1) sys = "down_correlated";
        if (systematic == 2) sys = "up_uncorrelated";
        if (systematic == -2) sys = "down_uncorrelated";

        bjet_SF = b_readers.b_reader.eval_auto_bounds(sys, BTagEntry::FLAV_B, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour - 4.) < 0.01){ //cjet

        if (systematic == 0 || systematic == 3 || systematic  == -3) sys = "central";
        if (systematic == 1) sys = "up_correlated";
        if (systematic == -1) sys = "down_correlated";
        if (systematic == 2) sys = "up_uncorrelated";
        if (systematic == -2) sys = "down_uncorrelated";

        bjet_SF = b_readers.c_reader.eval_auto_bounds(sys, BTagEntry::FLAV_C, eta, pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet

        //udsgs only have one uncertainty which is taken to be uncorrelated wrt
        //to others and the year
        if (systematic == 3) sys = "up";
        if (systematic == -3) sys = "down";
        if(abs(systematic) < 3) sys = "central";


        bjet_SF = b_readers.udsg_reader.eval_auto_bounds(sys, BTagEntry::FLAV_UDSG, eta,pt);
        weight = btag_weight_helper(pt, eta, bjet_SF, udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    if(bjet_SF == 0) printf("WARNING: Scale factor return 0 for Flavour %1.0f pt %.0f eta %1.1f \n!",
                            flavour, pt, eta);
    return weight;
}


float get_pos_btag_weight(int nJets, bool jet1_tagged, Double_t pt1, Double_t eta1, Float_t flavour1, 
        bool jet2_tagged, Double_t pt2, Double_t eta2, Float_t flavour2, BTag_effs btag_effs, BTag_readers b_readers){
    //compute weighting from btagging scale factors for tagging one of the jets
    //as a b

    Double_t bjet1_SF, bjet1_eff, bjet2_SF, bjet2_eff;
    float jet1_weight, jet2_weight;
    TH2D *b_eff = btag_effs.b_eff_ttbar;
    TH2D *c_eff = btag_effs.c_eff_ttbar;
    TH2D *udsg_eff = btag_effs.udsg_eff_ttbar;

    if(nJets ==0) return 1.;

    if(std::abs(flavour1 - 5.) < 0.01){ //bjet
        bjet1_SF = b_readers.b_reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta1, pt1);
        bjet1_eff= btag_eff(pt1, eta1, b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour1 - 4.) < 0.01){ //cjet
        bjet1_SF = b_readers.c_reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta1, pt1);
        bjet1_eff= btag_eff(pt1, eta1, c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet1_SF = b_readers.udsg_reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta1,pt1);
        bjet1_eff= btag_eff(pt1, eta1, udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }

    if(jet1_tagged) jet1_weight = bjet1_SF;
    else  jet1_weight = (1. - bjet1_eff * bjet1_SF) / (1. - bjet1_eff);

    if(nJets ==1) return jet1_weight;
    

    if(std::abs(flavour2 - 5.) < 0.01){ //bjet
        bjet2_SF = b_readers.b_reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta2, pt2);
        bjet2_eff= btag_eff(pt2, eta2, b_eff);
        //printf("B jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);

    }
    else if (std::abs(flavour2 - 4.) < 0.01){ //cjet
        bjet2_SF = b_readers.c_reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta2, pt2);
        bjet2_eff= btag_eff(pt2, eta2, c_eff);
        //printf("C jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }
    else{ //udsg jet
        bjet2_SF = b_readers.udsg_reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta2,pt2);
        bjet2_eff= btag_eff(pt2, eta2, udsg_eff);
        //printf("UDSG jet, SF is %0.3f, weight is %.4f \n", bjet_SF, weight);
    }

    if(jet2_tagged) jet2_weight = bjet2_SF;
    else  jet2_weight = (1. - bjet2_eff * bjet2_SF) / (1. - bjet2_eff);
    
    float result  = jet1_weight * jet2_weight;
    
    printf("tagged1  flav1 tagged weight %i  %.0f %.2f, tagged2, flav2, weight  %i %.0f %.2f, res %.2f \n", jet1_tagged, flavour1, jet1_weight, jet2_tagged, flavour2, jet2_weight, result);
    return result;
}

void setup_btag_SFs(BTag_readers *btag_r, BTag_effs *b_effs, int year, bool incl_systematics = true){
    printf("Setting up btag SFs... \n");
    TH1::AddDirectory(kFALSE);
    char file[80];
    BTagCalibration calib_bc, calib_udsg;
    
    if(year == 2016){
        calib_bc = BTagCalibration("DeepCSV", "../analyze/SFs/2016/DeepCSV_2016LegacySF_V1_YearCorrelation-V1.csv");
        calib_udsg = BTagCalibration("DeepCSV", "../analyze/SFs/2016/DeepCSV_2016LegacySF_V1.csv");
    }
    if(year == 2017){
        calib_bc = BTagCalibration("DeepCSV", "../analyze/SFs/2017/DeepCSV_94XSF_V4_B_F_YearCorrelation-V1.csv");
        calib_udsg = BTagCalibration("DeepCSV", "../analyze/SFs/2017/DeepCSV_94XSF_V4_B_F.csv");
    }
    if(year == 2018){
        calib_bc = BTagCalibration("DeepCSV", "../analyze/SFs/2018/DeepCSV_102XSF_V1_YearCorrelation-V1.csv");
        calib_udsg = BTagCalibration("DeepCSV", "../analyze/SFs/2018/DeepCSV_102XSF_V1.csv");
    }

    std::vector<string> sys_types, udsg_sys_types;
    if(incl_systematics){ 
        sys_types = {"up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"};
        udsg_sys_types = {"up", "down"};
    }

    btag_r->b_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", sys_types);
    btag_r->b_reader.load(calib_bc, BTagEntry::FLAV_B, "comb");
    btag_r->c_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", sys_types);
    btag_r->c_reader.load(calib_bc, BTagEntry::FLAV_C, "comb");
    btag_r->udsg_reader = BTagCalibrationReader (BTagEntry::OP_MEDIUM, "central", udsg_sys_types);
    btag_r->udsg_reader.load(calib_udsg, BTagEntry::FLAV_UDSG, "incl");
    
    TFile *f0_ttbar, *f0_dy, *f0_diboson;

    TFile *f0;
    if (year == 2016){
        f0_ttbar = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_ttbar.root");
        f0_dy = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_dy.root");
        f0_diboson = TFile::Open("../analyze/SFs/2016/Btag_eff_MC_2016_diboson.root");
    }
    else if (year == 2017){
        f0_ttbar = TFile::Open("../analyze/SFs/2017/Btag_eff_MC_2017_ttbar.root");
        f0_dy = TFile::Open("../analyze/SFs/2017/Btag_eff_MC_2017_dy.root");
        f0_diboson = TFile::Open("../analyze/SFs/2017/Btag_eff_MC_2017_diboson.root");
    }
    else if (year == 2018){
        f0_ttbar = TFile::Open("../analyze/SFs/2018/Btag_eff_MC_2018_ttbar.root");
        f0_dy = TFile::Open("../analyze/SFs/2018/Btag_eff_MC_2018_dy.root");
        f0_diboson = TFile::Open("../analyze/SFs/2018/Btag_eff_MC_2018_diboson.root");
    }

    f0_ttbar->cd();
    b_effs->b_eff_ttbar = (TH2D *) f0_ttbar->Get("b_eff")->Clone();
    b_effs->c_eff_ttbar = (TH2D *) f0_ttbar->Get("c_eff")->Clone();
    b_effs->udsg_eff_ttbar = (TH2D *) f0_ttbar->Get("udsg_eff")->Clone();
    b_effs->b_eff_ttbar->SetDirectory(0);
    b_effs->c_eff_ttbar->SetDirectory(0);
    b_effs->udsg_eff_ttbar->SetDirectory(0);


    f0_dy->cd();
    b_effs->b_eff_dy = (TH2D *) f0_dy->Get("b_eff")->Clone();
    b_effs->c_eff_dy = (TH2D *) f0_dy->Get("c_eff")->Clone();
    b_effs->udsg_eff_dy = (TH2D *) f0_dy->Get("udsg_eff")->Clone();
    b_effs->b_eff_dy->SetDirectory(0);
    b_effs->c_eff_dy->SetDirectory(0);
    b_effs->udsg_eff_dy->SetDirectory(0);


    f0_diboson->cd();
    b_effs->b_eff_diboson = (TH2D *) f0_diboson->Get("b_eff")->Clone();
    b_effs->c_eff_diboson = (TH2D *) f0_diboson->Get("c_eff")->Clone();
    b_effs->udsg_eff_diboson = (TH2D *) f0_diboson->Get("udsg_eff")->Clone();
    b_effs->b_eff_diboson->SetDirectory(0);
    b_effs->c_eff_diboson->SetDirectory(0);
    b_effs->udsg_eff_diboson->SetDirectory(0);


    printf(" Done \n");
}
#endif
