#ifndef DYAFB_TEMPMAKERC
#define DYAFB_TEMPMAKERC
#include "TempMaker.h"

TempMaker::TempMaker(TTree *t, bool isdata=false, int y  = 2016){
    is_data = isdata;
    year = y;
    t_in = t;
}

void TempMaker::setup(){



    nEntries  =  t_in->GetEntries();
    if(year == 2016) el_lumi=el_lumi16;
    if(year == 2017) el_lumi=el_lumi17;
    if(year == 2018) el_lumi=el_lumi18;


    t_in->SetBranchAddress("m", &m);
    t_in->SetBranchAddress("xF", &xF);
    t_in->SetBranchAddress("cost", &cost);
    t_in->SetBranchAddress("jet1_btag", &jet1_btag);
    t_in->SetBranchAddress("jet2_btag", &jet2_btag);
    t_in->SetBranchAddress("met_pt", &met_pt);
    t_in->SetBranchAddress("jet1_pt", &jet1_pt);
    t_in->SetBranchAddress("jet2_pt", &jet2_pt);
    t_in->SetBranchAddress("jet1_eta", &jet1_eta);
    t_in->SetBranchAddress("jet2_eta", &jet2_eta);
    t_in->SetBranchAddress("has_nobjets", &has_no_bjets);

    if(do_muons){
        t_in->SetBranchAddress("mu_p", &lep_p);
        t_in->SetBranchAddress("mu_m", &lep_m);
        t_in->SetBranchAddress("mu1_pt", &mu1_pt);
        t_in->SetBranchAddress("mu1_eta", &mu1_eta);
        t_in->SetBranchAddress("mu2_pt", &mu2_pt);
        t_in->SetBranchAddress("mu2_eta", &mu2_eta);
        t_in->SetBranchAddress("mu_p_SF", &mu_p_SF);
        t_in->SetBranchAddress("mu_m_SF", &mu_m_SF);
        t_in->SetBranchAddress("mu1_charge", &mu1_charge);
        t_in->SetBranchAddress("mu2_charge", &mu2_charge);
    }

    if(do_electrons){

        t_in->SetBranchAddress("el_p", &lep_p);
        t_in->SetBranchAddress("el_m", &lep_m);
        t_in->SetBranchAddress("el1_pt", &el1_pt);
        t_in->SetBranchAddress("el1_eta", &el1_eta);
        t_in->SetBranchAddress("el2_pt", &el2_pt);
        t_in->SetBranchAddress("el2_eta", &el2_eta);
        t_in->SetBranchAddress("el1_charge", &el1_charge);
        t_in->SetBranchAddress("el2_charge", &el2_charge);
    }
    if(do_emu){
        t_in->SetBranchAddress("el", &el);
        t_in->SetBranchAddress("mu", &mu);
        t_in->SetBranchAddress("mu1_charge", &mu1_charge);
        t_in->SetBranchAddress("met_pt", &met_pt);
        t_in->SetBranchAddress("el1_pt", &el1_pt);
        t_in->SetBranchAddress("el1_eta", &el1_eta);
        t_in->SetBranchAddress("mu1_pt", &mu1_pt);
        t_in->SetBranchAddress("mu1_eta", &mu1_eta);
    }

    if(is_one_iso){
        if(do_muons) t_in->SetBranchAddress("iso_mu", &iso_lep);
        if(do_electrons) t_in->SetBranchAddress("iso_el", &iso_lep);
        if(do_emu) t_in->SetBranchAddress("iso_lep", &iso_lep);

    }


    if(!is_data){
        if(is_gen_level){ 
            t_in->SetBranchAddress("inc_id1", &inc_id1);
            t_in->SetBranchAddress("inc_id2", &inc_id2);
            t_in->SetBranchAddress("cost_st", &cost_st);
            t_in->SetBranchAddress("pdf_weights", &pdf_weights);
            t_in->SetBranchAddress("nnpdf30_weight", &nnpdf30_weight);
        }
        t_in->SetBranchAddress("nJets", &nJets);
        t_in->SetBranchAddress("gen_weight", &gen_weight);
        t_in->SetBranchAddress("pu_SF", &pu_SF);
        t_in->SetBranchAddress("pu_SF_up", &pu_SF_up);
        t_in->SetBranchAddress("pu_SF_down", &pu_SF_down);
        t_in->SetBranchAddress("mu_R_up", &mu_R_up);
        t_in->SetBranchAddress("mu_R_down", &mu_R_down);
        t_in->SetBranchAddress("mu_F_up", &mu_F_up);
        t_in->SetBranchAddress("mu_F_down", &mu_F_down);
        t_in->SetBranchAddress("mu_RF_up", &mu_RF_up);
        t_in->SetBranchAddress("mu_RF_down", &mu_RF_down);
        t_in->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt);
        t_in->SetBranchAddress("jet1_btag_SF", &jet1_btag_SF);
        t_in->SetBranchAddress("jet2_btag_SF", &jet2_btag_SF);
        t_in->SetBranchAddress("jet1_flavour", &jet1_flavour);
        t_in->SetBranchAddress("jet2_flavour", &jet2_flavour);
        t_in->SetBranchAddress("alpha_up", &alphaS_up);
        t_in->SetBranchAddress("alpha_down", &alphaS_down);
        t_in->SetBranchAddress("prefire_SF", &prefire_SF);
        t_in->SetBranchAddress("prefire_SF_up", &prefire_SF_up);
        t_in->SetBranchAddress("prefire_SF_down", &prefire_SF_down);
        t_in->SetBranchAddress("top_ptrw", &top_ptrw);

        if(do_muons || do_emu){
            if(is_gen_level){
                t_in->SetBranchAddress("gen_mu_p", &gen_lep_p);
                t_in->SetBranchAddress("gen_mu_m", &gen_lep_m);
            }
            t_in->SetBranchAddress("era1_HLT_SF", &era1_HLT_SF);
            t_in->SetBranchAddress("era1_iso_SF", &era1_iso_SF);
            t_in->SetBranchAddress("era1_id_SF", &era1_id_SF);
            t_in->SetBranchAddress("era2_HLT_SF", &era2_HLT_SF);
            t_in->SetBranchAddress("era2_iso_SF", &era2_iso_SF);
            t_in->SetBranchAddress("era2_id_SF", &era2_id_SF);
            if(do_muons){
                t_in->SetBranchAddress("mu_p_SF_up", &mu_p_SF_up);
                t_in->SetBranchAddress("mu_m_SF_up", &mu_m_SF_up);
                t_in->SetBranchAddress("mu_p_SF_down", &mu_p_SF_down);
                t_in->SetBranchAddress("mu_m_SF_down", &mu_m_SF_down);
            }
        }

        if(do_electrons || do_emu){

            if(is_gen_level){
                t_in->SetBranchAddress("gen_el_p", &gen_lep_p);
                t_in->SetBranchAddress("gen_el_m", &gen_lep_m);
            }
            t_in->SetBranchAddress("el_id_SF", &el_id_SF);
            t_in->SetBranchAddress("el_reco_SF", &el_reco_SF);
            t_in->SetBranchAddress("el_HLT_SF", &el_HLT_SF);
            el_HLT_SF = 1.;
            if(do_elScale_sys || do_elSmear_sys){
                if(do_elScale_sys >0){
                    if(sys_label.find("elScaleStat") != string::npos){
                        t_in->SetBranchAddress("elp_scale_stat_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_stat_up", &elm_rescale);
                    }
                    if(sys_label.find("elScaleSyst") != string::npos){
                        t_in->SetBranchAddress("elp_scale_syst_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_syst_up", &elm_rescale);
                    }
                    if(sys_label.find("elScaleGain") != string::npos){
                        t_in->SetBranchAddress("elp_scale_gain_up", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_gain_up", &elm_rescale);
                    }
                }
                else if(do_elScale_sys < 0){
                    if(sys_label.find("elScaleStat") != string::npos){
                        t_in->SetBranchAddress("elp_scale_stat_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_stat_down", &elm_rescale);
                    }
                    if(sys_label.find("elScaleSyst") != string::npos){
                        t_in->SetBranchAddress("elp_scale_syst_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_syst_down", &elm_rescale);
                    }
                    if(sys_label.find("elScaleGain") != string::npos){
                        t_in->SetBranchAddress("elp_scale_gain_down", &elp_rescale);
                        t_in->SetBranchAddress("elm_scale_gain_down", &elm_rescale);
                    }
                }
                else if(do_elSmear_sys < 0){
                    t_in->SetBranchAddress("elp_smear_down", &elp_rescale);
                    t_in->SetBranchAddress("elm_smear_down", &elm_rescale);
                }
                else if(do_elSmear_sys > 0){
                    t_in->SetBranchAddress("elp_smear_up", &elp_rescale);
                    t_in->SetBranchAddress("elm_smear_up", &elm_rescale);
                }

            }

        }
    }
}


void TempMaker::setup_systematic(const string &s_label){
    sys_label = s_label;


    if(!sys_label.empty()){
        //doing some systematic


        if(sys_label.find("Up") != string::npos){
            sys_shift = 1;
        }
        else if(sys_label.find("Down") != string::npos){
            sys_shift = -1;
        }

        else{
            printf("systematic label not empty, but doesn't have Up or Down \n");
        }
    }

    if(sys_shift !=0){
        if(sys_label.find("BTAGCOR") != string::npos) do_btag_sys = sys_shift;
        else if(sys_label.find("BTAGUNCOR") != string::npos) do_btag_sys = 2*sys_shift;
        else if(sys_label.find("BTAGLIGHT") != string::npos) do_btag_sys = 3*sys_shift;
        else if(sys_label.find("Pu") != string::npos) do_pileup_sys = sys_shift;
        else if(sys_label.find("muRC") != string::npos) do_muRC_sys = sys_shift;
        else if(sys_label.find("muHLTBAR") != string::npos) do_muHLT_barrel_sys = sys_shift;
        else if(sys_label.find("muIDBAR") != string::npos) do_muID_barrel_sys = sys_shift;
        else if(sys_label.find("muISOBAR") != string::npos) do_muISO_barrel_sys = sys_shift;
        else if(sys_label.find("muHLTEND") != string::npos) do_muHLT_endcap_sys = sys_shift;
        else if(sys_label.find("muIDEND") != string::npos) do_muID_endcap_sys = sys_shift;
        else if(sys_label.find("muISOEND") != string::npos) do_muISO_endcap_sys = sys_shift;
        else if(sys_label.find("muIDSYS") != string::npos) do_muID_endcap_sys = sys_shift;
        else if(sys_label.find("muISOSYS") != string::npos) do_muISO_endcap_sys = sys_shift;

        else if(sys_label.find("elIDBAR") != string::npos) do_elID_barrel_sys = sys_shift;
        else if(sys_label.find("elHLTBAR") != string::npos) do_elHLT_barrel_sys = sys_shift;
        else if(sys_label.find("elRECOBAR") != string::npos) do_elRECO_barrel_sys = sys_shift;
        else if(sys_label.find("elIDEND") != string::npos) do_elID_endcap_sys = sys_shift;
        else if(sys_label.find("elHLTEND") != string::npos) do_elHLT_endcap_sys = sys_shift;
        else if(sys_label.find("elRECOEND") != string::npos) do_elRECO_endcap_sys = sys_shift;
        else if(sys_label.find("elIDSYS") != string::npos) do_elID_SYS_sys = sys_shift;
        else if(sys_label.find("elRECOSYS") != string::npos) do_elRECO_SYS_sys = sys_shift;

        else if(sys_label.find("elScaleStat") != string::npos) do_elScale_sys = sys_shift;
        else if(sys_label.find("elScaleSyst") != string::npos) do_elScale_sys = sys_shift;
        else if(sys_label.find("elScaleGain") != string::npos) do_elScale_sys = sys_shift;
        else if(sys_label.find("elSmear") != string::npos) do_elSmear_sys = sys_shift;





        else if(sys_label.find("prefire") != string::npos) do_prefire_sys = sys_shift;


        else if(sys_label.find("REFAC") != string::npos && sys_shift > 0) systematic = &mu_RF_up;
        else if(sys_label.find("REFAC") != string::npos && sys_shift < 0) systematic = &mu_RF_down;
        else if(sys_label.find("RENORM") != string::npos && sys_shift > 0) systematic = &mu_R_up;
        else if(sys_label.find("RENORM") != string::npos && sys_shift < 0) systematic = &mu_R_down;
        else if(sys_label.find("FAC") != string::npos && sys_shift > 0) systematic = &mu_F_up;
        else if(sys_label.find("FAC") != string::npos && sys_shift < 0) systematic = &mu_F_down;
        else if(sys_label.find("alphaS") != string::npos && sys_shift < 0) systematic = &alphaS_down;
        else if(sys_label.find("alphaS") != string::npos && sys_shift > 0) systematic = &alphaS_up;
        else if(sys_label.find("A0Den") != string::npos) do_A0_sys = sys_shift;

        else if(sys_label.find("MET") != string::npos){
            t_in->SetBranchAddress("met_pt", &dummy);
            if(sys_label.find("METJEC") != string::npos && sys_shift > 0) t_in->SetBranchAddress("met_jec_up", &met_pt);
            else if(sys_label.find("METJEC") != string::npos && sys_shift < 0) t_in->SetBranchAddress("met_jec_down", &met_pt);
            else if(sys_label.find("METJER") != string::npos && sys_shift > 0) t_in->SetBranchAddress("met_jer_up", &met_pt);
            else if(sys_label.find("METJER") != string::npos && sys_shift < 0) t_in->SetBranchAddress("met_jer_down", &met_pt);
            else if(sys_label.find("METHEM") != string::npos && sys_shift > 0) t_in->SetBranchAddress("met_hem_up", &met_pt);
            else if(sys_label.find("METHEM") != string::npos && sys_shift < 0) t_in->SetBranchAddress("met_hem_down", &met_pt);
            else printf("COULDN'T PARSE SYSTEMATIC %s !!! \n \n", sys_label.c_str());
        }



        else if(sys_label.find("pdf") != string::npos){
            if(sys_shift > 0) sscanf(sys_label.c_str(), "_pdf%iUp", &do_pdf_sys);
            else sscanf(sys_label.c_str(), "_pdf%iDown", &do_pdf_sys);
            printf("Doing pdf sys %i \n", do_pdf_sys);
        }
        else if(sys_label.find("ptrw") != string::npos){
            int foo;
            if(sys_shift > 0){
                sscanf(sys_label.c_str(), "_ptrw%ib%iUp", &do_ptrw_sys, &foo);
            }
            else{
                sscanf(sys_label.c_str(), "_ptrw%ib%iDown", &do_ptrw_sys, &foo);
                do_ptrw_sys *= -1;
            }
            printf("Doing ptrw sys %i \n", do_ptrw_sys);
        }
        else if(sys_label.find("emucostrw") != string::npos){
            int foo;
            if(do_emu_costrw){
                if(sys_shift > 0){
                    sscanf(sys_label.c_str(), "_emucostrw%ib%iUp", &do_emu_costrw_sys, &foo);
                }
                else{
                    sscanf(sys_label.c_str(), "_emucostrw%ib%iDown", &do_emu_costrw_sys, &foo);
                    do_emu_costrw_sys *= -1;
                }
                printf("Doing emu costrw sys %i \n", do_emu_costrw_sys);
            }
        }

        else printf("COULDN'T PARSE SYSTEMATIC %s !!! \n \n", sys_label.c_str());


        //extra parsing for lepton efficiency pt ranges
        if(sys_label.find("PTHIGH") != string::npos){
            printf("Doing high pt range for this systematic \n");
            el_SF_pt_range = 1;
        }
        else if(sys_label.find("PTLOW") != string::npos){
            printf("Doing low pt range for this systematic \n");
            el_SF_pt_range = -1;
        }


    }
}


void TempMaker::getEvent(int i){
    t_in->GetEntry(i);
    if(do_muons) not_cosmic = notCosmic(*lep_p, *lep_m);
    if(do_emu){
        cm = *el + *mu;
    }
    else cm = *lep_p + *lep_m;
    pt = cm.Pt();



}

void TempMaker::doCorrections(){
    bool did_something = false;
    if(do_muons){
        if(!do_RC || do_muRC_sys != 0){

            Float_t mu_p_corr, mu_m_corr;
            if(!do_RC){
                //undo RC 
                mu_p_corr = 1./mu_p_SF;
                mu_m_corr = 1./mu_m_SF;
            }
            else if(do_muRC_sys != 0){
                if(do_muRC_sys > 0){
                    mu_p_corr = mu_p_SF_up/mu_p_SF;
                    mu_m_corr = mu_m_SF_up/mu_m_SF;
                }
                else if(do_muRC_sys < 0){
                    mu_p_corr = mu_p_SF_down/mu_p_SF;
                    mu_m_corr = mu_m_SF_down/mu_p_SF;
                }
            }
            lep_p->SetPtEtaPhiE(lep_p->Pt() * mu_p_corr, lep_p->Eta(), lep_p->Phi(), mu_p_corr * lep_p->E());
            lep_m->SetPtEtaPhiE(lep_m->Pt() * mu_m_corr, lep_m->Eta(), lep_m->Phi(), mu_m_corr * lep_m->E());
            did_something = true;


        }
    }

    if(do_electrons){

        if(do_elScale_sys || do_elSmear_sys){
            lep_p->SetPtEtaPhiE(lep_p->Pt() * (elp_rescale), lep_p->Eta(), lep_p->Phi(), lep_p->E() *(elp_rescale));
            lep_m->SetPtEtaPhiE(lep_m->Pt() * (elm_rescale), lep_m->Eta(), lep_m->Phi(), lep_m->E() *(elm_rescale));
            did_something = true;
        }
    }
    if(did_something){
        cm = *lep_p + *lep_m;
        m = cm.M();
        xF = compute_xF(cm);
        cost = get_cost(*lep_p, *lep_m);
    }
    if(is_gen_level){
        gen_cm = *gen_lep_p + *gen_lep_m; 
        gen_m = gen_cm.M();
        gen_rap = gen_cm.Rapidity();
        gen_pt = gen_cm.Pt();
        gen_cost = cost_st;
    }


    if(std::isnan(cost_st)){
        printf("Gen level cost is Nan, using reco level \n");
        cost_st = cost;
    }

}

void TempMaker::fixRFNorm(TH2 *h, int mbin, int year){
    double avg = 1.;

    double *h_RF_up, *h_RF_down, *h_R_up, *h_R_down, *h_F_up, *h_F_down;


    if(sys_label.find("REFAC") != string::npos && sys_shift > 0) avg = RF_pdf_helper.h_RF_up->GetBinContent(mbin+1);
    else if(sys_label.find("REFAC") != string::npos && sys_shift < 0) avg = RF_pdf_helper.h_RF_down->GetBinContent(mbin+1);
    else if(sys_label.find("RENORM") != string::npos && sys_shift > 0) avg = RF_pdf_helper.h_R_up->GetBinContent(mbin+1);
    else if(sys_label.find("RENORM") != string::npos && sys_shift < 0) avg = RF_pdf_helper.h_R_down->GetBinContent(mbin+1);
    else if(sys_label.find("FAC") != string::npos && sys_shift > 0) avg = RF_pdf_helper.h_F_up->GetBinContent(mbin+1);
    else if(sys_label.find("FAC") != string::npos && sys_shift < 0) avg = RF_pdf_helper.h_F_down->GetBinContent(mbin+1);
    else if(sys_label.find("pdf") != string::npos && sys_shift > 0) avg = RF_pdf_helper.h_pdfs[do_pdf_sys-1]->GetBinContent(mbin+1);
    else if(sys_label.find("pdf") != string::npos && sys_shift < 0) avg = 1./RF_pdf_helper.h_pdfs[do_pdf_sys-1]->GetBinContent(mbin+1);

    if(avg != 1.){
        printf("Sys label was %s, mbin %i correcting average weight by %.5f \n", sys_label.c_str(), mbin, 1./avg);
    }
    h->Scale(1./avg);
}

float TempMaker::getEvtWeight(bool incl_btag_SFs = true){
    if(is_data){
        evt_weight = 1.;
        return 1.;
    }


    if(do_pileup_sys == -1) pu_SF = pu_SF_down;
    if(do_pileup_sys == 1) pu_SF = pu_SF_up;
    if(do_pdf_sys != 0){
        float pdf_weight = pdf_weights[do_pdf_sys-1];
        pdf_weight = std::max(std::min(pdf_weight, 3.0f), 0.333f);
        if(sys_shift > 0) *systematic = pdf_weight;
        if(sys_shift < 0) *systematic = 1./pdf_weight;
        //if(counter < 20) printf("%.4f \n", *systematic);
        //counter++;
    }
    if(std::isnan(*systematic)){
        //printf("sys is %.4f  \n", *systematic);
        *systematic = 1.;
    }
    *systematic = std::max(std::min(*systematic, 3.0f), 0.333f);
    //top_ptrw is 1 for non-ttbar samples
    double base_weight = gen_weight * (*systematic) * pu_SF * top_ptrw;
    if(incl_btag_SFs){
        if(do_btag_sys != 0){
#ifndef STAND_ALONE
            jet1_btag_SF = get_btag_weight(jet1_pt, jet1_eta, (Float_t) jet1_flavour , btag_effs, b_reader, do_btag_sys, btag_mc_eff_idx);
            jet2_btag_SF = get_btag_weight(jet2_pt, jet2_eta, (Float_t) jet2_flavour , btag_effs, b_reader, do_btag_sys, btag_mc_eff_idx);
#endif
        }
        if (nJets >= 1){
            base_weight *= jet1_btag_SF;
        }
        if (nJets >= 2){
            base_weight *= jet2_btag_SF;
        }
    }

    if(do_ptrw){
        base_weight *= get_ptrw_SF(ptrw_SFs, m, pt, do_ptrw_sys); 
    }

    if(year < 2018 && (do_electrons || do_emu)){
        if(do_prefire_sys == -1) base_weight *= prefire_SF_down;
        else if(do_prefire_sys == 1) base_weight *= prefire_SF_up;
        else base_weight *= prefire_SF;
    }

    if(do_emu_costrw){
        base_weight *= get_emu_costrw_SF(emu_costrw, cost, m, do_emu_costrw_sys);
    }


    if(do_muons){
        if(do_muHLT_barrel_sys || do_muHLT_endcap_sys){
            era1_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era1_SFs.HLT_SF, era1_SFs.HLT_MC_EFF, do_muHLT_barrel_sys, do_muHLT_endcap_sys);
            era2_HLT_SF = get_HLT_SF(mu1_pt, mu1_eta, mu2_pt, mu2_eta, era2_SFs.HLT_SF, era2_SFs.HLT_MC_EFF, do_muHLT_barrel_sys, do_muHLT_endcap_sys);
        }

        if(do_muISO_barrel_sys || do_muISO_endcap_sys){
            era1_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ISO_SF,  do_muISO_barrel_sys, do_muISO_endcap_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ISO_SF,  do_muISO_barrel_sys, do_muISO_endcap_sys);
            era2_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ISO_SF,  do_muISO_barrel_sys, do_muISO_endcap_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ISO_SF,  do_muISO_barrel_sys, do_muISO_endcap_sys);
        }
        else if(do_muISO_SYS_sys){
            era1_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ISO_SF_SYS,  do_muISO_SYS_sys, do_muISO_SYS_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ISO_SF,  do_muISO_SYS_sys, do_muISO_SYS_sys);
            era2_iso_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ISO_SF_SYS,  do_muISO_SYS_sys, do_muISO_SYS_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ISO_SF,  do_muISO_SYS_sys, do_muISO_SYS_sys);
        }


        if(do_muID_barrel_sys || do_muID_endcap_sys){
            era1_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ID_SF,  do_muID_barrel_sys, do_muID_endcap_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ID_SF,  do_muID_barrel_sys, do_muID_endcap_sys);
            era2_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ID_SF,  do_muID_barrel_sys, do_muID_endcap_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ID_SF,  do_muID_barrel_sys, do_muID_endcap_sys);
        }
        else if(do_muID_SYS_sys){
            era1_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era1_SFs.ID_SF_SYS,  do_muID_SYS_sys, do_muID_SYS_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era1_SFs.ID_SF,  do_muID_SYS_sys, do_muID_SYS_sys);
            era2_id_SF = get_mu_SF(mu1_pt, mu1_eta, year, era2_SFs.ID_SF_SYS,  do_muID_SYS_sys, do_muID_SYS_sys) * 
                get_mu_SF(mu2_pt, mu2_eta, year, era2_SFs.ID_SF,  do_muID_SYS_sys, do_muID_SYS_sys);
        
        }


        Float_t era1_weight = base_weight * era1_HLT_SF * era1_iso_SF * era1_id_SF ;
        Float_t era2_weight = base_weight * era2_HLT_SF * era2_iso_SF * era2_id_SF ;


        if(year==2016) evt_weight = 1000*(era1_weight*bcdef_lumi16 + era2_weight*gh_lumi16);
        if(year==2017) evt_weight = 1000*(era1_weight*mu_lumi17);
        if(year==2018) evt_weight = 1000*(era1_weight*mu_lumi18_era1 + era2_weight*mu_lumi18_era2);
    }
    else if(do_electrons){
        if(do_elID_barrel_sys || do_elID_endcap_sys) 
            el_id_SF = get_el_SF(el1_pt, el1_eta, el_SF.ID_SF, do_elID_barrel_sys, do_elID_endcap_sys, el_SF_pt_range) * 
                       get_el_SF(el2_pt, el2_eta, el_SF.ID_SF, do_elID_barrel_sys, do_elID_endcap_sys, el_SF_pt_range);
        if(do_elRECO_barrel_sys || do_elRECO_endcap_sys) 
            el_reco_SF = get_el_SF(el1_pt, el1_eta, el_SF.RECO_SF, do_elRECO_barrel_sys, do_elRECO_endcap_sys, el_SF_pt_range) * 
                         get_el_SF(el2_pt, el2_eta, el_SF.RECO_SF, do_elRECO_barrel_sys, do_elRECO_endcap_sys, el_SF_pt_range);
        if(do_elHLT_barrel_sys || do_elHLT_endcap_sys) 
            el_HLT_SF = get_HLT_SF(el1_pt, el1_eta, el2_pt, el2_eta, el_SF.HLT_SF,  el_SF.HLT_MC_EFF, do_elHLT_barrel_sys, do_elHLT_endcap_sys, el_SF_pt_range);




        evt_weight = 1000. * el_lumi * base_weight * el_id_SF *el_reco_SF * el_HLT_SF;;

    }

    else if(do_emu){
        //no systematics implemented for now
        base_weight *=  el_id_SF *el_reco_SF;

        Float_t era1_weight = base_weight * era1_HLT_SF * era1_iso_SF * era1_id_SF ;
        Float_t era2_weight = base_weight * era2_HLT_SF * era2_iso_SF * era2_id_SF ;


        if(year==2016) evt_weight = 1000*(era1_weight*bcdef_lumi16 + era2_weight*gh_lumi16);
        if(year==2017) evt_weight = 1000*(era1_weight*mu_lumi17);
        if(year==2018) evt_weight = 1000*(era1_weight*mu_lumi18_era1 + era2_weight*mu_lumi18_era2);


    }
    /*
       if(count < 100)
       printf("%.4f \n", evt_weight);
       count++;
       */
    return evt_weight;

}
float TempMaker::getLQReweightingDenom(int flag2 = 0){
    int flag1 = FLAG_MUONS;
    if(do_electrons) flag1 = FLAG_ELECTRONS;
    return get_LQ_reweighting_denom(LQ_helper, flag1, flag2, gen_m, gen_rap, cost_st);
}

float TempMaker::getReweightingDenom(){
    //printf("%.2f %.2f %.2f \n", cost, gen_cm.M(), gen_cm.Pt());
    return get_reweighting_denom(A0_helper, cost_st, gen_m, gen_pt, abs(gen_rap), do_A0_sys);
}



void TempMaker::finish(){
    t_in->ResetBranchAddresses();
}

TempMaker::~TempMaker(){}


void fakes_cost_reweight(TH2 *h_fakes, TH1 *h_rw, int systematic = 0){
    int n_var1_bins = h_fakes->GetNbinsX();
    int n_cost_bins = h_fakes->GetNbinsY();
    int n_rw_bins = h_rw->GetNbinsX();

    //printf("Systemtic is %i \n", systematic);

    if(n_cost_bins != n_rw_bins){
        printf("Fakes reweighting binning doesn't match! hist %i bins and rw_hist %i bins\n", n_cost_bins, n_rw_bins);
        exit(1);
    }
    for(int i=1; i<= h_fakes->GetNbinsX(); i++){
        for(int j=1; j<= h_fakes->GetNbinsY(); j++){
            float old_val = h_fakes->GetBinContent(i, j);
            float old_err = h_fakes->GetBinError(i, j);

            float correction = h_rw->GetBinContent(j);

            //one systematic for every |cos(theta)| bin (nbins / 2)
            if(systematic != 0){

                //error includes sys and stat errors
                float error = h_rw->GetBinError(j);

                int sys_bin = abs(systematic);
                int opp_bin = (h_rw->GetNbinsX() + 1) -sys_bin;
                if(j == sys_bin || j == opp_bin){
                    //shift the reweighting in this bin by the error
                    if(systematic >0) correction += error;
                    if(systematic <0) correction -= error;
                }

            }
            float new_val = old_val * correction;
            new_val = std::max(new_val, 1e-6f);

            float new_err = (old_err/old_val) * new_val; // stat error


            //printf("Old %.1f +/- %.1f, correction %.2f New %.1f +/- %.1f \n", old_val, old_err, correction, new_val, new_err);

            h_fakes->SetBinContent(i,j, new_val);
            h_fakes->SetBinError(i,j, new_err);
        }
    }
}

#endif
