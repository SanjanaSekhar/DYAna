#include "../../utils/NTupleReader.C"



void Select_gen_level(int nJobs =1, int iJob = 0, string fin = "", int year =-1)
{
    if(fin == "") fin = string("EOS_files/2017/DY_files_noskim.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu17_dy_gen.root", false);
    if (year == -1) year = 2017;
    nt.year = year;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.RC_from_gen = true;

    nt.print_gen_warning = false;

    int good = 0;
    int bad = 0;
    int nEl(0), nMu(0), nTau(0);

    nt.fout->cd();
    TTree *t_el = new TTree("T_gen_el", "");
    TTree *t_mu = new TTree("T_gen_mu", "");
    TTree *t_tau = new TTree("T_gen_tau", "");

    TLorentzVector gen_p, gen_m, cm;
    float cost_st,cost,m, gen_weight;

    float mu_R_up, mu_R_down, mu_F_up, mu_F_down, mu_RF_up, mu_RF_down;
   

    t_el->Branch("gen_p", "TLorentzVector", &nt.gen_lep_p_vec);
    t_el->Branch("gen_m", "TLorentzVector", &nt.gen_lep_m_vec);
    t_el->Branch("hard_p", "TLorentzVector", &nt.hard_lep_p_vec);
    t_el->Branch("hard_m", "TLorentzVector", &nt.hard_lep_m_vec);
    t_el->Branch("inc1", "TLorentzVector", &nt.inc1_vec);
    t_el->Branch("inc2", "TLorentzVector", &nt.inc2_vec);
    t_el->Branch("gen_weight", &gen_weight);
    t_el->Branch("m", &m);
    t_el->Branch("cost_st", &cost_st);
    t_el->Branch("cost", &cost);
    t_el->Branch("mu_R_up", &mu_R_up);
    t_el->Branch("mu_R_down", &mu_R_down);
    t_el->Branch("mu_F_up", &mu_F_up);
    t_el->Branch("mu_F_down", &mu_F_down);
    t_el->Branch("mu_RF_up", &mu_RF_up);
    t_el->Branch("mu_RF_down", &mu_RF_down);
    t_el->Branch("inc_id1", &nt.inc_id1);
    t_el->Branch("inc_id2", &nt.inc_id2);
    t_el->Branch("sig_event", &nt.signal_event);
    

    t_mu->Branch("gen_p", "TLorentzVector", &nt.gen_lep_p_vec);
    t_mu->Branch("gen_m", "TLorentzVector", &nt.gen_lep_m_vec);
    t_mu->Branch("hard_p", "TLorentzVector", &nt.hard_lep_p_vec);
    t_mu->Branch("hard_m", "TLorentzVector", &nt.hard_lep_m_vec);
    t_mu->Branch("inc1", "TLorentzVector", &nt.inc1_vec);
    t_mu->Branch("inc2", "TLorentzVector", &nt.inc2_vec);
    t_mu->Branch("gen_weight", &gen_weight);
    t_mu->Branch("m", &m);
    t_mu->Branch("cost_st", &cost_st);
    t_mu->Branch("cost", &cost);
    t_mu->Branch("mu_R_up", &mu_R_up);
    t_mu->Branch("mu_R_down", &mu_R_down);
    t_mu->Branch("mu_F_up", &mu_F_up);
    t_mu->Branch("mu_F_down", &mu_F_down);
    t_mu->Branch("mu_RF_up", &mu_RF_up);
    t_mu->Branch("mu_RF_down", &mu_RF_down);
    t_mu->Branch("inc_id1", &nt.inc_id1);
    t_mu->Branch("inc_id2", &nt.inc_id2);
    t_mu->Branch("sig_event", &nt.signal_event);

    t_tau->Branch("gen_p", "TLorentzVector", &nt.gen_lep_p_vec);
    t_tau->Branch("gen_m", "TLorentzVector", &nt.gen_lep_m_vec);
    t_tau->Branch("hard_p", "TLorentzVector", &nt.hard_lep_p_vec);
    t_tau->Branch("hard_m", "TLorentzVector", &nt.hard_lep_m_vec);
    t_tau->Branch("inc1", "TLorentzVector", &nt.inc1_vec);
    t_tau->Branch("inc2", "TLorentzVector", &nt.inc2_vec);
    t_tau->Branch("gen_weight", &gen_weight);
    t_tau->Branch("m", &m);
    t_tau->Branch("cost_st", &cost_st);
    t_tau->Branch("cost", &cost);
    t_tau->Branch("inc_id1", &nt.inc_id1);
    t_tau->Branch("inc_id2", &nt.inc_id2);

    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            nt.fillEvent();


            mu_F_up = nt.scale_Weights[0];
            mu_F_down = nt.scale_Weights[1];
            mu_R_up = nt.scale_Weights[2];
            mu_R_down = nt.scale_Weights[4];
            mu_RF_up = nt.scale_Weights[3];
            mu_RF_down = nt.scale_Weights[5];

            int gen_id = nt.selectAnyGenParts();
            gen_p = nt.gen_lep_p_vec;
            gen_m = nt.gen_lep_m_vec;
            cm = gen_p + gen_m;
            m = cm.M();
            if(m > 150.){
                cost_st = nt.cost_st;
                cost = get_cost(gen_p, gen_m);
                gen_weight = nt.gen_weight;
                nt.inc1_vec.Print();
                nt.inc2_vec.Print();


                if(gen_id == 11){
                    nEl++;
                    t_el->Fill();
                }
                else if(gen_id == 13){
                    nMu++;
                    t_mu->Fill();
                }
                else if(gen_id == 15){
                    nTau++;
                    t_tau->Fill();
                }

            }
        } 

        printf("moving on to next file, currently %i events %i Els %i Mus %i Taus %i fails \n\n", nt.nEvents,nEl, nMu, nTau, nt.nFailedID );


    }
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files."
            "There were also %i background events (%i qq and %i gg)"
            "There were %i Failed ID's \n" , 
            nt.nQQb, nt.nQGlu, nt.nTauTau, nt.nSignal, nt.fileCount, nt.nQQ + nt.nGluGlu, nt.nQQ, nt.nGluGlu, nt.nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    nt.fout->cd();

    t_el->Write();
    t_mu->Write();
    t_tau->Write();

    nt.finish();

    return;
}

