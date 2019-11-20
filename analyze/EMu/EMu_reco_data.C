#include "../../utils/NTupleReader.C"




void EMu_reco_data(int nJobs =1, int iJob = 0, string fin = "", int year = -1)
{


    if (fin == "") fin = string("EOS_files/2017/SingleMuon_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/EMu_data_test.root", true);
    if (year == -1) year = 2017;
    nt.year = year;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_emu = true;

    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");
    int iso_lep;
    nt.outTrees[1]->Branch("iso_lep", &iso_lep);


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.emu_ids  && nt.cm_m > 130.){
                nt.fillEvent();

                bool one_iso = nt.mu_tight_id0 ^ nt.el_iso0;

                //pick the category
                if(nt.opp_sign && nt.mu_tight_id0 && nt.el_iso0){ //signal region
                    nt.outTrees[0]->Fill();
                }
                else if(!nt.opp_sign && nt.mu_tight_id0 && nt.el_iso0){ //samesign region
                    nt.outTrees[3]->Fill();
                }
                else if(nt.opp_sign && one_iso){ //wjets control region
                    if(nt.mu_tight_id0) iso_lep = 0;
                    else          iso_lep = 1;
                    nt.outTrees[1]->Fill();
                }
                else if(nt.opp_sign && !nt.mu_tight_id0 && !nt.el_iso0){ //qcd control region
                    nt.outTrees[2]->Fill();
                }
            }


        }

        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    nt.finish();

    return;
}

