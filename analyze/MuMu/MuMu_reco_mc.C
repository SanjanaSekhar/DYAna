#include "../../utils/NTupleReader.C"




void MuMu_reco_mc(int nJobs =1, int iJob = 0, string fin = "", int year =-1)
{
    if(fin == "") fin = string("EOS_files/2016/DY_files_test.txt");
    NTupleReader nt(fin.c_str(),"output_files/MuMu_DY_test16.root", false);
    if (year == -1) year = 2016;
    nt.year = year;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = true;
    nt.do_RC = true;
    nt.RC_from_gen = true;
    nt.setupSFs();
    nt.setupRC();
    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");
    nt.setupOutputTree("T_DY_back");


    int iso_mu;
    nt.outTrees[1]->Branch("iso_mu", &iso_mu); 

    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.loose_dimuon_id &&
                    //nt.cm_m > 130. ){
                    nt.cm_m > 50. && nt.cm_m < 130. ){
                nt.fillEvent();
                nt.fillEventSFs();
                nt.parseGenParts();
                nt.fillEventRC();
                bool one_tight = nt.mu_tight_id0 ^ nt.mu_tight_id1;

                //pick the category
                if(nt.opp_sign && nt.tight_dimuon_id){ //signal region
                    if(nt.signal_event && !nt.failed_match){
                        nt.nSignal++;
                        nt.outTrees[0]->Fill();
                    }
                    else{
                        nt.outTrees[4]->Fill();
                    }
                }
                else if(!nt.opp_sign && nt.tight_dimuon_id){ //samesign region
                    nt.outTrees[3]->Fill();
                }
                else if(one_tight){ //wjets control region
                    if(nt.mu_tight_id0) iso_mu = 0;
                    else           iso_mu = 1;
                    nt.outTrees[1]->Fill();
                }
                else if(!nt.mu_tight_id0 && !nt.mu_tight_id1){ //qcd control region
                    nt.outTrees[2]->Fill();
                }

            }
        } 

        printf("moving on to next file, currently %i events %i Taus %i fails \n\n", nt.nEvents, nt.nTauTau, nt.nFailedID);


    }
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files."
            "There were also %i background events (%i qq and %i gg)"
            "There were %i Failed ID's \n" , 
            nt.nQQb, nt.nQGlu, nt.nTauTau, nt.nSignal, nt.fileCount, nt.nQQ + nt.nGluGlu, nt.nQQ, nt.nGluGlu, nt.nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    nt.finish();

    return;
}

