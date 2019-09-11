#include "../../utils/NTupleReader.C"




void ElEl_reco_data_batch(int nJobs =1, int iJob = 0, string fin = "", int year=-1)
{

    if (fin == "") fin = string("EOS_files/2016/SingleElectron_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_data_test.root", true);
    if (year == -1) nt.year = 2016;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");

    int iso_el;
    nt.outTrees[1]->Branch("iso_el", &iso_el); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);

            if(nt.good_trigger && nt.dielec_id && nt.cm_m > 130.){
                nt.fillEvent();
                bool one_iso = nt.el_iso0 ^ nt.el_iso1;

                //pick the category
                if(nt.opp_sign && nt.el_iso0 && nt.el_iso1){ //signal region
                    nt.outTrees[0]->Fill();
                }
                else if(!nt.opp_sign && nt.el_iso0 && nt.el_iso1){ //samesign region
                    nt.outTrees[3]->Fill();
                }
                else if(one_iso){ //wjets control region
                    if(nt.el_iso0) iso_el = 0;
                    else           iso_el = 1;
                    nt.outTrees[1]->Fill();
                }
                else if(!nt.el_iso0 && !nt.el_iso1){ //qcd control region
                    nt.outTrees[2]->Fill();
                }
            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    printf("Finished. There were %i events from %i files \n\n", nt.nEvents, nt.fileCount);
    nt.finish();

    return;
}

