#include "../../utils/NTupleReader.C"




void ElEl_reco_mc_batch(int nJobs =1, int iJob = 0, string fin = "", int year=-1)
{


    if(fin == "") fin = string("EOS_files/2017/DY_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_dy_test.root", false);
    if (year == -1) nt.year = 2017;


    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;
    nt.do_SFs = true;
    nt.setupSFs();
    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_DY_back");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");

    int iso_el;
    nt.outTrees[2]->Branch("iso_el", &iso_el); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.dielec_id && nt.cm_m > 130. ){
                nt.fillEvent();
                nt.fillEventSFs();

                bool one_iso = nt.el_iso0 ^ nt.el_iso1;


                //pick the category
                if(nt.opp_sign && nt.el_iso0 && nt.el_iso1){ //signal region
                    nt.parseGenParts();
                    if(nt.signal_event && !nt.failed_match){
                        nt.nSignal++;
                        nt.outTrees[0]->Fill();
                    }
                    else{
                        nt.outTrees[1]->Fill();
                    }
                }
                else if(!nt.opp_sign && nt.el_iso0 && nt.el_iso1){ //samesign region
                    nt.outTrees[4]->Fill();
                }
                else if(one_iso){ //wjets control region
                    if(nt.el_iso0) iso_el = 0;
                    else           iso_el = 1;
                    nt.outTrees[2]->Fill();
                }
                else if(!nt.el_iso0 && !nt.el_iso1){ //qcd control region
                    nt.outTrees[3]->Fill();
                }


            }
        } 

        printf("moving on to next file, currently %i events %i Taus %i fails \n\n", nt.nEvents, nt.nTauTau, nt.nFailedID);


    }
    printf("There were %i qqbar, %i qGlu (%i of them tautau) in %i kept events in %i files. "
            "There were also %i background events (%i qq and %i gg). "
            "There were %i Failed ID's \n" , 
            nt.nQQb, nt.nQGlu, nt.nTauTau, nt.nSignal, nt.fileCount, nt.nQQ + nt.nGluGlu, nt.nQQ, nt.nGluGlu, nt.nFailedID);
    //printf("Ran on MC data and produced templates with %i events\n", nEvents);
    nt.finish();

    return;
}

