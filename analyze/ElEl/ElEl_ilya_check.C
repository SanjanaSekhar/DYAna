#include "../../utils/NTupleReader.C"




void ElEl_ilya_check(int nJobs =1, int iJob = 0, string fin = "", int year=-1)
{

    if (fin == "") fin = string("EOS_files/2016/SingleElectron_files_test.txt");
    NTupleReader nt(fin.c_str(),"output_files/ElEl_data_test.root", true);
    if (year == -1) year = 2016;
    nt.year = year;

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_electrons = true;


    std::set<ULong64_t>  evts = {297058495, 297202109, 296572658}; 

    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(evts.find(nt.evt_EventNumber) != evts.end()){

            //if(nt.HLT_El27 && nt.el_iso0 && nt.el_iso1 && 
                    //nt.el_Pt[0] > 30. && nt.el_Pt[1] > 30. && nt.cm_m > 50.){
                nt.fillEvent();
                printf("Run %i Lumiblock %i Event Num %llu \n", nt.evt_RunNumber, nt.evt_LumiBlock, nt.evt_EventNumber);
                printf("Trigger %i \n", nt.HLT_El27);

                printf("Electron 1 TightID %i : pt eta phi: %.2f %.2f %.2f \n", nt.el_iso0, nt.el_Pt[0], nt.el_Eta[0], nt.el_Phi[0]);
                printf("Electron 2 TightID %i : pt eta phi: %.2f %.2f %.2f \n", nt.el_iso1, nt.el_Pt[1], nt.el_Eta[1], nt.el_Phi[1]);
                printf("Dielectron pt rap phi m: %.2f %.2f %.2f %.2f \n", nt.cm.Pt(), nt.cm.Rapidity(), nt.cm.Phi(), nt.cm.M());
                //if(nt.nEvents > 100) return;
            }
        } 


        printf("moving on to next file, currently %i events \n\n", nt.nEvents);

    }
    printf("Finished. There were %i events from %i files \n\n", nt.nEvents, nt.fileCount);
    nt.finish();

    return;
}

