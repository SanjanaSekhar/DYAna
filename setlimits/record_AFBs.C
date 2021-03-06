#include "madgraph_lhe_reader.C"
#include "../utils/bins.h"


void record_AFBs(int M_Zp = 2500, double cpl = 1.0){
    FILE *fout = fopen("AFBs.txt", "a");
    
    Double_t AFB_Zp;
    //printf("AFBs are: ");
    printf("Writing out AFBs for M %i kL %.2f \n", M_Zp, cpl);
    fprintf(fout, "M_Zp=%i kL=%.2f: ", M_Zp, cpl);
    int bin_start = 5;
    for(int idx=1; idx<4; idx++){

        char fin[120];
        sprintf(fin, "MG5_aMC_v2_6_2/Zp_to_mumu_temp/Events/run_0%i/unweighted_events.lhe", idx);
        int bin = idx + bin_start;

        AFB_Zp = get_AFB(fin, m_bins[bin-1], m_bins[bin], -1, true);
        //if can't get a good AFB, don't record
        if(AFB_Zp > 1.0 || AFB_Zp <= 0.01) exit(1);
        fprintf(fout, "%.3f ", AFB_Zp);
        printf("%.3f ", AFB_Zp);


    }
    printf("\n");
    fprintf(fout, "\n");


               

}
   




