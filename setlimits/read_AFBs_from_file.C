#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>


const int n_bins = 3;
Double_t AFB_SM[n_bins] = {0.570, 0.575, 0.590}; //same generator as Z'
Double_t AFB_exp[n_bins] = {0.608, 0.610, 0.612}; //amc@NLO values
Double_t AFB_measured[n_bins] = {0.617, 0.595, 0.639};
//Double_t AFB_unc[n_bins] = {0.022, 0.033, 0.054}; //measurement uncertainties
Double_t AFB_unc[n_bins] = {0.029, 0.04, 0.058}; //Total uncertainty, adding 0.02 theory uncertainty

void read_AFBs(FILE *f1, Double_t *AFBs, int M_Zp, Double_t kL){
    rewind(f1);
    char match_str[80];
    sprintf(match_str, "M_Zp=%i kL=%.2f", M_Zp, kL);
    char line[120];
    while(fgets(line, 120, f1)){
        //printf("%s \n", line);
        if(strstr(line,match_str) != NULL) {
            //printf("Line found! %s \n", line);
            sscanf(line, "M_Zp=%i kL=%lf: %lf %lf %lf  \n", &M_Zp, &kL, 
                    &AFBs[0], &AFBs[1], &AFBs[2]);
            return;
        }
    }
    printf("Couldn't find AFBs for M_Zp=%i kL=%.2f \n", M_Zp, kL);
    AFBs[0] = AFB_SM[0];
    AFBs[1] = AFB_SM[1];
    AFBs[2] = AFB_SM[2];
    return;
}

