#define STAND_ALONE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TLine.h"
#include "TTree.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "LHAPDF/LHAPDF.h"


using namespace LHAPDF;
using namespace std;

void add_pdf_weight(int year){



    float inc1_pdfweight, inc2_pdfweight, pdfweight;
    TLorentzVector *inc1(0), *inc2(0), *gen_lep_p(0), *gen_lep_m(0), cm;
    int inc_id1, inc_id2;

    float E_BEAM = 7500.;


    string pdfname;
    if(year == 2016)
        pdfname = string("NNPDF30_nlo_nf_5_pdfas");
    else
        pdfname = string("NNPDF31_nnlo_hessian_pdfas");
    const PDF* my_pdf = mkPDF(pdfname.c_str(),0);

    char mu_fname[100], el_fname[100];

    sprintf(mu_fname, "../../analyze/output_files/%i/MuMu%i_dy_aug05.root", year, year % 2000);
    sprintf(el_fname, "../../analyze/output_files/%i/ElEl%i_dy_aug05.root", year, year % 2000);


    TFile *f_elel_mc = new TFile(el_fname, "update");
    TTree *t_elel_mc = ((TTree *)f_elel_mc->Get("T_sig"));



    t_elel_mc->SetBranchAddress("inc1", &inc1);
    t_elel_mc->SetBranchAddress("inc2", &inc2);
    t_elel_mc->SetBranchAddress("gen_el_p", &gen_lep_p);
    t_elel_mc->SetBranchAddress("gen_el_m", &gen_lep_m);
    t_elel_mc->SetBranchAddress("inc_id1", &inc_id1);
    t_elel_mc->SetBranchAddress("inc_id2", &inc_id2);
    TBranch *inc1_pdfweight_el =  t_elel_mc->Branch("inc1_pdfweight", &inc1_pdfweight, "inc1_pdfweight/F");
    TBranch *inc2_pdfweight_el =  t_elel_mc->Branch("inc2_pdfweight", &inc2_pdfweight, "inc2_pdfweight/F");
    TBranch *pdfweight_el =  t_elel_mc->Branch("evt_pdfweight", &pdfweight, "evt_pdfweight/F");


    for (int i=0; i<t_elel_mc->GetEntries(); i++) {
        t_elel_mc->GetEntry(i);
        cm = *gen_lep_p + *gen_lep_m;

        float scale = cm.M() * cm.M() + cm.Pt() * cm.Pt();


        float x1 = std::abs(inc1->E() / E_BEAM);
        float x2 = std::abs(inc2->E() / E_BEAM);

        inc1_pdfweight = my_pdf->xfxQ2(inc_id1, x1, scale)/x1;
        inc2_pdfweight = my_pdf->xfxQ2(inc_id2, x2, scale)/x2;

        pdfweight = inc1_pdfweight * inc2_pdfweight;

        inc1_pdfweight_el->Fill();
        inc2_pdfweight_el->Fill();
        pdfweight_el->Fill();

    }

    f_elel_mc->cd();
    t_elel_mc->Write();
    f_elel_mc->Close();
    




    TFile *f_mumu_mc = new TFile(mu_fname, "update");
    TTree *t_mumu_mc = (TTree *)f_mumu_mc->Get("T_sig");

    t_mumu_mc->SetBranchAddress("inc1", &inc1);
    t_mumu_mc->SetBranchAddress("inc2", &inc2);
    t_mumu_mc->SetBranchAddress("gen_mu_p", &gen_lep_p);
    t_mumu_mc->SetBranchAddress("gen_mu_m", &gen_lep_m);
    t_mumu_mc->SetBranchAddress("inc_id1", &inc_id1);
    t_mumu_mc->SetBranchAddress("inc_id2", &inc_id2);
    TBranch *inc1_pdfweight_mu =  t_mumu_mc->Branch("inc1_pdfweight", &inc1_pdfweight);
    TBranch *inc2_pdfweight_mu =  t_mumu_mc->Branch("inc2_pdfweight", &inc2_pdfweight);
    TBranch *pdfweight_mu =  t_mumu_mc->Branch("evt_pdfweight", &pdfweight);

    for (int i=0; i<t_mumu_mc->GetEntries(); i++) {
        t_mumu_mc->GetEntry(i);
        cm = *gen_lep_p + *gen_lep_m;

        float scale = cm.M() * cm.M() + cm.Pt() * cm.Pt();


        float x1 = std::abs(inc1->E() / E_BEAM);
        float x2 = std::abs(inc2->E() / E_BEAM);

        inc1_pdfweight = my_pdf->xfxQ2(inc_id1, x1, scale)/x1;
        inc2_pdfweight = my_pdf->xfxQ2(inc_id2, x2, scale)/x2;

        pdfweight = inc1_pdfweight * inc2_pdfweight;

        inc1_pdfweight_mu->Fill();
        inc2_pdfweight_mu->Fill();
        pdfweight_mu->Fill();
    }

    f_mumu_mc->cd();
    t_mumu_mc->Write();
    f_mumu_mc->Close();

}





int main(int argc, char * argv[]){
    if (argc != 2){
        printf("Takes 1 argument (year) \n");
        exit(1);
    }
    int year = atoi(argv[1]);
    add_pdf_weight(year);
}
