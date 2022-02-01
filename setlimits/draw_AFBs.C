
#include "../plots/tdrstyle.C"
#include "../plots/CMS_lumi.C"
#include "pythia_lhe_reader.C"

Double_t AFB_utype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return SM_AFB(E, 0);
}

Double_t AFB_dtype(Double_t *x, Double_t *params){
    Double_t E = x[0];
    return SM_AFB(E, 1);
}

Double_t Zp_AFB_utype(Double_t *x, Double_t *params){
    Double_t M_Zp = 3000.; 
    Double_t E = x[0];
    return Zp_AFB(E, M_Zp, 0);
}

Double_t Zp_E6_AFB_utype(Double_t *x, Double_t *params){
    Double_t M_Zp = 3000.; 
    Double_t E = x[0];
    return Zp_AFB(E, M_Zp, 0, 1);
}

void draw_AFBs(){
    setTDRStyle();
    Double_t m_min = 10.;
    Double_t m_max = 2000.;
    TF1 *AFB_d = new TF1("AFB (d-type)", AFB_dtype, m_min,m_max);
    TF1 *AFB_u = new TF1("AFB (u-type)", AFB_utype, m_min,m_max);

    TF1 *Zp_AFB_u = new TF1("Zp_AFB_u", Zp_AFB_utype, m_min,m_max);

    TF1 *Zp_E6_AFB_u = new TF1("Zp_AFB_E6_u", Zp_E6_AFB_utype, m_min,m_max);
    AFB_u->SetLineColor(kBlue);
    Zp_AFB_u->SetLineColor(kGreen);
    Zp_E6_AFB_u->SetLineColor(kRed);

    AFB_u->SetLineWidth(4);
    Zp_AFB_u->SetLineWidth(4);
    Zp_E6_AFB_u->SetLineWidth(4);
    TCanvas *can = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    can->SetLogx();
    AFB_u->Draw();
    //AFB_d->Draw("same");
    Zp_AFB_u->Draw("same");
    //Zp_E6_AFB_u->Draw("same");


    AFB_u->GetYaxis()->SetRangeUser(-0.75, 0.75);
    AFB_u->GetYaxis()->SetTitle("A_{FB}");
    AFB_u->GetYaxis()->SetNdivisions(503);
    AFB_u->GetYaxis()->SetTitleSize(0.075);
    AFB_u->GetYaxis()->SetTitleOffset(1.0);
    AFB_u->GetYaxis()->SetLabelSize(0.05);
    AFB_u->GetYaxis()->CenterTitle();
    // X axis AFB_u plot settings
    AFB_u->GetXaxis()->SetTitle("Dilepton Mass (GeV)");
    AFB_u->GetXaxis()->SetTitleSize(0.05);
    AFB_u->GetXaxis()->SetTitleOffset(1.);
    AFB_u->GetXaxis()->SetLabelSize(0.05);
    AFB_u->GetXaxis()->SetTickLength(0.04);


    //can->SetLogx();

    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.6, 0.5, 0.85, 0.7);
    leg1->AddEntry(AFB_u, "SM (u-type)", "l");
    //leg1->AddEntry(AFB_d, "SM (d-type)", "l");
    leg1->AddEntry(Zp_AFB_u, "SM + Z'_{SSM} (3.0 Tev)", "l");
    //leg1->AddEntry(Zp_E6_AFB_u, "SM + Z'_{#chi} (1.5 Tev)", "l");
    leg1->SetTextSize(0.04);
    leg1->Draw();
}
