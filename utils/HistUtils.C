
void set_fakerate_errors(TH2D *h_errs, TH2D *h_fr, TH1F *h){
    float err_sum = 0.;
    for(int i=1; i<= h_errs->GetNbinsX(); i++){
        for(int j=1; j<= h_errs->GetNbinsY(); j++){
            float err = h_fr->GetBinError(i,j);
            float num = h_errs->GetBinContent(i,j);
            err_sum += err*err*num*num;
        }
    }
    //increase error sum by 50% to approximate subtraction of 2 fail template
    err_sum *= 1.50;
    float n_err_events = h_errs->Integral();
    for(int i=1; i<= h->GetNbinsX(); i++){
        float bin_num = h->GetBinContent(i);
        float scaling = pow(bin_num/n_err_events, 2);
        float num_err = pow(h->GetBinError(i),2);
        float weight_err = scaling * err_sum;
        float new_err = sqrt(num_err + weight_err);
        //printf("Err was %.1f, is now %.1f \n", sqrt(num_err), new_err);

        h->SetBinError(i, new_err);
    }
}

void set_fakerate_errors(TH2D *h_errs, TH2D *h_fr, TH2F *h){
    float err_sum = 0.;
    for(int i=1; i<= h_errs->GetNbinsX(); i++){
        for(int j=1; j<= h_errs->GetNbinsY(); j++){
            float err = h_fr->GetBinError(i,j);
            float num = h_errs->GetBinContent(i,j);
            err_sum += err*err*num*num;
        }
    }
    //increase error sum by 50% to approximate subtraction of 2 fail template
    err_sum *= 1.50;
    float n_err_events = h_errs->Integral();
    for(int i=1; i<= h->GetNbinsX(); i++){
        for(int j=1; j<= h->GetNbinsY(); j++){
            float bin_num = h->GetBinContent(i, j);
            float scaling = pow(bin_num/n_err_events, 2);
            float num_err = pow(h->GetBinError(i,j),2);
            float weight_err = scaling * err_sum;
            float new_err = sqrt(num_err + weight_err);
            //printf("Err was %.1f, is now %.1f \n", sqrt(num_err), new_err);

            h->SetBinError(i,j, new_err);
        }
    }
}
    



double get_cost(TLorentzVector lep_p, TLorentzVector lep_m, bool do_reco_flip = true){

    TLorentzVector cm = lep_p + lep_m;
    double root2 = sqrt(2);
    double lep_p_pls = (lep_p.E()+lep_p.Pz())/root2;
    double lep_p_min = (lep_p.E()-lep_p.Pz())/root2;
    double lep_m_pls = (lep_m.E()+lep_m.Pz())/root2;
    double lep_m_min = (lep_m.E()-lep_m.Pz())/root2;
    double qt2 = cm.Px()*cm.Px()+cm.Py()*cm.Py();
    double cm_m2 = cm.M2();
    //cost_p = cos(theta)_r (reconstructed collins soper angle, sign
    //may be 'wrong' if lepton pair direction is not the same as inital
    //quark direction)
    double cost = 2*(lep_m_pls*lep_p_min - lep_m_min*lep_p_pls)/sqrt(cm_m2*(cm_m2 + qt2));
    if(cm.Pz() <0. && do_reco_flip) cost = -cost;
    return cost;
}

bool notCosmic(TLorentzVector v1, TLorentzVector v2){
    Double_t ang = v1.Angle(v2.Vect());
    return ang < (TMath::Pi() - 0.005);
}
Double_t compute_xF(TLorentzVector cm){
    return  abs(2.*cm.Pz()/13000.); 
}

bool goodElEta(Float_t el1_eta){
    el1_eta = abs(el1_eta);
    return (el1_eta < 1.4442) || (el1_eta > 1.566 && el1_eta < 2.5); 
}

double get_cost_v2(TLorentzVector lep_p, TLorentzVector lep_m){

    double Ebeam = 6500.;
    double Pbeam = sqrt(Ebeam*Ebeam - 0.938*0.938);
    TLorentzVector cm = lep_p + lep_m;
    TLorentzVector p1(0., 0., Pbeam, Ebeam);
    TLorentzVector p2(0., 0., -Pbeam, Ebeam);

    if(cm.Pz() < 0. ){
        TLorentzVector p = p1;
        p1 = p2;
        p2 = p;
    }

    TVector3 beta = -cm.BoostVector();
    lep_m.Boost(beta);
    lep_p.Boost(beta);
    p1.Boost(beta);
    p2.Boost(beta);

    // Now calculate the direction of the new z azis

    TVector3 p1u = p1.Vect();
    p1u.SetMag(1.0);
    TVector3 p2u = p2.Vect();
    p2u.SetMag(1.0);
    TVector3 bisec = p1u - p2u;
    bisec.SetMag(1.0);
    TVector3 plep = lep_m.Vect();
    double cost = TMath::Cos(plep.Angle(bisec));
    return cost;
}

bool has_no_bjets(Int_t nJets, Double_t jet1_pt, Double_t jet2_pt, 
        Double_t jet1_cmva, Double_t jet2_cmva){
    Double_t med_btag = 0.4432;
    if(nJets ==0) return true;
    else if(nJets == 1){
        if(jet1_pt < 20.) return true;
        else return jet1_cmva < med_btag;
    }
    else{
        return (jet1_pt < 20. || jet1_cmva < med_btag) && (jet2_pt < 20. || jet2_cmva < med_btag);
    }
}

typedef struct {
    TH2D *h;
} FakeRate;
//
//static type means functions scope is only this file, to avoid conflicts
void setup_new_el_fakerate(FakeRate *FR, int year){
    TFile *f0;
    if (year == 2016)      f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron16_data_fakerate_nov1.root");
    else if (year == 2017) f0 = TFile::Open("../analyze/FakeRate/root_files/2017/SingleElectron17_data_fakerate_nov1.root");
    else if (year == 2018) f0 = TFile::Open("../analyze/FakeRate/root_files/2018/SingleElectron18_data_fakerate_nov1.root");
    else printf("Year was %i ?? \n", year);
    TH2D *h1 = (TH2D *) gDirectory->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    FR->h = h1;
    f0->Close();
}
void setup_new_mu_fakerate(FakeRate *FR, int year){
    TFile *f0;
    if (year == 2016)      f0 = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon16_data_fakerate_oct30.root");
    else if (year == 2017) f0 = TFile::Open("../analyze/FakeRate/root_files/2017/SingleMuon17_data_fakerate_oct18.root");
    else if (year == 2018) f0 = TFile::Open("../analyze/FakeRate/root_files/2018/SingleMuon18_data_fakerate_oct18.root");
    else printf("Year was %i ?? \n", year);
    TDirectory *subdir = gDirectory;
    TH2D *h1 = (TH2D *) subdir->Get("h_rate_new")->Clone();
    h1->SetDirectory(0);
    FR->h = h1;
    f0->Close();
}


Double_t get_new_fakerate_prob(Double_t pt, Double_t eta, TH2D *h){
    //pt=35;
    if (pt >= 150) pt =150 ;
    if(abs(eta) > 2.4) eta = 2.3;

    Double_t variation = 0.;
    Double_t err;


    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    Double_t prob = h->GetBinContent(xbin, ybin);
    err = h->GetBinError(xbin, ybin);
    //printf("prob err %.2f %.2f \n", prob, err);
    prob += err * variation;
    //printf("prob: %.2f \n", prob);
    if(prob < 0.001 || prob >= 0.99){
        printf("Warning: %.2f Rate for pt %.0f, eta %1.1f! \n"
                "err %.2f varr %.2f \n",
                prob, pt, eta, err, variation);
        ybin -=1;
        prob = h->GetBinContent(xbin, ybin);
        err = h->GetBinError(xbin, ybin);
        prob += err * variation;
        if(prob < 0.001) printf("Tried 1 lower pt bin and still 0 fakerate \n");
    }

    prob = min(prob, 0.98);
    prob = max(prob, 0.02);
    //printf("Efficiency is %f \n", eff);
    return prob;
}
Double_t get_new_fakerate_unc(Double_t pt, Double_t eta, TH2D *h){
    //pt=35;
    if (pt >= 150) pt =150 ;
    if(abs(eta) > 2.4) eta = 2.3;

    Double_t err;


    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin = x_ax->FindBin(std::abs(eta));
    int ybin = y_ax->FindBin(pt);

    err = h->GetBinError(xbin, ybin);
    return err;
}
