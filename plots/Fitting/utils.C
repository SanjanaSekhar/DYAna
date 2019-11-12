
TH1F* convert2d(TH2F *h_2d){
    int n_xf_bins = h_2d->GetNbinsX();
    int n_cost_bins = h_2d->GetNbinsY();

    char name[40];
    sprintf(name, "h1_%s", h_2d->GetName());
    TH1F *h_1d = new TH1F(name, "",  n_xf_bins * n_cost_bins, 0, n_xf_bins*n_cost_bins);
    for(int i=1; i<=n_xf_bins; i++){
        for(int j=1; j<= n_cost_bins; j++){
            float content = h_2d->GetBinContent(i,j);
            float error = h_2d->GetBinError(i,j);
            int gbin = (i-1)*n_cost_bins + j;
            //printf("gbin %i: i j %i %i \n", gbin, i, j);
            h_1d->SetBinContent(gbin, content);
            h_1d->SetBinError(gbin, error);
        }
    }
    return h_1d;
}
