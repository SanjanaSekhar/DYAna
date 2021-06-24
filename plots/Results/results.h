
const int n_m_bins = 7;
Double_t m[n_m_bins] = {185, 225., 285., 415., 605., 850., 1200.};
Double_t m_err[n_m_bins] = {15., 25., 35., 95., 95., 150., 200.};

Double_t y_powheg[n_m_bins] =      {0.611, 0.596, 0.597, 0.592, 0.590, 0.590, 0.590};
Double_t y_powheg_errs[n_m_bins] = {0.004, 0.005, 0.004, 0.004, 0.004, 0.004, 0.004};

Double_t y_amc[n_m_bins] =         {0.612, 0.608, 0.601, 0.605, 0.605, 0.607, 0.610};
Double_t y_amc_errs[n_m_bins] =    {0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002};

Double_t A0_amc[n_m_bins] =        {0.065, 0.057, 0.063, 0.050, 0.038, 0.026, 0.025};
Double_t A0_amc_errs[n_m_bins] =   {0.007, 0.004, 0.005, 0.006, 0.008, 0.005, 0.006};


Double_t *y_comb =      y_amc;
Double_t y_comb_errs[n_m_bins] = {0.013, 0.014, 0.014, 0.013, 0.023, 0.033, 0.052};

Double_t *y_mumu =      y_amc;
Double_t y_mumu_errs[n_m_bins] = {0.016, 0.016, 0.0175, 0.017, 0.028, 0.042, 0.065};

Double_t * y_elel=      y_amc;
Double_t y_elel_errs[n_m_bins] = {0.018, 0.018, 0.020, 0.019, 0.035, 0.049, 0.077};



Double_t A0_comb[n_m_bins] =        {0.056, 0.047, 0.055, 0.042, 0.030, 0.018, 0.012};
Double_t A0_comb_errs[n_m_bins] =   {0.015, 0.018, 0.022, 0.025, 0.055, 0.080, 0.130};

Double_t A0_mumu[n_m_bins] =        {0.056, 0.047, 0.055, 0.042, 0.030, 0.018, 0.012};
Double_t A0_mumu_errs[n_m_bins] =   {0.019, 0.019, 0.026, 0.029, 0.065, 0.100, 0.180};

Double_t A0_elel[n_m_bins] =        {0.056, 0.047, 0.055, 0.042, 0.030, 0.018, 0.012};
Double_t A0_elel_errs[n_m_bins] =   {0.019, 0.019, 0.026, 0.029, 0.065, 0.100, 0.180};
