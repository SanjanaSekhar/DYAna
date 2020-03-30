!/bin/bash
set -x 

hadd -f output_files/2016/MuMu16_comb_back_nov20.root output_files/2016/MuMu16_ttbar_nov20.root output_files/2016/MuMu16_diboson_nov20.root  output_files/2016/MuMu16_wt_nov20.root
hadd -f output_files/2016/MuMu16_fakes_contam_nov20.root output_files/2016/MuMu16_comb_back_nov20.root output_files/2016/MuMu16_dy_nov20.root output_files/2016/MuMu16_photInd_nov20.root
./prune_root_file.sh output_files/2016/MuMu16_fakes_contam_nov20.root

hadd -f output_files/2017/MuMu17_comb_back_nov20.root output_files/2017/MuMu17_ttbar_nov20.root output_files/2017/MuMu17_diboson_nov20.root  output_files/2017/MuMu17_wt_nov20.root
hadd -f output_files/2017/MuMu17_fakes_contam_nov20.root output_files/2017/MuMu17_comb_back_nov20.root output_files/2017/MuMu17_dy_nov20.root output_files/2017/MuMu17_photInd_nov20.root
./prune_root_file.sh output_files/2017/MuMu17_fakes_contam_nov20.root

hadd -f output_files/2018/MuMu18_comb_back_nov20.root output_files/2018/MuMu18_ttbar_nov20.root output_files/2018/MuMu18_diboson_nov20.root  output_files/2018/MuMu18_wt_nov20.root
hadd -f output_files/2018/MuMu18_fakes_contam_nov20.root output_files/2018/MuMu18_comb_back_nov20.root output_files/2018/MuMu18_dy_nov20.root output_files/2018/MuMu18_photInd_nov20.root
./prune_root_file.sh output_files/2018/MuMu18_fakes_contam_nov20.root


hadd -f output_files/2016/ElEl16_comb_back_nov20.root output_files/2016/ElEl16_ttbar_nov20.root output_files/2016/ElEl16_diboson_nov20.root  output_files/2016/ElEl16_wt_nov20.root
hadd -f output_files/2016/ElEl16_fakes_contam_nov20.root output_files/2016/ElEl16_comb_back_nov20.root output_files/2016/ElEl16_dy_nov20.root output_files/2016/ElEl16_photInd_nov20.root
./prune_root_file.sh output_files/2016/ElEl16_fakes_contam_nov20.root

hadd -f output_files/2017/ElEl17_comb_back_nov20.root output_files/2017/ElEl17_ttbar_nov20.root output_files/2017/ElEl17_diboson_nov20.root  output_files/2017/ElEl17_wt_nov20.root
hadd -f output_files/2017/ElEl17_fakes_contam_nov20.root output_files/2017/ElEl17_comb_back_nov20.root output_files/2017/ElEl17_dy_nov20.root output_files/2017/ElEl17_photInd_nov20.root
./prune_root_file.sh output_files/2017/ElEl17_fakes_contam_nov20.root

hadd -f output_files/2018/ElEl18_comb_back_nov20.root output_files/2018/ElEl18_ttbar_nov20.root output_files/2018/ElEl18_diboson_nov20.root  output_files/2018/ElEl18_wt_nov20.root
hadd -f output_files/2018/ElEl18_fakes_contam_nov20.root output_files/2018/ElEl18_comb_back_nov20.root output_files/2018/ElEl18_dy_nov20.root output_files/2018/ElEl18_photInd_nov20.root
./prune_root_file.sh output_files/2018/ElEl18_fakes_contam_nov20.root
