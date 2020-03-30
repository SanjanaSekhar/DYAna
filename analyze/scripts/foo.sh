#!/bin/bash

set -x



hadd -f output_files/2018/MuMu18_fakes_contam_feb17.root output_files/2018/MuMu18_dy_feb10.root output_files/2018/MuMu18_phot_ind_feb17.root output_files/2018/MuMu18_comb_back_feb18.root
./prune_root_file.sh output_files/2018/MuMu18_fakes_contam_feb17.root

hadd -f output_files/2018/ElEl18_fakes_contam_feb17.root output_files/2018/ElEl18_dy_feb10.root output_files/2018/ElEl18_phot_ind_feb17.root output_files/2018/ElEl18_comb_back_feb18.root
./prune_root_file.sh output_files/2018/ElEl18_fakes_contam_feb17.root


#hadd output_files/2016/MuMu16_comb_back_feb18.root output_files/2016/MuMu16_ttbar_feb17.root output_files/2016/MuMu16_diboson_feb17.root output_files/2016/MuMu16_WT_feb17.root
#hadd output_files/2016/ElEl16_comb_back_feb18.root output_files/2016/ElEl16_ttbar_feb17.root output_files/2016/ElEl16_diboson_feb17.root output_files/2016/ElEl16_WT_feb17.root
#
#
#hadd output_files/2017/MuMu17_comb_back_feb18.root output_files/2017/MuMu17_ttbar_feb17.root output_files/2017/MuMu17_diboson_feb17.root output_files/2017/MuMu17_WT_feb17.root
#hadd output_files/2017/ElEl17_comb_back_feb18.root output_files/2017/ElEl17_ttbar_feb17.root output_files/2017/ElEl17_diboson_feb17.root output_files/2017/ElEl17_WT_feb17.root
#
#hadd output_files/2018/MuMu18_comb_back_feb18.root output_files/2018/MuMu18_ttbar_feb17.root output_files/2018/MuMu18_diboson_feb17.root output_files/2018/MuMu18_WT_feb17.root
#hadd output_files/2018/ElEl18_comb_back_feb18.root output_files/2018/ElEl18_ttbar_feb17.root output_files/2018/ElEl18_diboson_feb17.root output_files/2018/ElEl18_WT_feb17.root
