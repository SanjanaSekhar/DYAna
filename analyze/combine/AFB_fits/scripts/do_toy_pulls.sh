#!/bin/bash
set -x 

mask_params="mask_Y16_emu16=1,mask_Y17_emu17=1,mask_Y18_emu18=1" 

for idx in 0 1 2 3 4 5
do
#idx=$1

    comb_card=cards/combined_fit_mbin${idx}.txt
    workspace=workspaces/combined_fit_toy_impacts_${idx}.root
    pars_corr="alphaS,alphaDen,RENORM,FAC,dy_xsec,bk_xsec,gam_xsec,"
    pars16="Pu16,BTAG16,elScaleStat16,elScaleSyst16,elScaleGain16,elSmear16,elHLT16,elID16,elRECO16,muRC16,muID16,muHLT16,lumi16,ee16_qcd,mu16_qcd,emu16_qcd,R_ee16_os_fakes,"
    pars17="Pu17,BTAG17,elScaleStat17,elScaleSyst17,elScaleGain17,elSmear17,elHLT17,elID17,elRECO17,muRC17,muID17,muHLT17,lumi17,ee17_qcd,mu17_qcd,emu17_qcd,R_ee17_os_fakes,"
    pars18="Pu18,BTAG18,elScaleStat18,elScaleSyst18,elScaleGain18,elSmear18,elHLT18,elID18,elRECO18,muRC18,muID18,muHLT18,lumi18,ee18_qcd,mu18_qcd,emu18_qcd,R_ee18_os_fakes"
    pars="${pars_corr}${pars16}${pars17}${pars18}"
    #pars="${pars_corr}"


    text2workspace.py $comb_card --keyword-value M_BIN=${idx} -P Analysis.DYAna.my_model:dy_AFB -o $workspace --channel-masks
    combine -M FitDiagnostics $workspace --saveWorkspace --skipBOnlyFit --robustFit 1 -t 1 --toysNoSystematics --setParameters ${mask_params},Afb=0,A0=0.05 #-v 10
    python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics.root -p Afb --skipFitB -g impacts/toy_pulls_mbin${idx}.root
done



