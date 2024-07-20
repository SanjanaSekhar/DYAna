To modify templates after jobs finish making them:

```
root -l analyze/LQ_merge_workspaces.root
python combine/combine_eehighlowpTsys.py
# To unify cos binning
python combine/unify_cos_bins.py
python combine/LQ_make_mcstat_sys_templates.py
python combine/LQ_rename_sys_corrs.py
# To replace fakedata with real data
python combine/LQ_merge_data_templates.py
# To delete last two rap bins
python combine/blind_rap_bins.py
```
