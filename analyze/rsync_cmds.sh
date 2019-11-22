#!/bin/bash
rsync -r -v --delete -e ssh oamram@cmslpc-sl6.fnal.gov:/uscms_data/d3/oamram/DY_analysis/src/Analysis/DYAna/analyze/output_files/2016/*.root  output_files/2016/
rsync -r -v --delete -e ssh oamram@cmslpc-sl6.fnal.gov:/uscms_data/d3/oamram/DY_analysis/src/Analysis/DYAna/analyze/output_files/2017/*.root  output_files/2017/
rsync -r -v --delete -e ssh oamram@cmslpc-sl6.fnal.gov:/uscms_data/d3/oamram/DY_analysis/src/Analysis/DYAna/analyze/output_files/2018/*.root  output_files/2018/
ls
