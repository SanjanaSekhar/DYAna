#!/bin/bash


set -x

ls -la
pwd
cd Analysis/DYAna/analyze/combine/AFB_fits
./do_impacts.sh $3
xrdcp -f impacts/mbin${3}.json $1

