#!/bin/bash

source /home/mhartz/.bashrc
cd /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter

./bin/cpPrecision -c1 input_hk.card -c2 input_hk.card -o /var/tmp/cp_prec_2det_hk_total_${1}.root -co osc_param.card -cs syst_param_total.card -d $1 &> /var/tmp/cp_prec_2det_hk_total_${1}.log

mv /var/tmp/cp_prec_2det_hk_total_${1}.root /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter/
mv /var/tmp/cp_prec_2det_hk_total_${1}.log /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter/

