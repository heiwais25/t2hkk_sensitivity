#!/bin/tcsh -f 
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":/reno/home/heiwais25/T2HKKSensitivity/Fitter/lib"

set workdir = /reno/home/heiwais25/T2HKKSensitivity/Fitter 
cd $workdir

set Dect1 = hk
set Dect2 = mt_bisul 
# // hk // kd_1p5 // kd_2p0 // kd_2p5 // mt_bisul

# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output 
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output_both_dect
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/for_kps
#set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_min
set outdir = /reno/home/heiwais25/T2HKKSensitivity/combine_source/
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_max
set Type_mode = sin23Precision
 #  // cpPrecision  // cpViolation // massHierarchy // sin23Precision

#set output = sin23_prec_mt_bisul.root

#./bin/$Type_mode -c1 input_$Dect2.card  -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs mod_syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect2.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs syst_param_total2.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova.card -cs syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs syst_param_total.card -d $1

./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/${Type_mode}2hk_output_run_$1.root -co osc_param.card -cs syst_param_total.card -s $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova_m.card -cs syst_param_total.card -d $1

