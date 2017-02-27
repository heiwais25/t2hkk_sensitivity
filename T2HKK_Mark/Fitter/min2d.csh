#!/bin/tcsh -f 
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":/reno/home/heiwais25/T2HKK_Mark/Fitter/lib"

set workdir = /reno/home/heiwais25/T2HKK_Mark/Fitter 
cd $workdir
#@ one = 80/40
# 40 bin #
#@ first = $1 / 40
#@ second = $1 % 40
# 100 bin #
@ first = $1 / 100
@ second = $1 % 100
@ runtime = 1645
#echo "$first + $second"
# first = $1
#calc "2*pi"
set Dect1 = hk
#echo "$one $other"
set Dect2 = hk
# // hk // kd_1p5 // kd_2p0 // kd_2p5 // mt_bisul

# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output 
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output_both_dect
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/for_kps
#set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_min
set outdir = /reno/home/heiwais25/T2HKK_Mark/combine_source/
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_max
set Type_mode = cpViolationPOT
 #  // cpPrecision  // cpViolation // massHierarchy // sin23Precision // sin23Violation // min2D
 #  // min2D_new	// cpViolation_time
#set output = sin23_prec_mt_bisul.root

#./bin/$Type_mode -c1 input_$Dect2.card  -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs mod_syst_param_total.card -d $1

# ./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/source_$1_run_1845.root -co osc_param.card -cs syst_param_total.card -y $1

# Two HK detector case
./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/source_$1_run_$runtime.root -co osc_param.card -cs syst_param_total_2HK.card -y $1

# Two HK + other detector case
# ./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/source_$1_run_$runtime.root -co osc_param.card -cs syst_param_total.card -y $1

# KD one detector
# ./bin/$Type_mode -c1 input_$Dect1.card -o $outdir/source_$1_run_$runtime.root -co osc_param.card -cs syst_param_total_KDonly.card -y $1

# JD one detector
# ./bin/$Type_mode -c1 input_$Dect1.card -o $outdir/source_$1_run_$runtime.root -co osc_param.card -cs syst_param_total_2HK.card -y $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova.card -cs syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/${Type_mode}_output_run_$1.root -co osc_param.card -cs syst_param_total.card -s $1

# change 2 oscillation parameter
# ./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect1/${Type_mode}_sin23_0.50_${first}_dm23_0.0024_${second}_run_gr_10000_mt_bisul_reactor_new.root -co osc_param.card -cs syst_param_total.card -p1 0.0024 -p2 0.50 -n1 dm32 -n2 ssqth23 -s1 ${first} -s2 ${second}

# change 1 oscillation parameter
#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/${Type_mode}_sin23_0.45_${first}_run_300_2hk_reactor.root -co osc_param.card -cs syst_param_total.card -p1 0.45 -n1 ssqth23 -n2 dm32 -s1 ${first} 


#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova_m.card -cs syst_param_total.card -d $1

