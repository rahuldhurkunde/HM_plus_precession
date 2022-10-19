#!/bin/bash

dir='/work/rahul.dhurkunde/HM_and_precession'
psd='aplus'
f_min=15.0
tau0_threshold=0.8
nsplits=10
nworkflow=10

HMs=1
precession=1

first=0
last=50000

exec_dir="$dir/scripts/$psd/HM_prec"

echo "Executing in directory $exec_dir"
echo "Injs from $first to $last"
cp $dir/scripts/submit.sh $exec_dir
cp $dir/scripts/remove_dax.sh $exec_dir
cp $dir/scripts/combine_FFs.py $exec_dir
cd $exec_dir
mkdir results

$dir/scripts/FF_parser --config-files $dir/Config/$psd/workflow.ini \
					--template_bank $dir/banks/$psd/sorted_bank.hdf \
					--injection_dir $dir/injections/final_injections/ \
					--psd_file $dir/psds/$psd.txt \
					--HMs $HMs \
					--precession $precession \
					--f_min $f_min \
					--first $first \
					--last $last \
					--tau0_threshold $tau0_threshold \
					--nsplits $nsplits \
					--nworkflow $nworkflow 
