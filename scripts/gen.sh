#!/bin/bash

dir='/work/rahul.dhurkunde/HM_and_precession'
psd='aplus'
f_min=15.0
tau0_threshold=0.8
nsplits=10
nworkflow=10

HMs=0
precession=0

first=0
last=50000

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	echo "Injs from $first to $last"
	cp submit.sh $1
	cp remove_dax.sh $1
	cp combine_FFs.py $1
	cd $1
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
	#./submit.sh
fi
#							--injection_dir $dir/injections/small_injections/ \
#							--injection_dir $dir/injections/final_injections/ \
