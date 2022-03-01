#!/bin/bash

dir='/work/rahul.dhurkunde/HM_and_precession'
psd='ZDHP'
f_min=15.0
tau0_threshold=0.4

HMs=0
precession=0

first=0
last=300

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	cp submit.sh $1
	cd $1
	rm gw.dax gw.ini tc.txt

	$dir/scripts/FF_parser --config-files $dir/Config/$psd/workflow.ini \
							--template_bank $dir/banks/$psd/combined_bank.hdf \
							--injection_dir $dir/injections/100000_inj/ \
							--output_dir output/ \
							--HMs $HMs \
							--precession $precession \
							--f_min $f_min \
							--start $first \
							--end $last \
							--ninj_per_file 10 \
							--tau0_threshold $tau0_threshold
	./submit.sh
fi
#							--template_bank $dir/banks/parallel/small_bank/combined_bank.hdf \
#							--injection_dir $dir/injections/small_injections/ \
