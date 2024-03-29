#!/bin/bash

dir='/work/rahul.dhurkunde/HM_and_precession'
psd='ZDHP'
f_min=15.0
tau0_threshold=1.4
tb_splits=10

HMs=0
precession=0

first=3000
last=50000

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	echo "Injs from $first to $last"
	cp submit.sh $1
	cd $1
	rm gw.dax gw.ini tc.txt

	$dir/scripts/FF_parser --config-files $dir/Config/$psd/workflow.ini \
							--template_bank $dir/banks/$psd/sorted_bank.hdf \
							--injection_dir $dir/injections/final_injections/ \
							--output_dir output/ \
							--HMs $HMs \
							--precession $precession \
							--f_min $f_min \
							--start $first \
							--end $last \
							--tau0_threshold $tau0_threshold \
							--tb_splits $tb_splits 
	./submit.sh
fi
#							--injection_dir $dir/injections/small_injections/ \
