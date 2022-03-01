#!/bin/bash
dir='/work/rahul.dhurkunde/HM_and_precession'

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	cp submit.sh $1
	cd $1
	rm gw.dax gw.ini tc.txt

	$dir/scripts/FF_parser --config-files $dir/Config/workflow.ini \
							--injection_dir $dir/injections/100000_inj/nonaligned_injections/ \
							--output_dir output/ \
							--start 0 \
							--end 300
	./submit.sh
fi
#							--injection_dir $dir/injections/small_injections/ \
