#!/bin/bash

startround=400
endround=499
nrounds=$(($endround+1-$startround))
echo "TOTAL rounds = $nrounds"
dir='/work/rahul.dhurkunde/HM_and_precession'

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	sed -n '8p' < $dir/Config/FF_config.ini
	sed -n '9p' < $dir/Config/FF_config.ini
	cp submit.sh $1/
	cd $1
	for i in $(seq $startround $endround); do
	echo $i
		
	mkdir $i
	cp submit.sh $i/
	cd $i/

	$dir/scripts/FF_parser --config-files $dir/Config/workflow.ini \
							--psd_file $dir/psds/ZDHP.txt \
							--template_bank $dir/banks/ZDHP/combined_bank.hdf \
							--injections $dir/injections/50000_inj/aligned_injections/$i.hdf \
							--tau_crawl $dir/injections/tau_files/tau_crawl_50000_aligned.txt \
							--output_dir output/	
	./submit.sh
	cd ../
	done
fi
