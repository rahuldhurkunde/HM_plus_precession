#!/bin/bash

dir='/work/rahul.dhurkunde/HM_and_precession'

if [ $# -eq 0 ]
then
	echo "Provide sub-directory name"
	exit 1
else
	echo "Executing in directory $1"
	cp submit.sh $1/
	cd $1
	for i in {0..0}; do
	echo $i
		
	mkdir $i
	cp submit.sh $i/
	cd $i/

	$dir/scripts/FF_parser --config-files $dir/Config/workflow.ini --template_bank $dir/banks/parallel/aligned_bank/combined_bank.hdf --injections $dir/injections/50000_inj/aligned_injections/$i.hdf --tau_crawl $dir/injections/tau_files/flat_tau_crawl_0.5.txt --output_dir output/	
	./submit.sh
	cd ../
	done
fi
