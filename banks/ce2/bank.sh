if [ $# -eq 0 ]
then
	echo "Provide the psd file"
	exit 1
else
		echo "Computing aligned bank using psd $1"
		for i in {275..275}; do
		start=`echo "1.315* $i + 1.0" | bc`
		end=`echo "1.315* ($i + 1) + 1.0" | bc`
		echo $i $start $end

		OMP_NUM_THREADS=1 \
		condor_run -a accounting_group=cbc.prod.search -a request_memory=60000 unbuffer pycbc_brute_bank \
		--verbose \
		--output-file bank-$i.hdf \
		--minimal-match 0.97 \
		--tolerance .001 \
		--buffer-length 4 \
		--sample-rate 2048 \
		--tau0-threshold 0.5 \
		--approximant SEOBNRv4_ROM \
		--tau0-crawl 4 \
		--tau0-start $start \
		--tau0-end $end \
		--params mass1 mass2 spin1z spin2z \
		--min 5   1  -0.99 -0.05 \
		--max 30  3  0.99  0.05 \
		--psd-file $1 \
		--seed 1 \
		--low-frequency-cutoff 7.0  >  $i.out &
		done
fi
