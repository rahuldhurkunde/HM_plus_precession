for i in {0..19}; do
start=`echo "10.0* $i + 1.0" | bc`
end=`echo "10.0* ($i + 1) + 1.0" | bc`
echo $i $start $end

OMP_NUM_THREADS=1 \
condor_run -a accounting_group=cbc.prod.search -a request_memory=40000 unbuffer pycbc_brute_bank \
--verbose \
--output-file bank-$i.hdf \
--minimal-match 0.97 \
--tolerance .001 \
--buffer-length 4 \
--sample-rate 2048 \
--tau0-threshold 0.8 \
--approximant SEOBNRv4_ROM \
--tau0-crawl 4 \
--tau0-start $start \
--tau0-end $end \
--params mass1 mass2 spin1z spin2z \
--min 5   1  -0.99 -0.05 \
--max 30  3  0.99  0.05 \
--psd-file ZERO_DET_high_P.txt \
--seed 1 \
--low-frequency-cutoff 30.0  >  $i.out &
done

