for i in {0..50}; do
echo $i

mkdir $i
cp submit.sh $i/
cd $i/
../FF_parser --config-files ../../Config/workflow.ini --template_bank ../../banks/parallel/aligned_bank/combined_bank.hdf --injections ../../injections/$i.hdf --output_dir output/
./submit.sh
cd ../
done
