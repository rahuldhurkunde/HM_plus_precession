if [ $# -eq 0 ]
then
	echo "Provide nworkflow"
	exit 1
else
	let end=$1-1
	echo "nworkflow $1 and last ind of loop $end"
	
	for i in $(seq 0 $end); do
		submit_dir="part_$i/"
		echo "Removing dax in sub dir $submit_dir"
		cd $submit_dir

		rm gw_$i.dax gw_$i.ini
		
		cd ../
	done
fi

