cat config_mc.txt | while read MAHMUT
	do
	echo $MAHMUT
	crab -create -cfg $MAHMUT
	crab -submit		
	done
