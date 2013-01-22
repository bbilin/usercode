

cat config_data.txt | while read MAHMUT
	do
	echo $MAHMUT
	crab -create -cfg $MAHMUT
	crab -submit		
	done




#int totdim=__td__;
#double mpl=__mp__;
#double mmin__mm__;
