#!/bin/bash
#type = "data"
#echo $type
rm files.txt
rm list_to_be_deleted.txt
eval `scramv1 runtime -sh`
#eoscms
cat bad_eos_folders.txt | while read bad
do
		echo $bad
 cmsLs /store/user/bbilin/ntuples/data/$bad >> files.txt

done

echo "Looplar sarıldı hacı! :)"

cat files.txt | while read fil
do
echo "/"${fil#*/}
	echo "/"${fil#*/}>>list_to_be_deleted.txt

done

cat list_to_be_deleted.txt | while read filess

do
echo $filess
cmsRm $filess
done

cat bad_eos_folders.txt | while read bad
do
		echo $bad
 cmsRmdir /store/user/bbilin/ntuples/$type/$bad >> files.txt

done
