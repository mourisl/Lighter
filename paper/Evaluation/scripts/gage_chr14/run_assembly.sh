#!/bin/sh

for prog in orig quake musket bless lighter
do
	cp SD2.config ${prog}/${prog}.config
	echo "q1=${prog}_read1.fq" >> ${prog}/${prog}.config
	echo "q2=${prog}_read2.fq" >> ${prog}/${prog}.config 
done

#exit

for prog in orig quake musket bless lighter
do
	cd ${prog} 
	for i in 67 73 77 #43 47 53 57 63
	do
		SOAPdenovo2-127mer all -s ${prog}.config -K $i -p 8 -o SD2_${prog}_$i > SD2_${prog}_$i.log 
	done
	cd .. 
done
#nohup SOAPdenovo2 all -s EC.config -K 75 -p 8 -o EC_75/EC > EC_75.log &
