#!/bin/sh
rm summd
for prog in orig quake musket bless lighter
do
	echo $prog >> summd
	perl ../../Tools/SumMD.pl ${prog}/${prog}.sam >> summd 
done
