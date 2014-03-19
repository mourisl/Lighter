#!/bin/sh
rm lighter.memusage
for cov in 35 70 140
do
	for err in 1 #3
	do
		if [ $cov = 35 ]
		then
			alpha=0.1
		elif [ $cov = 70 ]
		then 
			alpha=0.05
		else
			alpha=0.025
		fi
		
		#echo $alpha
		echo $cov >> lighter.memusage 
		~/Tools/memusage ./a.out -t 10 -r ~/data/lighter/simulate/ecoli_reads_${cov}_${err}.fq -k 17 5000000 $alpha -od simulate 2>> lighter.memusage
	done
done
