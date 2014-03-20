#!/bin/sh
rm summary.out
for prog in quake musket bless lighter
do
	for cov in 35 70 140
	do
		for err in 1 3
		do
			echo "# ${prog}_${cov}_${err}" 
			echo "# ${prog}_${cov}_${err}" >> summary.out
			../../verify ${prog}/${prog}_simulate_${cov}_${err}.fq | head -13 | tail -4 >> summary.out
		done
	done
done
