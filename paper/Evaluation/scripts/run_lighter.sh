#!/bin/sh

for cov in 35 70 140
do
	if [ $cov -eq 35 ]
	then
		alpha=0.1
	elif [ $cov -eq 70 ]
	then
		alpha=0.05
	else
		alpha=0.025
	fi

	for err in 1 3
	do
		./lighter -t 10 -r ~/data/Lighter/simulate/ecoli_reads_${cov}_${err}.fq -k 17 5000000 ${alpha} -od simulate
	done
done

./lighter -t 10 -r ~/data/Lighter/ERR022075_75x_trim/read_1.fq -r ~/data/Lighter/ERR022075_75x_trim/read_2.fq -k 17 5000000 0.05 -od ERR022075
cd ERR022075
mv read_1.cor.fq lighter_read1.fq
mv read_2.cor.fq lighter_read2.fq
cd ..

./lighter -t 10 -r ~/data/Lighter/gage_chr14/frag_1.fq -r ~/data/Lighter/gage_chr14/frag_2.fq -k 19 110000000 0.1 -od gage_chr14
cd gage_chr14
mv frag_1.cor.fq lighter_read1.fq
mv frag_2.cor.fq lighter_read2.fq
cd ..

