#!/bin/sh

cd simulate
for cov in 35 70 140
do
	for err in 1 3
	do
		../bless -prefix bless_simulate_${cov}_${err} -read ~/data/Lighter/simulate/ecoli_reads_${cov}_${err}.fq -kmerlength 17 
	done
done
cd ../

cd ERR022075
../bless -prefix bless_read -read1 ~/data/Lighter/ERR022075_75x_trim/read_1.fq -read2 ~/data/Lighter/ERR022075_75x_trim/read_2.fq -kmerlength 17 
cd ..

cd gage_chr14
../bless -prefix bless_read -read1  ~/data/Lighter/gage_chr14/frag_1.fastq -read2 ~/data/Lighter/gage_chr14/frag_2.fastq -kmerlength 19 
cd ..

