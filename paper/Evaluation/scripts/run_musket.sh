#!/bin/sh

cd simulate
for cov in 35 70 140
do
	for err in 1 3
	do
		../musket -p 10 -inorder -k 17 500000000 -o musket_simulate_${cov}_${err}.fq ~/data/Lighter/simulate/ecoli_reads_${cov}_${err}.fq
	done
done
cd ../

cd ERR022075
../musket -p 10 -inorder -k 17 500000000 -omulti musket ~/data/Lighter/ERR022075_75x_trim/read_1.fq ~/data/Lighter/ERR022075_75x_trim/read_2.fq
cd ..

cd gage_chr14
../musket -p 10 -inorder -k 19 11000000000 -omulti musket ~/data/Lighter/gage_chr14/frag_1.fastq ~/data/Lighter/gage_chr14/frag_2.fastq
cd ..

