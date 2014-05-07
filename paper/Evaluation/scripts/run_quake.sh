#!/bin/sh

for cov in 35 70 140
do
	for err in 1 3
	do
		quake.py -u -p 10 --hash_size=1000000 -k 17 -r ~/data/Lighter/simulate/ecoli_reads_${cov}_${err}.fq
		mv ~/data/Lighter/simulate/ecoli_reads_${cov}_${err}.cor.fq simulate/quake_simulate_${cov}_${err}.fq
	done
done
exit

cd ERR022075
quake.py -u -p 10 --hash_size=1000000 -k 17 -f reads_files
mv ~/data/Lighter/ERR022075_75x_trim/read_1.cor.fq quake_read1.cor.fq
mv ~/data/Lighter/ERR022075_75x_trim/read_2.cor.fq quake_read2.cor.fq
cd ..

cd gage_chr14
quake.py -u -p 10 --hash_size=1000000 -k 19 -f reads_files
mv ~/data/Lighter/gage_chr14/read_1.cor.fastq quake_read1.fq
mv ~/data/Lighter/gage_chr14/read_2.cor.fastq quake_read2.fq
cd ..
