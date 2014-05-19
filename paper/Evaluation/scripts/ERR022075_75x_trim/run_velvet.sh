#!/bin/sh


#exit

export OMP_NUM_THREADS=10
export OMP_THREAD_LIMIT=16

for prog in orig quake musket bless lighter
do
	cd ${prog} 
	for i in 43 47 49 51 53 55 57 63 67
	do
		velveth velvet_${i} ${i} -fmtAuto -separate -shortPaired ${prog}_read1.fq ${prog}_read2.fq
		velvetg velvet_${i} -ins_length 311 -exp_cov auto -scaffolding no
		perl ../../../Tools/VelvetN50.pl velvet_${i}/contigs.fa > velvet_${i}.N50
	done
	cd .. 
done
