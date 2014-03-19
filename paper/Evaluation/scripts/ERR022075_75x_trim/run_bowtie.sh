#!/bin/sh
rm summary.out
for prog in orig quake musket bless lighter
do
#echo "# ${prog}_${cov}_${err}" >> summary.out
#../../verify ${prog}/${prog}_simulate_${cov}_${err}.fq | head -13 | tail -4 >> summary.out
	cd $prog 
	bowtie2 -p 8 -x ~/data/lighter/ecoli_genome/e_coli_k12 -1 ${prog}_read1.fq -2 ${prog}_read2.fq > ${prog}.sam
	cd ..
done
