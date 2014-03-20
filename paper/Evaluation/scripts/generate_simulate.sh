#!/bin/sh

for i in 25 #10 20 280
do
	let cnt=$i*4700000/101
	#echo $cnt 
	mason illumina -i -n 101 -N $cnt -pi 0 -pd 0 -pmm 0.01 -pmmb 0.005 -pmme 0.03 -sq -hs 0 -hi 0 -o simulate/ecoli_reads_${i}_1.fq e_coli_k12.fa
done
exit 

#cov: 35x, 70x, 140x
for i in 35 70 140
do
	let cnt=$i*4700000/101
	#echo $cnt 
	mason illumina -i -n 101 -N $cnt -pi 0 -pd 0 -pmm 0.01 -pmmb 0.005 -pmme 0.03 -sq -hs 0 -hi 0 -o simulate/ecoli_reads_${i}_1.fq e_coli_k12.fa
done
for i in 35 70 140
do
	let cnt=$i*4700000/101
	mason illumina -i -n 101 -N $cnt -pi 0 -pd 0 -pmm 0.03 -pmmb 0.015 -pmme 0.09 -sq -hs 0 -hi 0 -o simulate/ecoli_reads_${i}_3.fq e_coli_k12.fa
done
