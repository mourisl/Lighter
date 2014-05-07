#!/bin/sh
cp ../e_coli_k12.fa ref1.fa
perl mutate.pl ref1.fa > ref2.fa

let cnt=35*4700000/101
mason illumina -i -n 101 -N $cnt -pi 0 -pd 0 -pmm 0.01 -pmmb 0.005 -pmme 0.03 -sq -hs 0 -hi 0 -o ref1_reads.fq ref1.fa
mason illumina -i -n 101 -N $cnt -pi 0 -pd 0 -pmm 0.01 -pmmb 0.005 -pmme 0.03 -sq -hs 0 -hi 0 -s 17 -o ref2_reads.fq ref2.fa
cat ref1_reads.fq ref2_reads.fq > diploid.fq
