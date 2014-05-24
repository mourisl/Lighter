First:
1. Change to the `paper/Evaluation` subdirectory of the repo
2. Obtain [Quake].  We used v0.3.
3. Obtain [Musket].  We used v1.1.
4. Obtain [BLESS].  We used v0p12.

[Quake]: http://www.cbcb.umd.edu/software/quake/
[Musket]: http://musket.sourceforge.net/homepage.htm#latest
[BLESS]: http://sourceforge.net/projects/bless-ec/

To generate simulated E. coli data set:
1. Download FASTA from refseq: http://www.ncbi.nlm.nih.gov/nuccore/NC_010473.1
2. Move `sequnce.fa` to current directory and rename 
3. Run `sh scripts/generate_simulate.sh`

Generating the diploid test case:
generate_diploid.sh 

Generating the 2x100bp 75x ERR022075 data set
perl Sample.pl 0.077 < ERR022075_1.fastq > read1.fq
perl Sample.pl 0.077 < ERR022075_1.fastq > read2_tmp.fq
perl Trim.pl read2_tmp.fq > read2.fq

Generating the kmers to test:
perl GetRefKmers.pl 17 < e_coli_k12.fa
perl GetKmersAroundpos.pl 17 ref1.fa ref2.fa mutate.pos

Softwares:
quake v0.3 : run_quake.sh
musket v1.1 : run_musket.sh
bless v0.12 : run_bless.sh
lighter : run_lighter.sh

Evaluating simulated data set:
Compile verify.cpp
simulate/run.sh

Evaluating real data set:
in the corresponding fold:
run_bowtie.sh
evaluate_bowite.sh
run_velvet.sh
(pick the parameter giving best N50)
run_quast.sh

Evaluating memory usage:
compile memusage.cpp (obtained from musket's project)
memusage (command to the program).
For musket, we tested different number of size of bloom filter and choose the result with minimum memory consumption.
e.g.
For gage_chr14 data set:
~/Tools/memusage ./musket -p 10 -inorder -k 19 2200000000 -omulti tmp ~/data/Lighter/gage_chr14/frag_1.fastq ~/data/Lighter/gage_chr14/frag_2.fastq
