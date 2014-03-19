rm snp_count
for prog in quake musket bless lighter
do
	echo $prog >> snp_count 
	cd $prog
	samtools view -bS ${prog}.sam > $prog.bam
	samtools sort $prog.bam $prog.sorted

	samtools mpileup -uD -f ~/data/lighter/e_coli_k12.fa $prog.sorted.bam | bcftools view -bvcg - > $prog.bcf
	bcftools view $prog.bcf | vcfutils.pl varFilter -D 500 | wc -l >> ../snp_count 
	cd ..
done
