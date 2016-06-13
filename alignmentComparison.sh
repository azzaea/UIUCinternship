
module load samtools
echo Total_reads percentage_mapped Mean_MAPQ > alignment_summary_changing_k.txt
for i in $(seq 10 2 30); do
	total=$(samtools view -c a$i.bam)
	mapped=$(samtools view -c -F0x4 a$i.bam)
	mapq=$(samtools view a$i.bam | awk '{sum+=$5} END { print "Mean MAPQ = ",sum/NR}')
	echo total mapped mapq >> alignment_summary_changing_k.txt
done

module load R/3.2.3



