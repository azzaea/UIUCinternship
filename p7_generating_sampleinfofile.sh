#!/bin/bash

dataset=synthetic_1_to_20_rare_to_common_individuals_ration
reads="/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/$dataset/generated_reads_$dataset"

samples=synthetic.sample*read1.fq.*

ls $reads/$samples > $reads/my_reads

if [ -s  $reads/$dataset.sampleinfo ]; then 
	truncate -s 0  $reads/$dataset.sampleinfo
fi

echo sample R1 R2

while read read_line; do
	R1=$(echo $read_line)
	R2=$(echo $R1| sed -e 's/read1.fq/read2.fq/')
	sample=${R1##*/} #${var##pattern}:delete longest match of pattern from beginning
	sample=${sample%_*} #${var%%pattern}      # delete longest match of pattern from the end
	echo $sample $R1 $R2 >> $reads/$dataset.sampleinfo
done < $reads/my_reads

rm $reads/my_reads

