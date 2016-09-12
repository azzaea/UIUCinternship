#!/bin/bash

################################################################################################ 
	#################### Torque preparation: PBS commands ###################
################################################################################################

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads"
results="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results"

module load picard-tools/2.4.1
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar
module load samtools/1.3.1

cd $results/p3_postalignment

java -Xmx32G -jar $picard MarkDuplicates INPUT=cleaned_bam.bam OUTPUT=markduplicates.bam METRICS_FILE=markduplicates_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$results/tmp


java -Xmx32G -jar $picard MarkDuplicatesWithMateCigar INPUT=cleaned_bam.bam OUTPUT=markduplicatesMateCigar.bam METRICS_FILE=markduplicatesMateCigar_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$results/tmp

##### Doing comparisons/ Quality checking:

ls markduplica*bam > list
readarray bams < list
rm list

cd findings
echo
echo "samtools' summary" >  average_quality.txt
echo bamfile Mean_MAPQ >>  average_quality.txt

for bam in "${bams[@]}"; do
	bam=$(echo $bam | tr -d ' ')
	mapq=$(samtools view ../$bam | awk '{sum+=$5} END { print sum/NR}')
	echo $bam $mapq >> average_quality.txt

	java -jar $picard CollectWgsMetrics I=../$bam O=coverage_metrics_$bam.txt R=$reference/ucsc.hg19.fasta
	# produces metrics relating to reads that pass base and mapping quality filters and coverage (read depth) levels (user defined)
	# to use this tool with WES, you need to give the coordinates of your genomic region!!

done



echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
