#!/bin/bash

################################################################################################ 
	#################### Torque preparation: PBS commands ###################
################################################################################################

### Use the bourne shell
#PBS -S /bin/bash 			

### Run in the queue named default
#PBS -q default				

#### #PBS -d /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/src-azza

### Specify the number of cpus for your job (nodes and processors)
#PBS -l nodes=1:ppn=23			

### the destination for reporting (defaults to script's directory): Seperate these to get a better sense!
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt		

### email me when my job aborts or ends
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

mkdir -p $results/p3_postalignment/findings
cd $results/p3_postalignment

##### Doing comparisons/ Quality checking:

ls *.?am > list
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

	java -jar $picard CollectAlignmentSummaryMetrics R=$reference/ucsc.hg19.fasta I=../$bam O=overall_summary_$bam.txt
	# produces metrics relating to the overall alignment of reads withing a SAM/BAM. ie it works per bam file
	# The reference file needs a companion dictionary here! it means, use the GATK reference, not the bwa indexed file

	java -jar $picard CollectWgsMetrics I=../$bam O=coverage_metrics_$bam.txt R=$reference/ucsc.hg19.fasta
	# produces metrics relating to reads that pass base and mapping quality filters and coverage (read depth) levels (user defined)
	# to use this tool with WES, you need to give the coordinates of your genomic region!!

	java -jar $picard CollectInsertSizeMetrics I=../$bam O=insert_size_metrics_$bam.txt H=insert_size_histogram_$bam.pdf M=0.5
	# produces metrics for validating library construction including the insert size distribution and read 	orientation of paired end libraries
	# setting the minimum percentage (M=0.5) is convineint for processing a small file

done



echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
