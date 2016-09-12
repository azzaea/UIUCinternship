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

#### Needed modules:
module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar 

module load picard-tools/2.4.1
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar

module load samtools/1.3.1

module load bedtools/2.24.0 

module load R/3.2.3

gunzip $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf*
gunzip $reference/1000G_phase1.indels.hg19.sites.vcf*


mkdir $results/p4_indelRealign_bqsr
cd $results/p4_indelRealign_bqsr

cp $results/p3_postalignment/markduplicates.* .

# Default run
java -jar $GenomeAnalysisTK \
	-T  RealignerTargetCreator \
	-R $reference/ucsc.hg19.fasta \
	-I  markduplicates.bam\
	-L chr1 \
	-known $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-known $reference/1000G_phase1.indels.hg19.sites.vcf \
	-o realigner.default.0.intervals

cut -d'-' -f2 realigner.default.0.intervals|cut -d':' -f2 >ends
cut -d':' -f2 realigner.default.0.intervals|cut -d':' -f2|cut -d'-' -f1 > starts

Rscript -e 'ends=read.table("ends",stringsAsFactors=F);starts=read.table("starts",stringsAsFactors=F);interval.count=nrow(ends);interval.mean.length=sum(ends-starts)/interval.count; interval=data.frame(interval.count,interval.mean.length);write.table(x=interval,file="tmp",append=T)'

# Try to get the genomic coverage from the interval file:
sed  's/:\|-/\t/gi'  realigner.default.0.intervals > intervals.bed
genomeCoverageBed -d -i intervals.bed -g /home/apps/bedtools/bedtools-2.10.0/genomes/human.hg19.genome > coverage.hist.withoutbam.txt
# this command doesnt work cause some intervals contain only a start position (ie the third column entry for some rows is empty... you can view a certain line using:
# sed -n 4p file --> line 4
# sed -i '4d' file --> delete line 4
# sed -i '4i txt' file --> insert txt @ line 4
# you need to find a way to dect if the entry is empty to repeat the previous entry)

java -Xmx4g -Djava.io.tmpdir=$results/tmp -jar $GenomeAnalysisTK \
	-T IndelRealigner \
	-I  markduplicates.bam\
	-R $reference/ucsc.hg19.fasta \
	-known $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-known $reference/1000G_phase1.indels.hg19.sites.vcf \
	-targetIntervals realigner.default.0.intervals \
	-o realignedBam.default.0.bam 
	#remeber to set --consensusDeterminationModel KNOWNS_ONLY when you are dealing with large dataset
# -Djava.io.tmpdir sets the folder to put temproary files. Equivalent to Picards' tmp_dir

genomeCoverageBed -d -ibam realignedBam.default.0.bam  > coverage.hist.frombam.txt
# -d: compute coverage per base >> output is chromosome   position    depth(per base)


if [ $realigned_reads_only_bam_is_needed ]; then
	#To subset realigned reads only into a valid BAM:
	samtools view realignedBam.default.0.bam | grep 'OC' | cut -f1 | sort |uniq > realigned_reads_count.txt # create a list of readnames. 

	# Queryname sort the input BAM
	java -jar $picard SortSam INPUT=realignedBam.default.0.bam OUTPUT=xyz_querynamesort.bam SORT_ORDER=queryname

	# Create a new BAM containing read sets
	java -jar $picard FilterSamReads INPUT=xyz_querynamesort.bam OUTPUT=extraced_realigned.bam \
	FILTER=includeReadList READ_LIST_FILE=realigned_OC.txt SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=$results/tmp
fi



# idxstat's output is not exactly what I needed: coverage per interval
#samtools index realignedBam.default.0.bam  #needed for idxstat to work
#samtools idxstats realignedBam.default.0.bam > detailed_coverage.txt
#sed -i -e '1iRef.Seq.Name       SequenceLength       MappedReads    UnmappedReads\' detailed_coverage.txt 




#echo -e "\n\n########################################################################################"
#echo -e "#############                CHECKING PARAMETERS                         ###############"
#echo -e "########################################################################################\n\n"

if [ $not_needed] ; then
	# Actually, these are the parameters that you would need to change for the RealignerTargetCreator. But then, this is meaningless if you don't know how to measure success!
	declare -a parameters=(maxIntervalSize minReadsAtLocus windowSize)
	declare -a min=(200 2 5)
	declare -a step=(150 2 5)
	declare -a max=(800 10 40)

	cd $results/p2_alignment/
	mkdir ${parameters[@]}

	echo The parameters being tested and their ranges are given below:
	echo paramters: ${parameters[@]}, 
	echo minimum  : ${min[@]} 
	echo maximum  : ${max[@]}

	pos=0
	while [ $pos -lt ${#parameters[@]} ]; do
	        par=${parameters[pos]}
		cd $results/p2_alignment/$par
	        for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
			START=$(date +%s)	
			bwa mem -t 12 -$par "$i" -M -R  '@RG\tID:foo\tSM:bar\tLB:library1'  $reference/human $reads/read1.fq $reads/read2.fq  > "a.$par.$i.sam"
			END=$(date +%s)
			DIFF=$(( $END - $START ))
			if [ -s "a.$par.$i.sam" ]; then
				echo "Alignment successeful! with -$par $i" 
			else 
				echo "BWA aligned :0: using default parameter (*=0=)"> "a.$par.$i.summary.txt"
				echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
				continue
			fi
			alignments=$(samtools view -c a.$par.$i.sam)
			if [ "$alignments" -eq 0]; then
				echo			
				echo " Unfortunately, I can NOT process the parameter $par = $i with bwa mem" > "a.$par.$i.summary.txt"
				echo
				echo "BWA aligned :0: using default parameter (*=0=)" > "a.$par.$i.summary.txt"
				echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
				continue
			fi

			echo 
			echo "BWA Mem aligned :$alignments: using parameter :$par: =$i" > "a.$par.$i.summary.txt"
			echo 
			echo "Execution time is :$DIFF: seconds" >> "a.$par.$i.summary.txt"
			echo 
			samtools view -bS "a.$par.$i.sam" >"a.$par.$i.bam"
			samtools flagstat a.$par.$i.bam >> "a.$par.$i.summary.txt"
		done
	let pos+=1
	done
	
	
#igv mAYBE
fi


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
