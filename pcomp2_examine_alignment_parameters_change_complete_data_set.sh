#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=2:ppn=23			
#PBS -o localhost:$HOME/outputs-bwacomp.log.txt
#PBS -e localhost:$HOME/errors-bwacomp.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m ae

# program input:
dataset=error_rates

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome"
reads="/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/$dataset/generated_reads_$dataset"
results="/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/results"

read1=synthetic.sample*_var_${dataset}_merged_read1.fq

ls $reads/$read1 > $reads/my_reads

if [ ! -d $results/$dataset ]; then
	mkdir -p $results/$dataset
fi

while read read_line ; do
	R1=$(echo $read_line)
	R2=$(echo $R1| sed -e 's/read1.fq/read2.fq/')
	sample=${R1##*/} #${var##pattern}:delete longest match of pattern from beginning
	sample=${sample%%_*} #${var%%pattern}      # delete longest match of pattern from the end
	if [ ! -d $results/$dataset/$sample ]; then
        	mkdir -p $results/$dataset/$sample
	fi

	######### Quality checking: 
	module load fastqc/0.11.5
	fastqc $reads/$read1 --outdir=$results/p1_quality
	fastqc $reads/$read2 --outdir=$results/p1_quality
	# In reality, to see this result, just type: firefox required.read_fastqc.html
	# Reports show data of good quality, except for some kmer content in read1 (Gloria: its near the end, so 	should be ok). Also, the encoding is Sanger, so the reads are neat..


######### Alignment: The default settings
if [ ! -d $results/p2_alignment/ ]; then
	mkdir -p $results/p2_alignment/
fi

cd $results/p2_alignment/

mkdir default
cd default

module load bwa/0.7.15


	
	START=$(date +%s)
	bwa mem -M -t 12 -R  '@RG\tID:H7.5.R1\tPL:illumina\tPU:H7LH3CCXX.5.0\tLB:R1\tPI:0\tDT:2016-7-1\tSM:NA12878-Garvan' $reference/human $reads/$read1 $reads/$read2 > a.default.0.sam
# note that it is the indexed file(s) that was needed in this stage. It has a different name than the 	reference (remember the -p argument), so I need to use its name (human)
	END=$(date +%s)
	[ -s a.default.0.sam ] && echo "Default alignment successeful!" || exit
	alignments=$(samtools view -c a.default.0.sam)
	if [ "$alignments" -eq 0]; then
		echo			
		echo " Unfortunately, I can NOT process the your request with default parameters" 
		echo
		exit
	fi
	samtools view -bS a.default.0.sam > a.default.0.bam

	DIFF=$(( $END - $START ))

	echo 
	echo "BWA Mem aligned :$alignments: using parameter :default: (*=0=)"  > a.default.0.summary.txt
	echo 
	echo "Execution time is :$DIFF: seconds" >> a.default.0.summary.txt
	echo 
	samtools flagstat a.default.0.bam >> a.default.0.summary.txt  # Generating summary statistics


done < $reads/my_reads
######### Alignment: The combinatorial settings: changing a variable at a time, with the objective of having a sense of how things work


#echo -e "\n\n########################################################################################"
#echo -e "#############                CHECKING PARAMETERS                         ###############"
#echo -e "########################################################################################\n\n"

declare -a parameters=(k r w d c D m W A B O E L U T)
declare -a min=(3 .5 20 20 300 .1 20 0 1 1 1 1 1 1 10)
declare -a step=(3 .5 20 20 300 .1 20 3 2 2 2 2 2 3 10)
declare -a max=(60 4 200 200 10000 1 200 30 20 20 20 20 20 40 80)

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
		bwa mem -t 12 -$par "$i" -M -R  '@RG\tID:H7.5.R1\tPL:illumina\tPU:H7LH3CCXX.5.0\tLB:R1\tPI:0\tDT:2016-7-1\tSM:NA12878-Garvan'  $reference/human $reads/read1.fastq $reads/read2.fastq  > "a.$par.$i.sam"
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

done #ending the loop that goes through each pair of files


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
