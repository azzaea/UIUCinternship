#!/bin/bash

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

module load R/3.2.3


if [ ! -d $results/p5_Haplotypecaller/param_sweep ]; then
	mkdir -p $results/p5_Haplotypecaller/param_sweep
fi

cd $results/p5_Haplotypecaller/param_sweep

echo Run using the default parameters: > hcparams.txt
echo  "parameters=(stand_emit_conf stand_call_conf mbq kmerSize minDanglingBranchLength )" >> hcparams.txt
echo  "default=   (       30                30     10   10,25            4 ) "  >> hcparams.txt

echo "-------------------------------------------------------------------------------------------------" >> hcparams.txt

# 1. $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
	-T HaplotypeCaller\
	-R $reference/ucsc.hg19.fasta\
	-I $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam\
	-L chr1\
	--dbsnp $reference/dbsnp_138.hg19.vcf\
	-o vars_default.0.g.vcf\
	-ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (BQSR default)" >> "hcparams.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_0.g.vcf\
        -o vars_default_0.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (BQSR default)" >> "hcparams.txt"

declare -a parameters=(stand_emit_conf stand_call_conf mbq kmerSize minDanglingBranchLength )
declare -a default=(30 30 10 10,25 4 ) 
declare -a min=(0 0 4 10 1 )
declare -a step=(10 10 4 15 2)
declare -a max=(50 50 20 26 10)

mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below: >> hcparams.txt

echo paramters: ${parameters[@]}, >> hcparams.txt

echo minimum  : ${min[@]} >> hcparams.txt

echo maximum  : ${max[@]} >> hcparams.txt

echo >> hcparams.txt
 
echo  "------------------------------------------------------------------------------------------------------------" >> hcparams.txt


pos=0
echo start while loop ;

while [ $pos -lt ${#parameters[@]} ]
do
	par=${parameters[pos]}
	cd $results/p5_Haplotypecaller/param_sweep/$par
	for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T HaplotypeCaller\
		        -R $reference/ucsc.hg19.fasta\
		        -I $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam\
		        -L chr1\
		        --dbsnp $reference/dbsnp_138.hg19.vcf\
		        -o vars_$par.$i.g.vcf\
		        -ERC GVCF\
			-$par "$i"
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
		echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (BQSR default)" >> "hcparams.txt"

		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T GenotypeGVCFs\
		        -R $reference/ucsc.hg19.fasta\
		        -V vars_$par.$i.g.vcf\
		        -o vars_$par_$i.vcf
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
		echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (BQSR default)" >> "hcparams.txt"
	done
let pos+=1
done


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
