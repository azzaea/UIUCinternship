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

module load picard-tools/2.4.1
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar

if [ $AddingReadGroupInfoAndIndexing ]; then
	cd $results/p4_indelRealign_bqsr

	java -jar $picard AddOrReplaceReadGroups\
	      I=markduplicates.bam\
	      O=markduplicates_RG_added.bam\
	      RGID=foo\
	      RGLB=library1\
	      RGPL=illumina\
	      RGPU=unit1\
	      RGSM=bar

	# You need to re-index your bam, for the GATK. It would have been better to add RG and index using novosort as it does both!
        module load novocraft
        novosort markduplicates_RG_added.bam -i -c 12 -o markduplicates_RG_added2.bam

	# Some GATK tools don't like gzipped files. unzip for safety!
	gunzip $reference/dbsnp_138.hg19.vcf*
	gunzip $reference/1000G_phase1.indels.hg19.sites.vcf*
	gunzip $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf*	
fi


# Default run
cd $results/p4_indelRealign_bqsr
if [ $defaultrun ];
then 
	module load gatk/3.6
	GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

	module load R/3.2.3

	java -jar $GenomeAnalysisTK\
		-T BaseRecalibrator\
		-R $reference/ucsc.hg19.fasta\
		-I markduplicates_RG_added2.bam\
		-L chr1\
		-knownSites $reference/dbsnp_138.hg19.vcf\
	  	-knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
	 	-knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
		-o recal.table.default.0

	exitcode=$?
	echo exit code from BaseRecalibrator = $exitcode

	java -jar $GenomeAnalysisTK\
        	-T BaseRecalibrator\
	        -R $reference/ucsc.hg19.fasta\
	        -I markduplicates_RG_added2.bam\
	        -L chr1\
	        -knownSites $reference/dbsnp_138.hg19.vcf\
	        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
	        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
	        -BQSR recal.table.default.0\
	        -o after_recal.table.default.0

	exitcode=$?
	echo exit code from the rerun of BaseRecalibrator = $exitcode

	java -jar $GenomeAnalysisTK\
		-T AnalyzeCovariates\
		-R $reference/ucsc.hg19.fasta\
		-L chr1\
		-before recal.table.default.0\
		-after after_recal.table.default.0\
		-csv my-report.csv\
		-plots recalibration_plots.default.0.pdf

	exitcode=$?
	echo exit code from AnalyzeCovariates = $exitcode

	java -jar $GenomeAnalysisTK\ 
		-T PrintReads\ 
		-R $reference/ucsc.hg19.fasta\ 
		-I markduplicates_2_RG.bam\ 
		-L chr1\
		-BQSR recal.table.default.0\
		-o recal.default.0.bam
	exitcode=$?
	echo exit code from PrintReads = $exitcode
fi

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

module load R/3.2.3


#if [ $varyparameters ] ; then
	declare -a parameters=(ics maxCycle mcs bqsrBAQGOP ddq idq lqt mdq)
	declare -a defaults=(3 500 2 40 45 45 2 -1)
	declare -a min=(1 250 1 10 10 10 1 2)
	declare -a step=(2 150 2 10 10 10 2 2)
	declare -a max=(13 1000 13 70 70 70 12 12) 
	declare -a exceptional_run="-ddq=-1  -idq=-1"

	cd $results/p4_indelRealign_bqsr/bqsr
	mkdir ${parameters[@]}

	echo The parameters being tested and their ranges are given below:
	echo paramters: ${parameters[@]}, 
	echo minimum  : ${min[@]} 
	echo maximum  : ${max[@]}

	pos=0
	while [ $pos -lt ${#parameters[@]} ]; do
	        par=${parameters[pos]}
		cd $results/p4_indelRealign_bqsr/bqsr/$par
	        for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
			START=$(date +%s)
			java -jar $GenomeAnalysisTK\
			        -T BaseRecalibrator\
			        -R $reference/ucsc.hg19.fasta\
			        -I $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam\
			        -L chr1\
			        -knownSites $reference/dbsnp_138.hg19.vcf\
			        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
			        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
				-$par "$i"\
			        -o recal.table.$par.$i

			END=$(date +%s)
                        DIFF=$(( $END - $START ))
                        if [ -s "recal.table.$par.$i" ]; then
				echo >> ../"bqsr.summary.txt"
                        else
                                echo "BQSR failed using -$par $i">> ../"bqsr.summary.txt"
                                echo "Execution time is :DIFF: seconds" >> ../"bqsr.summary.txt"
				echo
                                continue
                        fi
			
			exitcode=$?
		        echo exit code from BaseRecalibrator using -$par $i = $exitcode

			java -jar $GenomeAnalysisTK\
			        -T BaseRecalibrator\
			        -R $reference/ucsc.hg19.fasta\
			        -I $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam\
			        -L chr1\
			        -knownSites $reference/dbsnp_138.hg19.vcf\
			        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
			        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
			        -BQSR recal.table.$par.$i\
			        -o after_recal.table.$par.$i
			exitcode=$?
		        echo exit code from the rerun of  BaseRecalibrator with -$par $i = $exitcode

			java -jar $GenomeAnalysisTK\
			        -T AnalyzeCovariates\
			        -R $reference/ucsc.hg19.fasta\
			        -L chr1\
			        -before recal.table.$par.$i\
			        -after after_recal.table.$par.$i\
			        -plots recalibration_plots.$par.$i.pdf

			exitcode=$?
		        echo exit code from AnalyzeCovariates using $par $i = $exitcode

			START=$(date +%s)
			java -jar $GenomeAnalysisTK\
			        -T PrintReads\
			        -R $reference/ucsc.hg19.fasta\
			        -I $results/p4_indelRealign_bqsr/markduplicates_2_RG.bam\
			        -L chr1\
			        -BQSR recal.table.$par.$i\
			        -o recal.$par.$i.bam
			END=$(date +%s)
			DIFF=$(( $END -$START + $DIFF ))
			echo "Base Recalibration successeful! with -$par $i" >> ../"bqsr.summary.txt"
			echo "Execution time is : $DIFF : seconds" >> ../"bqsr.summary.txt"
			echo

		exitcode=$?
	        echo exit code from PrintReads with -$par $i = $exitcode

		done
	let pos+=1
	done
#fi


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
