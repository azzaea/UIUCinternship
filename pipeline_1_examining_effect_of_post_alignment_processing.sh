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

# HC needs to be called for 3 different files. While at it, let's compare the run time in each case. For fun, let's do it also for the file that is just aligned, just realigned, in addition to the best scoring files from BQSR. These are: 
# 1. $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam
# 2. $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.5.bam
# 3. $results/p4_indelRealign_bqsr/bqsr/default/recal.ics.7.bam
# 4. $results/p4_indelRealign_bqsr/realignedBam.default.0.bam
# 5. $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam  (instead of $results/a.default.0.bam  as there are no duplicates in the dataset anyway, and it would save time for putting the ReadGroup info) 
# Note that the last couple of files will need their read group info added first!

if [ ! -d $results/p5_Haplotypecaller ]; then
	mkdir -p $results/p5_Haplotypecaller
fi

cd $results/p5_Haplotypecaller
if [ $not_interested ]; then
# 1. $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
	-T HaplotypeCaller\
	-R $reference/ucsc.hg19.fasta\
	-I $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam\
	-L chr1\
	--dbsnp $reference/dbsnp_138.hg19.vcf\
	-o vars_default_bqsr.g.vcf\
	-ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (BQSR default)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_bqsr.g.vcf\
        -o vars_default_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (BQSR default)" >> "hc.summary-various.stages.txt"

# 2. $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.5.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.5.bam\
        -L chr1\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_ics5_bqsr.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (BQSR ics5)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_ics5_bqsr.g.vcf\
        -o vars_ics5_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (BQSR ics5)" >> "hc.summary-various.stages.txt"
fi

# 3. $results/p4_indelRealign_bqsr/bqsr/default/recal.ics.7.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.7.bam\
        -L chr1\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_ics7_bqsr.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (BQSR ics7)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_ics7_bqsr.g.vcf\
        -o vars_ics7_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (BQSR ics7)" >> "hc.summary-various.stages.txt"

if [ $not_interested ]
# 4. $results/p4_indelRealign_bqsr/realignedBam.default.0.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/realignedBam.default.0.bam\
        -L chr1\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_default_realigned.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (Realigned default)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_realigned.g.vcf\
        -o vars_default_realigned.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (Realigned default)" >> "hc.summary-various.stages.txt"

# 5. $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam\
        -L chr1\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_default_markDuplicates.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (MarkDuplicates default)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_markDuplicates.g.vcf\
        -o vars_default_markDuplicates.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (MarkDuplicates default)" >> "hc.summary-various.stages.txt"
fi

echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
