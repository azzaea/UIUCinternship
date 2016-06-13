#!/bin/bash

################################################################################################ 
# Program to calculate raw variants from human samples of WES short reads
# In order to run this pipeline please type at the command line
# start.sh <runfile>
################################################################################################

#PBS -S /bin/bash 
#PBS -q default
#PBS -d /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/src-azza
#PBS -j oe
#PBS -l nodes=1:ppn=23



echo -e "\n\n########################################################################################"
echo -e "#############                All looks good, proceed to alignement              ###############"
echo -e "########################################################################################\n\n"


# module load fastqc/0.11.4
# fastqc /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read1.fq

# fastqc /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read2.fq

# Download the reference for chromosome 1, from NCBI:
# wget ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/CHR_01/hs_ref_GRCh38.p2_chr1.fa.gz

# unzip:
# gunzip hs_ref_GRCh38.p2_chr1.fa.gz

# The reference is thus in: /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/genome
# and its name is: hs_ref_GRCh38.p2_chr1.fa

#################3 preparing the reference: I did this kind of late, so expect the date info to be messy!
module load samtools/1.2
samtools faidx /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/hs_ref_GRCh38.p2_chr1.fa


#module load bwa/0.7.10
#bwa index -a is /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/hs_ref_GRCh38.p2_chr1.fa

#bwa mem -M -t 16 /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/hs_ref_GRCh38.p2_chr1.fa /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read1.fq /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read2.fq > /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/aligned.sam

#module load samtools/1.2
#samtools flagstat /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/aligned.sam > alignmenet_summary.txt

#samtools view -bS /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/aligned.sam > /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/aligned.bam



module load picard-tools/1.141

# I actually ran this command on my node itself!!!
#java -jar picard.jar FastqToSam F1=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read1.fq F2=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read2.fq   O=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/fastq_to_bam.bam     SM=for_tool_testing 

java -jar picard.jar MergeBamAlignment ALIGNED=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/aligned.bam UNMAPPED=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/unmapped.bam \fastq_to_bam.bam    O=merge_alignments.bam      R=reference_sequence.fasta

#################################################################################################### 
################################ What are the group read info??? ################################
####################################################################################################


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
