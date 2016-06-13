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
echo -e "#############                BEGIN Quality checking              ###############"
echo -e "########################################################################################\n\n"

module load fastqc/0.11.4
fastqc /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read1.fq

fastqc /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads/H3A_NextGen_assessment.Chr1_50X.set3_read2.fq

echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
