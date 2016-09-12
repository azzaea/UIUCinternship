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

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar 

cd $results/p3_postalignment

java -jar $GenomeAnalysisTK \
              -T DiagnoseTargets \
              -R $reference/ucsc.hg19.fasta \
              -I markduplicates.bam \
              -L chr1 \
              -o output.vcf




echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
