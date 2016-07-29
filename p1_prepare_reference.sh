#!/bin/bash

#PBS -S /bin/bash 
#PBS -q default
#PBS -d /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/src-azza
#PBS -l nodes=1:ppn=23
#PBS -M aeahmed@illinois.edu
#PBS -M bea
#PBS -o localhost:$HOME/outputs-reference.log.txt
#PBS -e localhost:$HOME/errors-reference.log.txt 

echo -e "\n\n########################################################################################"
echo -e "##############                 Starting NOW                            ##################"      
echo -e "########################################################################################\n\n"

reference=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome
redmine=hpcbio-redmine@igb.illinois.edu
reportticket=5747
email=aeahmed@illinois.edu


# 1. Generating bwa index:
module load bwa/0.7.15
bwa index -a bwtsw $reference/HG19_GATKbundle2.8_noDecoys.fa
bwaexit=$?
echo BWA index generation exited with status=$bwaexit



# 2. Generating fasta index:
module load samtools/1.3.1
samtools faidx $reference/HG19_GATKbundle2.8_noDecoys.fa
faidxexit=$?
echo faidx index generation exited with status=$faidxexit


# 3. Generating sequence dictionary:
module load picard-tools/2.4.1
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar
java -jar $picard CreateSequenceDictionary\
	REFERENCE=$reference/HG19_GATKbundle2.8_noDecoys.fa\
	OUTPUT=$reference/HG19_GATKbundle2.8_noDecoys.dict
picardexit=$?
echo Dictionary generation exited with status=$picardexit

# 4. Reporting success to redmine and email
echo "Preparing new reference genome for synthetic dataset successful!" | mail -s "[Task #${reportticket}]" "$redmine,$email"

echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"      
echo -e "########################################################################################\n\n"
                                                                                            42,1          Bot


