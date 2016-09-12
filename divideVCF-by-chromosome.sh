#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=3
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

referencedir="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome/genome"

for i in seq 1 22 ; do
  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar\
    -R ${referencedir}/ucsc.hg19.fasta\
    -T SelectVariants\
    -v ${referencedir}/1000G_phase1.indels.hg19.sites.vcf\
    -L chr${i}\
    -o ${referencedir}/IndelsByChr/1000G.chr${referencedir}.vcf
done

for i in seq 1 22 ; do
  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar\
    -R ${referencedir}/ucsc.hg19.fasta\
    -T SelectVariants\
    -v ${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
    -L chr${i}\
    -o ${referencedir}/IndelsByChr/1000G.chr${referencedir}.vcf
done
