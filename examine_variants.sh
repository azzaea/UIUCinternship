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

module load picard-tools/2.4.1
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar

cd $results/p4_indelRealign_bqsr
module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

############## The real code:
#Apply calibration to SNPs

#Lets apply the recalibration to our SNP calls. Lets try with a threshold of 99.0

java -Xmx3g -jar GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R chr21.fasta \
   -input HG00418.snps.raw.vcf.gz \
   --ts_filter_level 99.0 \
   -tranchesFile var_recal/HG00418.tranches \
   -recalFile var_recal/HG00418.recal \
   -o HG00418.snps.recal.vcf

#Open the recalibrated file and check the "Filter" column, you should now see that all SNPs that were classified at 99.0 or better are all termed "PASS" and the rest are termed "LowQual", "TruthSensitivityTranche99.00to99.90" or "TruthSensitivityTranche99.90to100.00". We can now get all the trusted SNPs by taking the ones with "PASS".

grep "#" HG00418.snps.recal.vcf > header.vcf
grep "PASS" HG00418.snps.recal.vcf | cat header.vcf - > HG00418.snps.recal.pass.vcf
grep -c "PASS" HG00418.snps.recal.vcf

#Q3. How many SNPs do we have that PASS the filtering?

#However keep in mind that the even at the 99.0% sensitivity tranche we expect quite a number of false positive calls (the orange and orange scribbled bars in the plot). Lets try to divide the calls into known and novel (based on dbSNP) and plot the quality distributions using R.

perl -ane 'if ($_ =~ m/^#/) {} else {if ($F[2] eq ".") {print $F[5], "\n"}}' HG00418.snps.recal.pass.vcf > novel.qual.txt
perl -ane 'if ($_ =~ m/^#/) {} else {if ($F[2] ne ".") {print $F[5], "\n"}}' HG00418.snps.recal.pass.vcf > known.qual.txt

R
n = read.table("novel.qual.txt")
k = read.table("known.qual.txt")
plot(density(k[,1]), xlim=c(0,500), ylim=c(0,0.01), xlab="Variant Quality", main="Known vs Novel SNPs")
points(density(n[,1]), type="l", col="red")
legend("topright", legend=c("known", "novel"), col=c("black", "red"), lty=1)
dev.print("variants.known_vs_novel.pdf", device=pdf)

#Q4. Based on this, would you trust all of the novel SNPs? A: not the ones with low quality

#We could remove all novel SNPs with lower quality than eg 70 (this would be to apply hard-filtering).

perl -ane 'if ($_ =~ m/^#/) {print $_} else {if ($F[2] eq "." and $F[5] < 70) {} else {print $_}}' \
HG00418.snps.recal.pass.vcf > HG00418.snps.recal.pass.novel_filt.vcf
grep -v "#" -c HG00418.snps.recal.pass.novel_filt.vcf

