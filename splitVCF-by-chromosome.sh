#! /bin/bash


module load bcftools/1.3.1
module load tabix

cd /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome

# Generating a vcf.gz file if not already present:
bgzip -c 1000G_phase1.indels.hg19.sites.vcf > 1000G_phase1.indels.hg19.sites.vcf.gz
tabix -p vcf 1000G_phase1.indels.hg19.sites.vcf.gz

bgzip -c Mills_and_1000G_gold_standard.indels.hg19.sites.vcf > Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz 
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

# Splitting up the reference by chromosome/ contig:

cut -f1 ucsc.hg19.fasta.fai | xargs -i echo tabix 1000G_phase1.indels.hg19.sites.vcf.gz {} \| bgzip \> IndelsByChr/1000G.{}.vcf.gz \&|sh
cut -f1 ucsc.hg19.fasta.fai | xargs -i echo tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz {} \| bgzip \> IndelsByChr/Mills.{}.vcf.gz \&|sh

# Adding the header to each newly generated file
# first unzip your files!
ls IndelsByChr/*.gz > list
readarray chrs <list
rm list
for file in "${chrs[@]}" ; do
	bgzip -d IndelsByChr/${chrs[$file]}
done

ls IndelsByChr/1000G*vcf> list
readarray chrs < list
rm list
 
for file in "${chrs[@]}" ; do
	bcftools view -h 1000G_phase1.indels.hg19.sites.vcf.gz > header.txt
	cat header.txt IndelsByChr/${chrs[${file}]} > IndelsByChr/complete.${chrs[${file}]}
	mv IndelsByChr/complete.${chrs[${file}]} IndelsByChr/${chrs[${file}]}
done
	   
ls IndelsByChr/Mills* > list
readarray chrs < list
rm list

for file in "${chrs[@]}" ; do
	bcftools view -h Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz  > header.txt
	cat header.txt IndelsByChr/${chrs[${file}]} > IndelsByChr/complete.${chrs[${file}]}
	mv IndelsByChr/complete.${chrs[${file}]} IndelsByChr/${chrs[${file}]}
done

module unload bcftools/1.3.1
module unload tabix
