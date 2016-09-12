#! /bin/bash

module load bcftools/1.3.1
module load tabix

cd /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome

# Splitting up the reference by chromosome/ contig:

cut -f1 ucsc.hg19.fasta.fai | xargs -i echo tabix 1000G_phase1.indels.hg19.sites.vcf.gz {} \| bgzip \> {}.1000G.vcf.gz \&|sh
cut -f1 ucsc.hg19.fasta.fai | xargs -i echo tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz {} \| bgzip \> {}.1000G.vcf.gz \&|sh

# Adding the header to each newly generated file

# to this for each chromosome, you need to create a proper for loop. This is suffficient though to test a single chromosome:

bcftools view -h 1000G_phase1.indels.hg19.sites.vcf.gz > test.txt
cat test.txt IndelsByChr/chr1.1000G.vcf >bigchr1.1000G.vcf
cp bigchr1.1000G.vcf IndelsByChr/


bcftools view -h Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz > test.txt
cat test.txt IndelsByChr/chr1.Mills.vcf >tmp.txt
mv tmp.txt IndelsByChr/chr1.Mills.vcf

module unload bcftools/1.3.1
module unload tabix
