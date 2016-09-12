#! /bin/bash

module load bcftools/1.3.1
module load tabix

cd /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome
: << 'comment'
# adding header info to each newly generated file:
ls IndelsByChr/* >list
readarray chrs < list
rm list

for file in "${chrs[@]}" ; do
	bgzip -d ${file}
done 
comment

ls IndelsByChr/1000G*vcf> list
readarray chrs < list
rm list
 
for file in "${chrs[@]}" ; do
	bcftools view -h 1000G_phase1.indels.hg19.sites.vcf.gz > header.txt
	echo ${file} 
	cat header.txt ${file} > tmp.vcf 
	mv tmp.vcf ${file}
done

: << 'END'	   
ls IndelsByChr/Mills* > list
readarray chrs < list
rm list

for file in "${chrs[@]}" ; do
	bcftools view -h Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz  > header.txt
	cat header.txt ${file} > ${file}.complete
	mv ${file}.complete ${file}
done

END

module unload bcftools/1.3.1
module unload tabix
