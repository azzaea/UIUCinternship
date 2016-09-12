#! /bin/bash


module load bcftools/1.3.1
module load tabix

cd /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome

# It wasn't necessary to gzip the files in the last step, so I'm unzipping here!!!
ls IndelsByChr/* >list
readarray chrs < list
rm list

for file in "${chrs[@]}" ; do
	bgzip -d ${file}
done 

# Now, add the header to the 1000G vcf files:
ls IndelsByChr/1000G*vcf> list
readarray chrs < list
rm list
 
bcftools view -h 1000G_phase1.indels.hg19.sites.vcf.gz > header.txt
for file in "${chrs[@]}" ; do
	cat header.txt ${file} > tmp.vcf 
	mv tmp.vcf ${file}
done

# Now, add the header to the Mills vcf files:
ls IndelsByChr/Mills.* > list
readarray chrs < list
rm list

bcftools view -h Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz  > header.txt
for file in "${chrs[@]}" ; do
	cat header.txt ${file} > tmp.vcf
	mv tmp.vcf ${file}
done

rm header.txt

module unload bcftools/1.3.1
module unload tabix
