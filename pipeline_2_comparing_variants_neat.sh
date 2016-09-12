#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-neat_var_comp.log.txt
#PBS -e localhost:$HOME/errors-neat_var_comp.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads"
results="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results"
TargetedRegions="/home/groups/hpcbio_shared/azza/TargetedRegions"

#### Needed modules:

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

module load R/3.2.3

# HC was called for 5 different files, and it would be nice to compare the variants called in each case. These files are:
# 1. $results/p5_Haplotypecaller 
# 2. $results/p5_Haplotypecaller
# 3. $results/p5_Haplotypecaller
# 4. $results/p5_Haplotypecaller
# 5. $results/p5_Haplotypecaller
# The objectives of this exercise:
## 1) See the effect of different preprocessig steps in the called variants
## 2) See how the different programs compare variants, and understand which is approapriate in this study. What is special here is that were are particularily interested in rare variants, and so the comparison needs to be very exact. In and ideal world, the comparison would further tell me where the 2 files agree and disagree (besides the stats on sensitivity and specifity)

# ######################################################################################################################## #
# The following tools should be tested:
# 1. vcf_compare from UIUC
# 2. varEval from GATK
# 3. hap.py from Illumina (https://github.com/Illumina/hap.py/tree/ee33640b218b5c575e999be619884992781ab9ac)
# 4. rtg-tools from RealTimeGenomics (https://github.com/RealTimeGenomics/rtg-tools/tree/fdeeb3554badf1286fffc413d57913a9b7dfaa46) 
#               note that the latest updates on hap.py and rtg-tools date back to 2015 (Sept, and Oct respectively)
# ######################################################################################################################## #


cd $results/p5_Haplotypecaller/comparison_stages/

##
#1. Using vcf_compare.py from UIUC

module load python/2.7.9
vcf_compare_dir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/builds/NEAT/neat-genreads/utilities

if [ $done ]; then

START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf -w $results/p5_Haplotypecaller/vars_default_bqsr.vcf -o default_bqsr_neat $results/tmp --incl-homs --incl-fail  -t $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR default)" >> "../variants_comparison.summary.txt"

###
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf -w $results/p5_Haplotypecaller/vars_ics5_bqsr.vcf -o ics5_bqsr_neat $results/tmp --incl-homs --incl-fail  -t $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed -T 90 --vcf-out  --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR ics5)" >> "../variants_comparison.summary.txt"
##
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf -w $results/p5_Haplotypecaller/vars_ics7_bqsr.vcf -o ics7_bqsr_neat $results/tmp --incl-homs --incl-fail  -t $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed -T 90 --vcf-out  --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR ics7)" >> "../variants_comparison.summary.txt"
##
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf -w $results/p5_Haplotypecaller/vars_default_realigned.vcf -o default_realigned_neat $results/tmp --incl-homs --incl-fail  -t $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed -T 90 --vcf-out  --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (Realignement default)" >> "../variants_comparison.summary.txt"

fi
##
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf -w $results/p5_Haplotypecaller/vars_default_markDuplicates.vcf -o default_markDuplicates_neat $results/tmp --incl-homs --incl-fail  -t $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed -T 90 --vcf-out  --no-plot 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (markDuplicates default)" >> "../variants_comparison.summary.txt"

####################################################################### generated vcfs
# 1. $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam
#$results/p5_Haplotypecaller/vars_default_bqsr.vcf

# 2. $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.5.bam
#$results/p5_Haplotypecaller/vars_ics5_bqsr.vcf

# 3. $results/p4_indelRealign_bqsr/bqsr/default/recal.ics.7.bam
#$results/p5_Haplotypecaller/vars_ics7_bqsr.vcf

# 4. $results/p4_indelRealign_bqsr/realignedBam.default.0.bam
#$results/p5_Haplotypecaller/vars_default_realigned.vcf

# 5. $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam
#$results/p5_Haplotypecaller/vars_default_markDuplicates.vcf


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
