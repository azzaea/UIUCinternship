#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-VariantEval.log.txt
#PBS -e localhost:$HOME/errors-neat_VariantEval.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/reads"
results="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results"
TargetedRegions="/home/groups/hpcbio_shared/azza/TargetedRegions-Azza-has-permission/"

#### Needed modules:

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

module load R/3.2.3


module load hap.py/0.3.0

# HC was called for 5 different files, and it would be nice to compare the variants called in each case. These files are:
# 1. results/p5_Haplotypecaller/vars_default_bqsr.vcf 
# 2. $results/p5_Haplotypecaller/vars_ics5_bqsr.vcf
# 3. $results/p5_Haplotypecaller/vars_ics7_bqsr.vcf
# 4. $results/p5_Haplotypecaller/vars_default_realigned.vcf 
# 5. $results/p5_Haplotypecaller/vars_default_markDuplicates.vcf

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

if [ ! -d $results/p5_Haplotypecaller/comparison_stages/hap.py ]; then
        mkdir -p $results/p5_Haplotypecaller/comparison_stages/hap.py
fi

cd $results/p5_Haplotypecaller/comparison_stages/hap.py


##
#2. Using VariantEval from GATK

START=$(date +%s)
hap.py \
 $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf \
 $results/p5_Haplotypecaller/vars_default_bqsr.vcf \
 -f  $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed \
 -o default_bqsr

END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR default)" >> "../variants_comparison.summary.txt"

###
START=$(date +%s)
hap.py \
 $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf \
 $results/p5_Haplotypecaller/vars_ics5_bqsr.vcf \
 -f $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed \
 -o ics5_bqsr
END=$(date +%s)
IFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR ics5)" >> "../variants_comparison.summary.txt"
##
START=$(date +%s)
hap.py \
 $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf \
 $results/p5_Haplotypecaller/vars_ics7_bqsr.vcf \
 -o ics7_bqsr \
 -f $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed 
 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (BQSR ics7)" >> "../variants_comparison.summary.txt"
##
START=$(date +%s)
hap.py \
  $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf \
 $results/p5_Haplotypecaller/vars_default_realigned.vcf \
 -o default_realigned_ \
 -f $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed 

END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (Realignement default)" >> "../variants_comparison.summary.txt"
##
START=$(date +%s)
hap.py \
 $reads/H3A_NextGen_assessment.Chr1_50X.set3_golden.vcf \
 $results/p5_Haplotypecaller/vars_default_markDuplicates.vcf \
 -f $TargetedRegions/Intersect.TruseqTargeted_with_PlatinumConfident.chr1.bed \
 -o default_markduplicates
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (markDuplicates default)" >> "../variants_comparison.summary.txt"

echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
