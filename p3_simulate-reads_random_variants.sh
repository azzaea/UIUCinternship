#!/bin/bash
set -x

# This file takes care of all the steps necessary to generate a synthetic sample

#inputs:
DatasetName=July17_LudasDebuggingTest

## set some paths
# Tools:
NEAT2_path="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/builds/NEAT/neat-genreads"
samtools_path="/home/apps/samtools/samtools-1.3.1/bin"
reference="/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/HG19_GATKbundle2.8_noDecoys.fa"

# Simulation-specific
targeted_region=/home/groups/hpcbio_shared/azza/TargetedRegions-Azza-has-permission/ZachUses_Illumina_truseq_exome_targeted_regions.hg19.chr.bed
known_variants=/home/groups/hpcbio_shared/azza/TargetedRegions-Azza-has-permission/Myrare_variants.vcf

OutputDir=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/$DatasetName
mkdir $OutputDir

NumberOfChunks=5

NumberOfSamples=20

#echo "rng seed value in the files are": > $OutputDir/simulation_parameters.txt
echo "sample-name rng-seed" > $OutputDir/simulation_parameters.txt

echo "NEAT parameters are:" > $OutputDir/neat_parameters.txt

for i in $(seq 1 ${NumberOfSamples})
do
        # Set random error rate each time you run the program. However, we need this to be predictable:
        random_number=$(awk -v i=$i 'BEGIN{srand(i); print int(rand()*i^6+1)}')
	echo ${i} ${random_number} >> $OutputDir/simulation_parameters.txt
	NeatParams="-R 100 -c 10 -E 0.008 -p 2 -t $targeted_region --pe 300 30 --bam --vcf --rng ${random_number} "

        if [ $i -eq 1 ]; then
                NeatParams="-R 100 -c 10 -E 0.008 -p 2 -t $targeted_region -v $known_variants --pe 300 30 --bam --vcf --rng ${random_number} "
        fi

	echo ${NeatParams} >> $OutputDir/neat_parameters.txt

	#create output folders
	OutputFolder=$OutputDir/${i}
	OutputLogsFolder=${OutputFolder}/Logs
	if [ ! -d ${OutputLogsFolder} ]
	then
	   mkdir -p ${OutputLogsFolder}
	fi
	OutputChunksFolder=${OutputFolder}/Chunks
	if [ ! -d ${OutputChunksFolder} ]
	then
	   mkdir -p ${OutputChunksFolder}
	fi
	OutputMergedFolder=${OutputFolder}/Merged
	if [ ! -d ${OutputMergedFolder} ]
	then
	   mkdir -p ${OutputMergedFolder}
	fi
	 
	echo "#!/bin/bash" > ${OutputLogsFolder}/GenericQsub.header
	echo "#PBS -l nodes=1:ppn=5" >> ${OutputLogsFolder}/GenericQsub.header
	echo "#PBS -l walltime=24:00:00" >> ${OutputLogsFolder}/GenericQsub.header
	echo "#PBS -M aeahmed@illinois.edu" >> ${OutputLogsFolder}/GenericQsub.header
	echo "#PBS -m ae" >> ${OutputLogsFolder}/GenericQsub.header
	 
	#### generate and submit the chunks jobs
	####
	AllChunksJobIds=""
	for (( jobnumber=1; jobnumber<=${NumberOfChunks}; jobnumber++ ))
	do
	   cat ${OutputLogsFolder}/GenericQsub.header > ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub
	   echo "#PBS -N ${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub
	   echo "#PBS -e ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.er" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub
	   echo "#PBS -o ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.ou" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub

	   echo " " >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}
		
	   echo "module load python/2.7.9" >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub	  
 
	   echo "python -u ${NEAT2_path}/genReads.py -r ${reference} -o ${OutputChunksFolder}/${DatasetName} ${NeatParams} --job ${jobnumber} ${NumberOfChunks} >  ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.log " >> ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub
 
	   # ChunkJobId is a temporary variable, used to fill out the total list of all job ids
	   ChunkJobId=`qsub ${OutputLogsFolder}/${DatasetName}.GenerateChunks.job_${jobnumber}_of_${NumberOfChunks}.qsub`
           ChunkJobId=${ChunkJobId%%.*}
	   AllChunksJobIds=${AllChunksJobIds}":"${ChunkJobId}
	done
 
 
	#### generate and submit the merge job
	####
	cat ${OutputLogsFolder}/GenericQsub.header > ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
	echo "#PBS -N ${DatasetName}.MergeChunks" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
	echo "#PBS -e ${OutputLogsFolder}/${DatasetName}.MergeChunks.er" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
	echo "#PBS -o ${OutputLogsFolder}/${DatasetName}.MergeChunks.ou" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
	echo "#PBS -W depend=afterok${AllChunksJobIds}" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
 
	echo " " >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub

	echo "module load python/2.7.9" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
 
	echo "python ${NEAT2_path}/mergeJobs.py -i ${OutputChunksFolder}/${DatasetName} -o ${OutputMergedFolder}/${DatasetName}_merged -s  ${samtools_path}/samtools" >> ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub
 
	qsub ${OutputLogsFolder}/${DatasetName}.MergeChunks.qsub


done

