#!/bin/bash

# This script gets the names of the files in the specified parent directory, finds the sample names, and then for each, appends the sample name to names of the files in the merged directory.

Datasetname=synthetic_1_to_20_rare_to_common_individuals_ration


parent_dir=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/$Datasetname

cd $parent_dir

dir_name=${PWD##*/}

mkdir generated_reads_$dir_name

dirs=(*/)

for s in $(seq 1 ${#dirs[@]}); do 
	cd $s
	sample=${PWD##*/}
	cd $parent_dir/${s}/Chunks
	rename synthetic synthetic.sample$sample *
	# also fix this: for ease later on, make the name:
	# synthetic.sample<#>_var_$dir_name_merged_read<#>.fq
	# use something lik:
	#       rename _var_err_rate _var_$dir_name *
	cp $parent_dir/${s}/Chunks/* $parent_dir/generated_reads_$dir_name
	cd $parent_dir
done

