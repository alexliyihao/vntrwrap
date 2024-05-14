#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=covid_VNTR_count_chr_6_GRCh37

# Memory requirement
#SBATCH --mem=2G

# Change directory to the current
#SBATCH --chdir=./

# Specify that bash shell should be used to process this script
#SBATCH #!/bin/bash

# Specify the output error file
#SBATCH --output=/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SLURM_logs/VNTR_count/%x.out
#SBATCH --error=/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SLURM_logs/VNTR_count/%x.err

# Partition, use yours
#SBATCH --partition=GEN

## Section 2: Path Settings & module loading ----------------------------------------------------
# modify this if you are not using module
module load SAMTOOLS/1.17
# IBD2R file, format provided in the source_file
IBD2R_path="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/sources/WES_read_depth/GRCh37_chr6_possible_coding_vntr_regions.IBD2R_gt_0.25.txt"
# name of bam files, format provided in the source_file
SOURCE_LIST="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/covid_dataset/data_inflow/index_list.txt"
# the parent path saving all the bam files
SOURCE_BAM_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/covid_dataset/covid_bam"
# output path
OUT_READS="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/analysis_results/VNTR_counts/VNTR_counts_chr6_GRCh37_covid_dataset.txt"

start_time=$(date +%s)

# read the whole list
readarray -t INPUTFILES < $SOURCE_LIST
echo "Setting: ${#INPUTFILES[@]} files detected"

for bam_file in "${INPUTFILES[@]}"
do
    # get bam path
    BAM_PATH=$SOURCE_BAM_PATH/$bam_file
    # get the corresponding bai path
    # modify this if your file is not one .bam.bai but .bai(check the SGE version)
    BAI_PATH="$BAM_PATH.bai"
    bam_basename="$(basename $bam_file)"

    if [ -f $BAM_PATH -a -f $BAI_PATH ]; then
	echo $bam_basename # for logging
	echo -n $bam_file >> $OUT_READS
	while read LINE
	do
	    ARR=( $LINE )
	    CHR=${ARR[0]}
	    BP_START_HG19=${ARR[1]}
	    BP_END_HG19=${ARR[2]}
	    # echo "${CHR}:${BP_START_HG19}-${BP_END_HG19}"
	    samtools view $BAM_PATH $CHR:$BP_START_HG19-$BP_END_HG19 --verbosity 0 \
		| awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {printf("\t%d",ctr)}' \
		>> $OUT_READS
	done < <(awk 'NR>1 {print $1,$4,$5}' $IBD2R_path)
	echo >> $OUT_READS
    fi
done

echo "Successfully completed VNTR counting, time taken: $(($(date +%s)-$start_time))"
