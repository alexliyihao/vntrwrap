## Section 1: Header ----------------------------------------------------
#!/bin/sh

# Specify the virtual memory usage
#$ -l h_vmem=2G
#$ -l csg_pri=TRUE
# Specify name to be used to identify this run
#$ -N VNTR_count_chr_6_GRCh37

# Specify the output error file
#$ -o /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SGE_logs/VNTR_count/VNTR_count_chr6_GRCh37.outerr
#$ -e /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SGE_logs/VNTR_count/VNTR_count_chr6_GRCh37.outerr
#$ -j y

# Change directory to the current
#$ -cwd

# Specify source your bash profile. Helpful if you have environment variables set up in your .bash_profile
#$ -S /bin/bash

## Section 2: Path Settings & module loading ----------------------------------------------------
# modify this if you are not using module
module load SAMTOOLS/1.17
# IBD2R file, format provided in the source_file
IBD2R_path="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/sources/WES_read_depth/GRCh37_chr6_possible_coding_vntr_regions.IBD2R_gt_0.25.txt"
# name of bam files, format provided in the source_file
SOURCE_LIST="/mnt/mfs/hgrcgrid/shared/LPA_analysis/coassin_pipeline/data_inflow/bam_list_VNTR_N=3915.txt"
# the parent path saving all the bam files
SOURCE_BAM_PATH="/mnt/mfs/hgrcgrid/data/whicap/WHICAP_WES/BAM/washeiDragenBamsList/washeiBamsUpdate2/BQSR/bqsrRGbam"
# output path
OUT_READS="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/VNTR_counts/VNTR_counts_chr6_GRCh37_second_run.txt"

start_time=$(date +%s)

# read the whole list
readarray -t INPUTFILES < $SOURCE_LIST
echo "Setting: ${#INPUTFILES[@]} files detected"

for bam_file in "${INPUTFILES[@]}"
do
    # get bam path
    BAM_PATH=$SOURCE_BAM_PATH/$bam_file
    # get the corresponding bai path
    BAI_PATH=${BAM_PATH%?}i

    bam_basename="$(basename $bam_file)"

    if [ -f $BAM_PATH -a -f $BAI_PATH ]; then
	echo $bam_basename # for logging

	echo -n $bam_file >> $OUT_READS

	while read LINE
	do
	    ARR=( $LINE )
	    CHR=${ARR[0]}
	    BP_START=${ARR[1]}
	    BP_END=${ARR[2]}
	    # echo $CHR:$BP_START-$BP_END
	    samtools view $BAM_PATH $CHR:$BP_START-$BP_END --verbosity 0 \
		| awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {printf("\t%d",ctr)}' \
		>> $OUT_READS
	done < <(awk -v IBD2R_path=$(IBD2R_path) 'NR>1 {print $1,$4,$5}' IBD2R_path)
	echo >> $OUT_READS
    fi
done

echo "Successfully completed VNTR counting, time taken: $(($(date +%s)-$start_time))"
