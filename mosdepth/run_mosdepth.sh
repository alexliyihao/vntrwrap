## Section 1: Header ----------------------------------------------------
#!/bin/sh

# Specify the virtual memory usage
#$ -l h_vmem=15G
#$ -l csg_pri=TRUE
# Specify name to be used to identify this run
#$ -N mosdepth_WHICAP_GRCh37

# Specify the output error file
#$ -o /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SGE_logs/mosdepth/$JOB_NAME_$TASK_ID.outerr
#$ -e /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/SGE_logs/mosdepth/$JOB_NAME_$TASK_ID.outerr
#$ -j y

# This sets the task range in the array (3915 subjects /25 subjects per batch)
#$ -t 1-157:1

# Change directory to the current
#$ -cwd

# Specify source your bash profile. Helpful if you have environment variables set up in your .bash_profile
#$ -S /bin/bash

## Section 2: Path Settings & module loading ----------------------------------------------------

# name of bam files, format provided in the count_read/source_file
SOURCE_LIST="/mnt/mfs/hgrcgrid/shared/LPA_analysis/coassin_pipeline/data_inflow/bam_list_VNTR_N=3915.txt"
# the parent path saving all the bam files
SOURCE_BAM_PATH="/mnt/mfs/hgrcgrid/data/whicap/WHICAP_WES/BAM/washeiDragenBamsList/washeiBamsUpdate2/BQSR/bqsrRGbam"
# the path saving your mosdepth
EXTERNAL_SOFTWARE_PATH="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software"
# the path saving the reference genome, check line 72 for some modification
EXTERNAL_DATA_PATH="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/dataset"
# output path
OUTPUT_PATH="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth"
PROJECT_NAME="covid_GRCh37"
## Section 3: Running ----------------------------------------------------

start_time=$(date +%s)

# read the whole list
readarray -t INPUTFILES < $SOURCE_LIST
echo "Setting: ${#INPUTFILES[@]} files detected"

BATCH_SIZE=25
OFFSET=$(((SGE_TASK_ID - 1) * BATCH_SIZE))
echo "Setting: working on Batch started from $OFFSET, batch size is $BATCH_SIZE"

AGGREGATE_OUTPUT=${OUTPUT_PATH}/mosdepth_${PROJECT_NAME}_batch_${SGE_TASK_ID}.txt
echo "Setting: the output file will be passed to $AGGREGATE_OUTPUT"

for bam_file in "${INPUTFILES[@]:OFFSET:BATCH_SIZE}"
do
    # get bam path
    BAM_PATH=$SOURCE_BAM_PATH/$bam_file
    # get corresponding bai path
    BAI_PATH=${BAM_PATH%?}i

    bam_basename="$(basename $bam_file)"
    # get the save path
    mkdir $OUTPUT_PATH/${bam_basename}_${PROJECT_NAME}

    if [ -f $BAM_PATH -a -f $BAI_PATH ]; then
      	echo "Running: running mosdepth on $bam_file"
        #-n: dont output per-base depth. skipping this output will speed execution
        #--fast-mode # dont look at internal cigar operations or correct mate overlaps
        #--by window size # window size
      	$EXTERNAL_SOFTWARE_PATH/mosdepth_v0.2.5 \
      	    -n \
          	--fast-mode \
      	    --by 1000 \
      	    -f $EXTERNAL_DATA_PATH/GRCh37_reference_genome/human_g1k_v37.fasta \
                  $OUTPUT_PATH/${bam_basename}_${PROJECT_NAME}/$bam_basename \
                  $BAM_PATH

      	echo "Running: Adding $bam_file mosdepth result to $AGGREGATE_OUTPUT"
      	echo -n $bam_file >> $AGGREGATE_OUTPUT
        # Only used chr 6 here, the source file is
        # awk 'length($1)<=2 && $1>=1 && $1<=22 { printf("\t%d",int(100*$4+0.5)) } END {printf("\n")}'
      	zcat $OUTPUT_PATH/${bam_basename}_${PROJECT_NAME}/${bam_basename}.regions.bed.gz \
      	    | awk '$1==6 { printf("\t%d",int(100*$4+0.5)) } END {printf("\n")}' \
      	    >> $AGGREGATE_OUTPUT
    else
      	echo "ERROR: File(s) not found for $bam_file"
    fi
done

gzip -f $AGGREGATE_OUTPUT

echo "Successfully completed mosdepth batch $SGE_TASK_ID, time taken: $(($(date +%s)-$start_time))"
