#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=mosdepth_covid_g1k

# This sets the task range in the array (ceil(N/25))
#SBATCH --array 1-82:1

# Memory requirement
#SBATCH --mem=20G

# Change directory to the current
#SBATCH --chdir=./

# Specify that bash shell should be used to process this script
#SBATCH #!/bin/bash

# Specify the outerr file
#SBATCH --output=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/SLURM_logs/mosdepth/%x/%j.out
#SBATCH --error=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/SLURM_logs/mosdepth/%x/%j.err

# Partition, use yours
#SBATCH --partition=GEN

## Section 2: Path Settings & module loading ----------------------------------------------------

# name of bam files, format provided in the count_read/source_file
SOURCE_LIST="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/covid_dataset/data_inflow/complete_files_first_25.txt"
# the parent path saving all the bam files
SOURCE_BAM_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/covid_dataset/covid_bam"
# the path saving your mosdepth
EXTERNAL_SOFTWARE_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/external_software"
# the path saving the reference genome
EXTERNAL_DATA_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset"
REFERENCE_PATH="${EXTERNAL_DATA_PATH}/hs37d5/hs37d5.fa"
# output path
OUTPUT_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth"
# prefix for output
PROJECT_NAME="covid_g1k"

## Section 3: Running ----------------------------------------------------

start_time=$(date +%s)

# read the whole list
readarray -t INPUTFILES < $SOURCE_LIST
echo "Setting: ${#INPUTFILES[@]} files detected"

BATCH_SIZE=25
OFFSET=$(((SGE_TASK_ID - 1) * BATCH_SIZE))
echo "Setting: working on Batch started from $OFFSET, batch size is $BATCH_SIZE"

AGGREGATE_OUTPUT=${OUTPUT_PATH}/${PROJECT_NAME}/mosdepth_${PROJECT_NAME}_batch_${SLURM_ARRAY_TASK_ID}.txt
echo "Setting: the output file will be passed to $AGGREGATE_OUTPUT"

for bam_file in "${INPUTFILES[@]:OFFSET:BATCH_SIZE}"
do
    # get bam path
    BAM_PATH=$SOURCE_BAM_PATH/$bam_file
    # get the corresponding bai path
    # modify this if your file is not one .bam.bai but .bai(check the SGE version)
    BAI_PATH="${BAM_PATH}.bai"

    bam_basename="$(basename $bam_file)"
    # get the save path
    mkdir $OUTPUT_PATH/${PROJECT_NAME}/${bam_basename}_${PROJECT_NAME}

    if [ -f $BAM_PATH -a -f $BAI_PATH ]; then
      	echo "Running: running mosdepth on $bam_file"
      	$EXTERNAL_SOFTWARE_PATH/mosdepth_v0.2.5 \
      	    -n \
          	--fast-mode \
      	    --by 1000 \
      	    -f $REFERENCE_PATH \
                  $OUTPUT_PATH/${PROJECT_NAME}/${bam_basename}_${PROJECT_NAME}/$bam_basename \
                  $BAM_PATH

      	echo "Running: Adding $bam_file mosdepth result to $AGGREGATE_OUTPUT"
      	echo -n $bam_file >> $AGGREGATE_OUTPUT
        # Only used chr 6 here, the source file is
        # awk 'length($1)<=2 && $1>=1 && $1<=22 { printf("\t%d",int(100*$4+0.5)) } END {printf("\n")}'
      	zcat $OUTPUT_PATH/${PROJECT_NAME}/${bam_basename}_${PROJECT_NAME}/${bam_basename}.regions.bed.gz \
      	    | awk '$1==6 { printf("\t%d",int(100*$4+0.5)) } END {printf("\n")}' \
      	    >> $AGGREGATE_OUTPUT

          else
      	echo "ERROR: File(s) not found for $bam_file"
    fi
done

gzip -f $AGGREGATE_OUTPUT

echo "Successfully completed mosdepth batch $SLURM_ARRAY_TASK_ID, time taken: $(($(date +%s)-$start_time))"
